nextflow.enable.dsl=2

def lanes = channel.from(1..params.lanes)
def run_dir_path = params.run_dir_path
def run_dir_name = new File(run_dir_path).getName()
def parts = run_dir_name.split('_')
def seq_id = ""
def fcid = ""

seq_id = parts[1]

def fcidPart = parts[3]
if (fcidPart.startsWith("A") || fcidPart.startsWith("B")) {
	fcid = fcidPart.substring(1)
} else {
	fcid = fcidPart
}

println "seq_id: $seq_id"
println "fcid: $fcid"

process tar {
	afterScript 'echo "Tar completed." | mail -s "Process Complete" your_email@example.com'

	input: 
	val tarName

	output:
	path '*.tar.bz2', emit: file

	script:
	"""
	tar cvjf ${tarName}.tar.bz2 ${params.run_dir_path}
	"""
}

process rsyncToArchive {
	input:
	path SRC 
	val seq_id

	shell:
	"""
	echo "rsync --copy-links --progress -e \"ssh -i ${HOME}/.ssh/id_rsa\" !{SRC} mk5636@babylon.bio.nyu.edu:/mnt/gencore/hpcpipe/archive/!{seq_id}"
	"""
}

process _basecall {
	echo true

	input:
	val lane

	output:
	val(lane), emit: lane

	shell:
	"""
	read_structure=\$(python3 -c "
	import xml.dom.minidom

	read_structure = ''
	runinfo = xml.dom.minidom.parse('${params.run_dir_path}/RunInfo.xml')
	nibbles = runinfo.getElementsByTagName('Read')

	for nib in nibbles:
		read_structure += nib.attributes['NumCycles'].value + 'T'

	print(read_structure)
	")

	run_barcode=\$(python3 -c "
	print('${params.run_dir_path}'.split('_')[-2].lstrip('0'))
	")

	out_path="${alpha}lane/default/${fcid}/${lane}"
	mkdir -p \$out_path

	tmp_work_dir="${tmp_dir}${fcid}/${lane}"
	mkdir -p \$tmp_work_dir

	module load $PICARD_MODULE
	module load $JDK_MODULE

	java -jar -Xmx58g \$PICARD_JAR IlluminaBasecallsToFastq \
		LANE=${lane} \
		READ_STRUCTURE=\${read_structure} \
		BASECALLS_DIR=${params.run_dir_path}/Data/Intensities/BaseCalls \
		OUTPUT_PREFIX=\${out_path}/${fcid}_l0${lane} \
		RUN_BARCODE=\${run_barcode} \
		MACHINE_NAME=${seq_id} \
		FLOWCELL_BARCODE=${fcid} \
		NUM_PROCESSORS=${task.cpus} \
		APPLY_EAMSS_FILTER=false \
		INCLUDE_NON_PF_READS=false \
		MAX_READS_IN_RAM_PER_TILE=200000 \
		MINIMUM_QUALITY=2 \
		COMPRESS_OUTPUTS=true \
		TMP_DIR=\${tmp_work_dir}

	rm -rf \${tmp_work_dir}
	"""
}

process make_pheniqs_config {
	publishDir "${alpha}/${fcid}/pheniqs_conf/${lane}", mode:'copy'

	input:
	val lane

	output:
	tuple(val(lane), path('demux.json'))

	when:
	!params.no_demux

	shell
	"""
	pheniqs_config.py \
		$seq_id \
		$fcid \
		${task.cpus} \
		$lane \
		${params.run_dir_path} \
		$alpha
	"""
}

process run_pheniqs {
	echo true

	input:
	tuple(val(lane), file(pheniqs_conf))

	output:
	tuple(env(lane), path('demux.out'))

	when:
	!params.no_demux

	shell
	"""
	module load $PHENIQS_MODULE
	rm -rf ${alpha}/sample/default/${fcid}/${lane}/*
	pheniqs demux -C $pheniqs_conf > 'demux.out'
	lane=${lane}
	"""

}

process _qc {
	echo true

	input:
	val pheniqs_out

	output:
	//env(success), emit: success

	shell
	"""
	demux_report.py "${pheniqs_out}" $run_dir_path $fcid $params.no_demux
	merge.sh $fcid $params.no_demux $alpha
	do_fastqc.sh </path/to/data: either alpha/sample/default/fcid, alpha/lane/default/fcid/, or alpha/merged/fcid>
	//rsync_qc_to_web 
	"""

}

process _deliver {
	echo true

	input:
	val qc_success

	output:

	when:
	qc_success == 'true'

	shell
	"""
	echo "test delivery, qc_success: ${qc_success}"
	"""

}

process _test {
	echo true

	shell:
	'''
	lane_num=4
	run_dir_name=$(basename !{params.run_dir_path})
	fcid=${run_dir_name##*_}
	if [[ ${fcid} != *"-"* ]] && [ ${#fcid} -gt 5 ];then
		fcid=${fcid:1}
	fi
	URL="!{tw_api_root}flowcell/illumina/${fcid}/?username=!{tw_user}&api_key=!{tw_api_key}"
	o=$(curl -s "$URL")
	echo "o = " $o
	id=$(echo "$response" | python3 -c '
	import json

	response = input()
	data = json.loads(response)
	lane_num = '$lane_num'
	id = [lane["id"] for lane in data["lanes"] if lane["lane_number"] == int(lane_num)][0]
	print(id)
	')
	'''
}

workflow test{
	_test()
}

workflow _demux{
	take:
		lane
	main:
		make_pheniqs_config(lane)
		run_pheniqs(make_pheniqs_config.out)
	emit:
		run_pheniqs.out
}

workflow archive{
	tar(run_dir_name)
	rsyncToArchive(tar.out.file, seq_id)
}

workflow basecall{
	_basecall(lanes)
	_demux(_basecall.out.lane) 
	def qc_dep = params.no_demux ? _basecall.out.lane.collect() : _demux.out.o.view().collect()
	_qc(qc_dep) 
	_deliver(_qc.out.success)
}

workflow demux{
	_demux(lanes) 
	def qc_dep = params.no_demux ? 'x' : _demux.out.o.view().collect()
	_qc(qc_dep) 
	_deliver(_qc.out.success)
}

workflow qc{
	_qc('x') 
	_deliver(_qc.out.success)
}

workflow deliver{
	_deliver('true')
}

workflow pheniqs_conf{
	make_pheniqs_config(lanes)
}

workflow {
	//archive()
	_basecall(lanes)
	_demux(_basecall.out.lane) 
	def qc_dep = params.no_demux ? _basecall.out.lane.collect() : _demux.out.collect()
	_qc(qc_dep) 
	//_deliver(_qc.out.success)
}