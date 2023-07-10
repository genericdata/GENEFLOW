nextflow.enable.dsl=2

def run_dir_path = params.run_dir_path
def run_dir_name = new File(run_dir_path).getName()
def parts = run_dir_name.split('_')
def seq_id = parts[1]
def fcidPart = parts[3]
def fcid = fcidPart.matches("^[AB].*") ? fcidPart.substring(1) : fcidPart

def num_lanes = new File(run_dir_path,'/Data/Intensities/BaseCalls')
	.listFiles()
	.findAll { it.name ==~ /L[0-9]{3}/ }
	.size()

def lanes = channel.from(1..num_lanes)

println "run_dir_path: $run_dir_path"
println "seq_id: $seq_id"
println "fcid: $fcid"
println "num_lanes: $num_lanes"

process tar {
	//afterScript 'echo "Tar completed." | mail -s "Process Complete" your_email@example.com'

	output:
	path '*.tar.bz2', emit: file

	script:
	"""
	tar cvjf ${run_dir_name}.tar.bz2 ${run_dir_path}
	"""
}

process rsyncToArchive {
	input:
	path SRC 

	shell:
	"""
	rsync --copy-links --progress -e \"ssh -i \${HOME}/.ssh/id_rsa\" ${SRC} core2:${archive_path}/${seq_id}/.
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
	env(success), emit: success

	beforeScript "export PYTHONPATH=$PYTHONPATH:${workflow.projectDir}/bin"

	shell
	"""
	# create reports
	demux_report.py "${pheniqs_out}" $run_dir_path $fcid $params.no_demux

	# check if we should merge, and merge if so
	do_merge=\$(python3 -c "from slime import check_do_merge;r=check_do_merge('HHWYGAFX5');print(str(r).lower())")
	echo "do_merge: \$do_merge"
	
	if [ "\$do_merge" == "true" ];then
		merge.sh $fcid $params.no_demux $alpha
	fi

	# get the path to the deliverable fastq files and run fastqc
	if [ "\$do_merge" = "true" ]; then
		path="$alpha/merged/$fcid"
	elif [ "\$do_merge" = "false" ] && [ "$params.no_demux" = "false" ]; then
		path="$alpha/sample/default/$fcid"
	elif [ "\$do_merge" = "false" ] && [ "$params.no_demux" = "true" ]; then
		path="$alpha/lane/default/$fcid"
	else
		echo "Invalid states. Check your do_merge and no_demux variables."
	fi

	echo "Path is set to: \$path"
	do_fastqc.sh \$path ${workflow.projectDir}/bin/mqc_config.yaml

	# rsync to core-fastqc server (mkdir first)
	ssh -i $HOME/.ssh/id_rsa core2 "mkdir -p ${fastqc_path}/${fcid}/"
	rsync -r -e "ssh -i $HOME/.ssh/id_rsa" --exclude=".*" . core2:${fastqc_path}/${fcid}/.

	# check qc and deliver
	summary_report=\$(ls */*summary_mqc.txt)
	qc.py \$summary_report \$path $run_dir_path

	success=false
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
	tar()
	rsyncToArchive(tar.out.file)
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
	archive()
	_basecall(lanes)
	_demux(_basecall.out.lane) 
	def qc_dep = params.no_demux ? _basecall.out.lane.collect() : _demux.out.collect()
	_qc(qc_dep) 
	//_deliver(_qc.out.success)
}
