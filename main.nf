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

process check_no_demux {
	echo true

	tag "${fcid}"

	beforeScript "export PYTHONPATH=$PYTHONPATH:${workflow.projectDir}/bin"

	input:
	val lane

	output:
	tuple val(lane), env(no_demux), emit: lane

	shell:
	"""
	no_demux=\$(python3 -c "from slime import check_demux;r=check_demux('${fcid}', ${lane});print(str(r).lower())")
	echo "check_no_demux: lane: $lane, no_demux: \$no_demux"
	"""
}

process check_do_merge {
	echo true

	tag "${fcid}"

	beforeScript "export PYTHONPATH=$PYTHONPATH:${workflow.projectDir}/bin"

	output:
	env(do_merge), emit: do_merge

	shell:
	"""
	do_merge=\$(python3 -c "from slime import check_do_merge;r=check_do_merge('${fcid}');print(str(r).lower())")
	echo "check_do_merge: do_merge: \$do_merge"
	"""
}


process tar {
	publishDir "${alpha}/logs/${fcid}/archive/tar/", mode:'copy', failOnError: true, pattern: '.command.*'

	tag "${fcid}"

	output:
	path("*.tar.bz2"), emit: file
	path(".command.*")

	script:
	"""
	# need the h flag to follow symlinks (?)
	# tar hcvjf ${run_dir_name}.tar.bz2 ${run_dir_path}
	module load pbzip2/1.1.13
	tar -c ${run_dir_path} | pbzip2 -c -p${task.cpus} -m2000 > ${run_dir_name}.tar.bz2 
	"""
}

process rsyncToArchive {
	echo true 

	publishDir "${alpha}/logs/${fcid}/archive/rsync/", mode:'copy', failOnError: true
	
	tag "${fcid}"

	input:
	path SRC
	val destination

	output:
	path(".command.*")
	env(exit_code), emit: exit_code // used as dependency for QC

	shell:
	"""
	echo "rsyncToArchive: SRC: ${SRC}, destination: ${destination}"
	ssh -i \$HOME/.ssh/id_rsa core2 "mkdir -p ${destination}"
	rsync --copy-links --progress -r -e \"ssh -i \${HOME}/.ssh/id_rsa\" ${SRC} core2:${destination}/.
	exit_code=\$?
	"""
}

process _basecall {
	publishDir "${alpha}/logs/${fcid}/basecall/${lane}", mode:'copy', failOnError: true, pattern: '.command.*'

	tag "${fcid}"

	input:
	val lane

	output:
	val(lane), emit: lane
	path(".command.*")

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

	out_path="${alpha}/lane/${fcid}/${lane}"
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
	echo true
	
	publishDir "${alpha}/pheniqs_conf/${fcid}/${lane}", mode:'copy', pattern: 'demux.json'
	publishDir "${alpha}/logs/${fcid}/demux/make_config/${lane}", mode:'copy', failOnError: true

	tag "${fcid}"

	input:
	tuple val(lane), val(no_demux)

	output:
	tuple val(lane), path('demux.json'), emit: pheniqs_conf
	path(".command.*")

	when:
	no_demux == "false"

	script:
	"""
	pheniqs_config.py \
		$fcid \
		$lane \
		${task.cpus}
	"""
}

process run_pheniqs {
	echo true

	publishDir "${alpha}/logs/${fcid}/demux/pheniqs/${lane}", mode:'copy', failOnError: true
	publishDir "${alpha}/pheniqs_out/${fcid}/${lane}", mode:'copy', failOnError: true, pattern: 'demux.out'

	tag "${fcid}"

	input:
	tuple(val(lane), file(pheniqs_conf))

	output:
	val(lane), emit: lane	// used as trigger for qc
	path("demux.out")
	path(".command.*")

	shell
	"""
	module load $PHENIQS_MODULE
	rm -rf ${alpha}/sample/${fcid}/${lane}/*
	pheniqs demux -C $pheniqs_conf > 'demux.out'
	"""
}

process demux_reports {
    echo true

    publishDir "${alpha}/logs/${fcid}/qc/${workflow.runName}/demux_reports/", mode:'copy', failOnError: true

    tag "${fcid}"

    module "anaconda3/2020.07"
    conda "/scratch/cgsb/gencore/mk5636/conda/slime"

    input:
    val lanes //because data might be merged, need to wait for all lanes

    output:
    path('*/*_summary_mqc.txt'), emit: summary_report
    path('*/*_mqc.txt'), emit: reports
    path(".command.*")

    shell:
    """
    demux_report.py $fcid
    """
}

process merge_lanes {
    publishDir "${alpha}/logs/${fcid}/qc/${workflow.runName}/merge/", mode:'copy', failOnError: true

    tag "${fcid}"

    beforeScript "export PYTHONPATH=$PYTHONPATH:${workflow.projectDir}/bin"

    input:
    val lanes //need to wait for all lanes
    val do_merge

    output:
    env(exit_code), emit: exit_code
    path(".command.*")

    when:
    do_merge == "true"

    shell:
    """
    merge.sh $fcid $alpha
    exit_code=\$?
    """
}

process get_lane_paths {
	echo true

	publishDir "${alpha}/logs/${fcid}/qc/${workflow.runName}/get_lane_paths/", mode:'copy', failOnError: true

	tag "${fcid}"

	beforeScript "export PYTHONPATH=$PYTHONPATH:${workflow.projectDir}/bin"

	input:
	val do_merge

	output:
	path(".command.*")
	path("lanes.txt"), emit: paths

	shell:
	"""
	if [[ "$do_merge" == "true" ]]; then
		echo "${alpha}/merged/${fcid}/merged" > lanes.txt
	else
		# need to set path to lane to get all lanes
		# then check each lane individually to see
		# if it has been demuxed
		path="${alpha}/lane/${fcid}"
		
		for dir in \${path}/*
		do
			dir=\$(basename "\$dir")

			# Call Python function and get the value
			no_demux=\$(python3 -c "from slime import check_demux;r=check_demux('${fcid}', \${dir});print(str(r).lower())" 2>&1)
			
			# Check if Python script executed successfully
			if [ \$? -ne 0 ]; then
				echo "ERROR: Failed to assign no_demux: \$no_demux; lane: \$dir"
				exit 1
			fi

			#echo "lane: \$dir, no_demux: \$no_demux"

			# Convert the boolean result to 'lane' or 'sample'
			if [[ "\$no_demux" == "true" ]]; then
				value="lane"
			else
				value="sample"
			fi

			echo "lane \$dir path:"
			#echo "${alpha}/\${value}/${fcid}/\${dir}" >> lanes.txt
			echo ""${alpha}/\${value}/${fcid}/\${dir}"" | tee -a lanes.txt
		done
	fi
	"""
}

process fastqc {
    echo true

    publishDir "${alpha}/logs/${fcid}/qc/${workflow.runName}/fastqc/${lane}", mode:'copy', failOnError: true

    tag "${fcid}"

    beforeScript "export PYTHONPATH=$PYTHONPATH:${workflow.projectDir}/bin"

    input:
    val path
    val merge_exit_code
    val qc_dep // need to wait for basecalling/demuxing to finish

    output:
    path(".command.*")
    tuple(val(path), path("*"), emit: delivery_path_and_fastqc_files)
 
    when:
    merge_exit_code == "0"

    script:
    lane = new File(path).getName().trim()
    """
    echo "path: $path"
    mkdir $lane
    cd $lane
    do_fastqc.sh $path
    cd ..
    """
}

process multiqc {
    echo true

    publishDir "${alpha}/logs/${fcid}/qc/${workflow.runName}/multiqc/${lane}", mode:'copy', failOnError: true

    tag "${fcid}"

    module "multiqc/1.9"

    input:
    tuple(val(path_to_data), path(fastqc))
    path reports

    output:
    path("${lane}/"), emit: files
    tuple(val(path_to_data), path("${lane}/${fcid}_${lane}_summary_mqc.txt"), emit: delivery_path_and_summary_report)
    path(".command.*")

    shell:
    lane = fastqc.getName()
    """
    echo "fastqc: $fastqc"
    echo "reports: $reports"
    echo "lane: $lane"
    echo "path: $path_to_data"

    cp ${fcid}_${lane}_*_mqc.txt ${lane}/.
    cd ${lane}
    multiqc -f -c ${workflow.projectDir}/bin/mqc_config.yaml .
    """
}

process deliver {
    echo true

    publishDir "${alpha}/logs/${fcid}/qc/${workflow.runName}/deliver/${lane}", mode:'copy', failOnError: true

    tag "${fcid}"

    input:
    val archive_exit_code
    val rsync_exit_code
    tuple(val(path), path(summary_report))

    output:
    path(".command.*")

    when:
    archive_exit_code.toInteger() == 0 && rsync_exit_code.toInteger() == 0

    shell:
    lane = new File(path).getName().trim()
    """
    echo "path: $path"
    echo "summary_report: $summary_report"
    echo "archive_exit_code: $archive_exit_code"

    # check qc and deliver
    qc_deliver.py ${path.trim()} ${summary_report}
    """
}

workflow _demux{
    take:
	lane
    main:
	check_no_demux(lane)
	make_pheniqs_config(check_no_demux.out.lane)
	run_pheniqs(make_pheniqs_config.out.pheniqs_conf)
    emit:
	run_pheniqs.out.lane
}

workflow _qc{
    take:
	archive_exit_code
        qc_dep
    main:
        demux_reports(qc_dep)
	check_do_merge()
        do_merge = check_do_merge.out.do_merge
        merge_lanes(qc_dep, do_merge)
	get_lane_paths(do_merge)
	path = get_lane_paths.out.paths.splitText()
        fastqc_dep = merge_lanes.out.exit_code.ifEmpty("0")
        fastqc(path, fastqc_dep, qc_dep)
        multiqc(fastqc.out.delivery_path_and_fastqc_files, demux_reports.out.reports)
        rsyncToArchive(multiqc.out.files, fastqc_path + "/${fcid}/")
	deliver(archive_exit_code, rsyncToArchive.out.exit_code, multiqc.out.delivery_path_and_summary_report)
}

workflow archive{
    main:
	tar()
	rsyncToArchive(tar.out.file, archive_path + "/${seq_id}/")
    emit:
	rsyncToArchive.out.exit_code
}

workflow basecall{
    _basecall(lanes)
    _demux(_basecall.out.lane) 
    qc_dep = _demux.out.collect().ifEmpty(_basecall.out.lane.collect())
    _qc(0, qc_dep)
}

workflow demux{
    _demux(lanes) 
    qc_dep = _demux.out.collect().ifEmpty('')
    _qc(0, qc_dep)
}

workflow qc{
    _qc(0, '') 
}

workflow pheniqs_conf{
    check_no_demux(lanes)
    make_pheniqs_config(check_no_demux.out.lane)
}

workflow {
    archive()
    _basecall(lanes)
    _demux(_basecall.out.lane) 
    qc_dep = _demux.out.collect().ifEmpty(_basecall.out.lane.collect())
    _qc(archive.out, qc_dep)
}

workflow.onComplete {
	def status = "NA"
	if(workflow.success) {
		status = "SUCCESS"
	} else {
		status = "FAILED"
	}

	sendMail {
		to ""
		subject "${fcid} ${status}"

		"""
		Pipeline execution summary
		---------------------------
		Success           : ${workflow.success}
		exit status       : ${workflow.exitStatus}
		Launch time       : ${workflow.start.format('dd-MMM-yyyy HH:mm:ss')}
		Ending time       : ${workflow.complete.format('dd-MMM-yyyy HH:mm:ss')} (duration: ${workflow.duration})
		Launch directory  : ${workflow.launchDir}
		Work directory    : ${workflow.workDir.toUriString()}
		Project directory : ${workflow.projectDir}
		Run directory     : ${params.run_dir_path}
		Script ID         : ${workflow.scriptId ?: '-'}
		Workflow session  : ${workflow.sessionId}
		Nextflow run name : ${workflow.runName}
		Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
		
		Command:
		${workflow.commandLine}

		Errors:
		Error Message:    : ${workflow.errorMessage}
		Error Report      : ${workflow.errorReport}
		"""
    }
}

workflow.onError {
	sendMail {
		to ""
		subject "${fcid} Error: ${workflow.errorMessage}"

		"""
		Pipeline execution summary
		---------------------------
		Error Message:    : ${workflow.errorMessage}
		Error Report      : ${workflow.errorReport}

		Command:
		${workflow.commandLine}
		"""
    }
}
