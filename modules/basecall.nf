process BASECALL {
    publishDir "${params.alpha}/logs/${params.fcid}/basecall/${lane}", mode:'copy', failOnError: true, pattern: '.command.*'
    tag "${params.fcid}"

    module params.PICARD_MODULE
    module params.JDK_MODULE

    input:
    val lane

    output:
    val(lane), emit: lane
    path(".command.*")

    shell:
    '''
    read_structure=$(python3 -c "
    import xml.dom.minidom

    read_structure = ''
    runinfo = xml.dom.minidom.parse('!{params.run_dir_path}/RunInfo.xml')
    nibbles = runinfo.getElementsByTagName('Read')

    for nib in nibbles:
        read_structure += nib.attributes['NumCycles'].value + 'T'
    
    print(read_structure)
    ")

    run_barcode=$(python3 -c "
    print('!{params.run_dir_path}'.split('_')[-2].lstrip('0'))
    ")

    out_path="!{params.alpha}/lane/!{params.fcid}/!{lane}"
    mkdir -p $out_path

    tmp_work_dir="!{params.tmp_dir}!{params.fcid}/!{lane}"
    mkdir -p $tmp_work_dir

    java -jar -Xmx58g $PICARD_JAR IlluminaBasecallsToFastq \
        LANE=!{lane} \
        READ_STRUCTURE=${read_structure} \
        BASECALLS_DIR=!{params.run_dir_path}/Data/Intensities/BaseCalls \
        OUTPUT_PREFIX=${out_path}/!{params.fcid}_l0!{lane} \
        RUN_BARCODE=${run_barcode} \
        MACHINE_NAME=!{params.seq_id} \
        FLOWCELL_BARCODE=!{params.fcid} \
        NUM_PROCESSORS=!{task.cpus} \
        APPLY_EAMSS_FILTER=false \
        INCLUDE_NON_PF_READS=false \
        MAX_READS_IN_RAM_PER_TILE=200000 \
        MINIMUM_QUALITY=2 \
        COMPRESS_OUTPUTS=true \
        TMP_DIR=${tmp_work_dir}

    rm -rf ${tmp_work_dir}
    '''
}