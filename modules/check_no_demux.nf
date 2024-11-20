process CHECK_NO_DEMUX {
    tag "${params.fcid}"
    echo true

    input:
    val lane

    output:
    tuple val(lane), env(no_demux), emit: lane

    shell:
    '''
    export PYTHONPATH=$PYTHONPATH:!{workflow.projectDir}/bin
    no_demux=$(python3 -c "from slime import check_demux;r=check_demux('!{params.fcid}', !{lane});print(str(r).lower())")
    echo "check_no_demux: lane: !{lane}, no_demux: $no_demux"
    '''
}