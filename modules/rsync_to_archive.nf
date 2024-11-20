process RSYNC_TO_ARCHIVE {
    publishDir "${params.alpha}/logs/${params.fcid}/archive/rsync/", mode:'copy', failOnError: true
    tag "${params.fcid}"
    echo true

    input:
    path SRC
    val destination

    output:
    path(".command.*")
    env(exit_code), emit: exit_code

    shell:
    '''
    echo "rsyncToArchive: SRC: !{SRC}, destination: !{destination}"
    ssh -i $HOME/.ssh/id_rsa core2 "mkdir -p !{destination}"
    rsync --copy-links --progress -r -e "ssh -i ${HOME}/.ssh/id_rsa" !{SRC} core2:!{destination}/.
    exit_code=$?
    '''
}