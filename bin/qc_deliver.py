#!/usr/bin/env python3

import sys
import os
import requests
import json
import glob
import shutil
import subprocess
from slime import *
from config import tw_api_root, tw_api_key, tw_user, delivery_folder_root, raw_run_dir_delivery_root, alpha


def parse_summary(summary_report):
    with open(summary_report, 'r') as f:
        lines = f.readlines()
        
    header_line = lines[3].split("\t")
    data_line = lines[4].split("\t")

    stats = {header_line[i].strip(): float(data_line[i].replace(',','')) for i in range(len(header_line))}
    return stats


def check_pool_errors(stats):
    success = True
    error = ''
    if "% Undetermined" in stats and stats["% Undetermined"] > 50:
        error = "undetermined % > 50 error: undetermined % = {}".format(
            stats["% Undetermined"]
        )
        success = False
    return {"success": success, "message": error}


def get_qc_messages(stats):
    message = ""
    if ("% PhiX Aligned" in stats and stats["% PhiX Aligned"] > 15) \
        or stats["Total # PF Reads"] / stats["Total # of Single-End Reads"] < 0.60:
        message += "<br>Unfortunately, your lane(s) received less reads than expected. This is likely related to an issue with metadata, quantification and/or pooling. We would be happy to discuss this further if desired.<br>"

    if (
        "% Undetermined" in stats
        and stats["% Undetermined"] - stats["% PhiX Aligned"] > 15
    ):
        message += "<br>There is a large number of 'undetermined reads' in your lane(s), those whose index cannot be appropriately identified as any of the ones expected for the libraries. This is likely related to an issue with metadata, pooling, and/or library prep. We would be happy to discuss this further if desired. <br>"

    return message


def set_run_status(fcid, status):
    url = f"{tw_api_root}flowcell/illumina/{fcid}/status/{status}"
    params = {
        'username': tw_user,
        'api_key': tw_api_key
    }
    response = requests.get(url, params=params)
    return response.json()


def get_delivery_fcid(fcid):
    return fcid.split('-')[1] if '-' in fcid else fcid


def update_lane_stats(name, lane_num, stat_name, stat_value):
    url = f"{tw_api_root}scanner/update_lane_stats/{name}/{lane_num}/{stat_name}/{stat_value}?username={tw_user}&api_key={tw_api_key}"
    response = requests.get(url)
    data = response.json()
        
    if 'error_message' in data:
        return {'success': False, 'msg': data['error_message']}
    
    return {'success': True, 'msg': ''}


def get_delivery_email(fcid, delivery_dir, run_dir_user_path, mqc_report_url, message, no_demux, allowed_barcode_mismatch, manufacturer):
    demux_text = f"Following basecalling, the reads were demultiplexed using Pheniqs version 2.1.0<sup>2</sup>." if not no_demux else ""
    basecall_text = "Picard IlluminaBasecallsToFastq version 2.23.8<sup>1</sup>, with APPLY_EAMSS_FILTER set to false" if manufacturer == "Illumina" else "Bases2Fastq version 2.1.0<sup>1</sup>"

    basecall_citation = '"Picard Toolkit.” 2019. Broad Institute, GitHub Repository. https://broadinstitute.github.io/picard/; Broad Institute.' if manufacturer == "Illumina" else "Bases2Fastq. Element Biosciences. https://docs.elembio.io/docs/bases2fastq/; Element Biosciences."
    
    delivery_template = f'''
<p>Dear GenCore Users,</p>

<p>Your data has been successfully processed by the <a href="https://github.com/gencorefacility/GENEFLOW">GENEFLOW</a> pipeline.</p>

<p>We are providing a suggested methods section for your publications:</p>

<p>The reads were basecalled using {basecall_text}. {demux_text} The entire process was executed using a custom nextflow pipeline, GENEFLOW<sup>3</sup>.</p>

<p>References</p>

<ol>
    <li>{basecall_citation}</li>
    <li>Galanti, L., Shasha, D. & Gunsalus, K.C. Pheniqs 2.0: accurate, high-performance Bayesian decoding and confidence estimation for combinatorial barcode indexing. BMC Bioinformatics 22, 359 (2021). https://doi.org/10.1186/s12859-021-04267-5.</li>
    <li>"GENEFLOW." 2023. New York University Center for Genomics and System Biology Genomics Core, GitHub Repository. https://github.com/gencorefacility/GENEFLOW; New York University.</li>
</ol>

<p>We kindly request that you include the following acknowledgements in your publications if this data contributed to your research:</p>
<p>This work was supported in part through the NYU IT High Performance Computing resources, services, and staff expertise. We acknowledge the Zegar Family Foundation for their generous support. We thank the NYU Center for Genomics and System Biology Genomics Core for their assistance and resources.</p>

<p>You must have the required permissions to access data on the HPC. If this is your first time sequencing, please visit: <a href="https://gencore.bio.nyu.edu/bioinformatics/getting-started/">https://gencore.bio.nyu.edu/bioinformatics/getting-started/</a></p>
<hr>
<p><b>Results for your recently completed sequencing run on flowcell {fcid} are available here:</b><br>
{delivery_dir}<br>
{run_dir_user_path}</p>

<p><b>All sequencing run and library statistics can be viewed in the interactive MultiQC report here:</b><br>
<a href="{mqc_report_url}">{mqc_report_url}</a><br>
{message}</p>

<p>Please let us know if you have any questions.</p>

<p>Best,<br>
GenCore Team</p>
    '''    
    return delivery_template

def deliver_raw_run_dir(run_dir_path, group):
    run_dir_name = os.path.basename(os.path.normpath(run_dir_path))
    raw_run_delivery_folder = raw_run_dir_delivery_root + group + "/" + run_dir_name
    run_dir_user_path = "<br>Raw Run Directory:<br>{}<br>".format(raw_run_delivery_folder)
    
    if os.path.exists(raw_run_delivery_folder):
        return run_dir_user_path
    
    #shutil.copytree(run_dir_path, raw_run_delivery_folder)
    #mode = 0o555 # This is the octal representation for r-xr-xr-x
    #change_permissions_recursive(raw_run_delivery_folder, mode)
    copy_command = subprocess.getoutput('mkdir -p {}; cp -rv {}/* {}/.'.format(raw_run_delivery_folder, run_dir_path, raw_run_delivery_folder))
    return run_dir_user_path


def deliver_data(fcid, path, lane_num, group, scheduled_date):
    delivery_dir = delivery_folder_root + "/" + group + "/" + scheduled_date + "_" + get_delivery_fcid(fcid) + "/" + lane_num
    #if os.path.exists(delivery_dir):
    #    shutil.rmtree(delivery_dir)
    #shutil.copytree(path, delivery_dir)
    #mode = 0o555 # This is the octal representation for r-xr-xr-x
    #change_permissions_recursive(delivery_dir, mode)
    copy_command = subprocess.getoutput('mkdir -p {}; rm -f {}/*; cp -v {}/* {}/.'.format(delivery_dir, delivery_dir, path, delivery_dir))
    return delivery_dir

def check_qc_and_deliver(path, summary_report_path):

    stats = parse_summary(summary_report_path)
    result = check_pool_errors(stats)
    success = result.get("success")
    error = result.get("message")

    summary_report_filename = os.path.basename(summary_report_path)
    fcid, lane_num = summary_report_filename.split("_")[:2]

    run = get_run_info(fcid)
    manufacturer = run['sequencer']['manufacturer']

    if success:
        print("qc_delivery.py: Passed QC")
        
        lanes = get_lanes(run["id"])["lanes"]
        
        # Get Lane and Pool Info
        i = lane_num
        if i == 'merged':
            i = 1
        
        lane = next(lane for lane in lanes if lane["lane_number"] == int(i))
        pool = get_pool(lane["id"])

        # Deliver Raw Run Directory if requested
        raw_run_dir_path = ''
        if run['deliver_run_dir']:
            raw_run_dir_path = deliver_raw_run_dir(get_run_dir(fcid)['run_dir'], pool['group']) 
            print("raw run dir delivered to: ", raw_run_dir_path)
            
        # Delivery Data
        scheduled_date = run['scheduled_date'].split("T")[0]
        delivery_dir = deliver_data(fcid, path, lane_num, pool['group'], scheduled_date)
        print("data delivered to: ", delivery_dir)
        
        # get delivery email
        message = ''
        if manufacturer == "Illumina":
            message = get_qc_messages(stats)
        mqc_report_url = "http://core-fastqc.bio.nyu.edu/{}/{}/multiqc_report.html".format(fcid, lane_num)
        delivery_email = get_delivery_email(fcid, delivery_dir, raw_run_dir_path, mqc_report_url, message, pool['no_demux'], pool['allowed_barcode_mismatch'])
        
        # update lane stats in tuboweb
        update_lane_stats(fcid, i, 'total_num_reads', stats["Total # of Single-End Reads"])
        update_lane_stats(fcid, i, 'total_num_pf_reads', stats["Total # PF Reads"])
        print("lane stats updated in tuboweb")
        
        # mark as devlivered in tuboweb (4 is to mark as delivered (is_ready_for_delivery=True))
        r = set_run_status(fcid, 4)
        print(r)
        print("run status updated in tuboweb")

        # send email
        pool_owner_email = f"{pool['created_by']}@nyu.edu"
        pi_email = f"{pool['pi_netid']}@nyu.edu" if pool['pi_netid'] else ''
        recipients = ['mk5636@nyu.edu', 'gencore-group@nyu.edu', pool_owner_email] + ([pi_email] if pi_email else [])
        #recipients = ['mk5636@nyu.edu']
        subject = "Data For " + fcid
        send_email(recipients, subject, delivery_email)
        print("email sent to: ", recipients)

    else:
        message = ''
        sequencer_name = run['sequencer']['name']
        first_demux_undetermined_pct = get_first_demux_undetermined_pct(fcid, 1)
        if (run['is_revcom_index2'] or first_demux_undetermined_pct is not None) \
            and ('NextSeq' in sequencer_name or 'NovaSeq' in sequencer_name):
            if first_demux_undetermined_pct is None:
                send_email(["mk5636@nyu.edu"], "ERROR For {}".format(fcid), 'first_demux_undetermined_pct is None') 
                print("first_demux_undetermined_pct is None")
                set_first_demux_undetermined_pct(fcid, 1, stats["% Undetermined"])
                flip_index2_revcom(fcid)
                run_redemux(fcid)
                message = "Resubmitted for demultiplexing. Undetermined after Demultiplex Attempt 1 = {}%".format(stats["% Undetermined"])
                send_email(["mk5636@nyu.edu"], "ERROR For {}".format(fcid), message)
            else:
                send_email(["mk5636@nyu.edu"], "ERROR For {}".format(fcid), "first_demux_undetermined_pct is not None it's: " + str(first_demux_undetermined_pct))
                print("first_demux_undetermined_pct is not None it's: " + str(first_demux_undetermined_pct))
                message = "Undetermined\nDemultiplex Attempt 1 = {}%\nDemultiplex Attempt 2 = {}%".format(first_demux_undetermined_pct, stats["% Undetermined"])
        else:
            message = error + "\nNo Auto Redemultiplexing\nhttp://core-fastqc.bio.nyu.edu/" + fcid

        send_email(["mk5636@nyu.edu"], "ERROR For {}".format(fcid), message)
    
def main():
    # qc_deliver.py <path_to_data_to_deliver> <path_to_sumary_report>
    # ex: qc_deliver.py /path/to/lane/FCID/1 /path/to/qc/FCID/1/sumary_report.txt
    # ex: qc_deliver.py /path/to/merged/FCID/merged path/to/qc/FCID/merged/sumary_report.txt
    path = sys.argv[1]
    summary_report_path = sys.argv[2]
    check_qc_and_deliver(path, summary_report_path)
        
if __name__ == "__main__":
    main()
