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
        message += "\nUnfortunately, your lane(s) received less reads than expected. This is likely related to an issue with metadata, quantification and/or pooling. We would be happy to discuss this further if desired.\n"

    if (
        "% Undetermined" in stats
        and stats["% Undetermined"] - stats["% PhiX Aligned"] > 15
    ):
        message += "\nThere is a large number of 'undetermined reads' in your lane(s), those whose index cannot be appropriately identified as any of the ones expected for the libraries. This is likely related to an issue with metadata, pooling, and/or library prep. We would be happy to discuss this further if desired. \n"

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


def get_delivery_email(fcid, delivery_dir, run_dir_user_path, mqc_report_url, message):
    delivery_template = f'''
Dear GenCore Users,

Results for your recently completed sequencing run on flowcell {fcid} are available here:
{delivery_dir}
{run_dir_user_path}
All sequencing run and library statistics can be viewed in the interactive MultiQC report here:
{mqc_report_url}
{message}
Please let us know if you have any questions.

Best,
GenCore Team

---
Note: You must have the required permissions to access data on the HPC. 
If this is your first time sequencing, please visit: https://gencore.bio.nyu.edu/bioinformatics/getting-started/
    '''
    
    return delivery_template


def deliver_raw_run_dir(run_dir_path, group):
    run_dir_name = os.path.basename(os.path.normpath(run_dir_path))
    raw_run_delivery_folder = raw_run_dir_delivery_root + group + "/" + run_dir_name
    run_dir_user_path = "\nRaw Run Directory:\n{}\n".format(raw_run_delivery_folder)
    
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
        message = get_qc_messages(stats)
        mqc_report_url = "http://core-fastqc.bio.nyu.edu/{}/{}/multiqc_report.html".format(fcid, lane_num)
        delivery_email = get_delivery_email(fcid, delivery_dir, raw_run_dir_path, mqc_report_url, message)
        
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
        if run['is_revcom_index2'] and ('NextSeq' in sequencer_name or 'NovaSeq' in sequencer_name):
            first_demux_undetermined_pct = get_first_demux_undetermined_pct(fcid, 2)
            if first_demux_undetermined_pct is None:
                set_first_demux_undetermined_pct(fcid, 2, stats["% Undetermined"])
                flip_index2_revcom(fcid)
                run_redemux(fcid)
                message = "Resubmitted for demultiplexing. Undetermined after Demultiplex Attempt 1 = {}%".format(stats["% Undetermined"])
            else:
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
