#!/usr/bin/env python3

import os
import sys
import requests
import json
import glob
import smtplib
from xml.dom import minidom
from shutil import copyfile, copytree
from config import tw_api_root, tw_user, tw_api_key, gmail_user, gmail_pwd, raw_run_root


def send_email(recipient, subject, body):
    FROM = gmail_user
    TO = recipient if isinstance(recipient, list) else [recipient]
    SUBJECT = "GENEFLOW: " + subject
    TEXT = body

    # Prepare actual message
    message = f"""From: {FROM}\nTo: {", ".join(TO)}\nSubject: {SUBJECT}\n\n{TEXT}"""

    try:
        with smtplib.SMTP("smtp.gmail.com", 587) as server:
            server.ehlo()
            server.starttls()
            server.login(gmail_user, gmail_pwd)
            server.sendmail(FROM, TO, message)

    except Exception as e:
        print("Failed to send email:", str(e))


def get_run_info(fcid):
    run_url = f"{tw_api_root}flowcell/illumina/{fcid}/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(run_url, params=params)
    run_data = response.json()
    return run_data


def get_num_lanes(fcid):
    run_url = f"{tw_api_root}flowcell/illumina/{fcid}/num_lanes"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(run_url, params=params)
    run_data = response.json()
    return run_data


def get_lanes(run_id):
    lanes_url = f"{tw_api_root}flowcell/{run_id}/lanes/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(lanes_url, params=params)
    lanes_data = response.json()
    return lanes_data


def get_pool(lane_id):
    pools_url = f"{tw_api_root}lane/{lane_id}/librarypools/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(pools_url, params=params)
    pools_data = response.json()
    
    if len(pools_data["librarypools"]) == 0:
        # todo: check if it's PHiX
        print(f"There is no Library Pool for lane: {lane_id}")
        return None

    pool_id = pools_data["librarypools"][0]["id"]
    pool_url = f"{tw_api_root}librarypool/{pool_id}/"
    response = requests.get(pool_url, params=params)
    pool_data = response.json()

    return pool_data


def get_libraries(pool_id):
    lpc_url = f"{tw_api_root}librarypool/{pool_id}/librarypoolconcentrations/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(lpc_url, params=params)
    lpc_data = response.json()
    return lpc_data["librarypoolconcentrations"]


def get_library(lib_id):
    url = f"{tw_api_root}library/{lib_id}/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(url, params=params)
    data = response.json()
    return data


def get_lib_barcodes(lib_id):
    lbc_url = f"{tw_api_root}library/{lib_id}/barcodes/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(lbc_url, params=params)
    lbc_data = response.json()

    barcodes = []
    for bc in lbc_data["barcodes"]:
        url = f"{tw_api_root}barcode/{bc['id']}/"
        response = requests.get(url, params=params)
        data = response.json()
        barcodes.append(data)

    return barcodes


def get_run_dir(fcid):
    query = f"{raw_run_root}/*/*{fcid}"
    print("Locating run dir: ", query)
    run_dir = glob.glob(os.path.join(query))
    
    if len(run_dir) != 1:
        print({'status': 'error', 'msg': f'len(run_dir) != 1, found {len(run_dir)}; fcid = {fcid}'})
        sys.exit(1)
    
    print("Found run dir: ", run_dir[0])
    return {'status': 'success', 'run_dir': run_dir[0]}


def check_do_merge(fcid):
    run = get_run_info(fcid)
    run_type = run["run_type_name"]
    return (
        "NextSeq" in run["sequencer"]["name"] or "NovaSeq" in run["sequencer"]["name"]
    ) and not run_type.startswith("XP")


# Check if the pool in this lane needs to be demultiplexed or not
def check_demux(fcid, lane_num):
    try:
        run_id = get_run_info(fcid)["id"]
        lanes = get_lanes(run_id)['lanes']
        lane = next(lane for lane in lanes if lane["lane_number"] == int(lane_num))
        pool = get_pool(lane["id"])

        # "pool.get("no_demux", True) is not False" checks if pool.no_demux is True, None, 
        # or if it does not exist in the pool. It defaults to True.
        # Using "is not False" ensures that we cover the cases where pool.no_demux is True or None.
        # If you simply use "is True", it won't cover the scenario when pool.no_demux is None or doesn't exist in pool.
        # Thus, return True if pool.no_demux is True, None, if there's an error, or if pool.no_demux doesn't exist.
        return pool is None or pool.get("no_demux", True) is not False
    except Exception as e:
        print(f"An error occurred: {e}")
        # In case of error, return True
        #return True  
        sys.exit(1)


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement = ''.join(complement[base] for base in reversed(seq))
    return reverse_complement


def set_first_demux_undetermined_pct(fcid, num, val):
    url = f"{tw_api_root}flowcell/illumina/{fcid}/set_first_demux_undetermined_pct/{num}"
    params = {"username": tw_user, "api_key": tw_api_key, "value": val}
    response = requests.get(url, params=params)
    data = response.json()
    return data
    
    
def get_first_demux_undetermined_pct(fcid, num):
    url = f"{tw_api_root}flowcell/illumina/{fcid}/get_first_demux_undetermined_pct/{num}"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(url, params=params)
    data = response.json()
    return data

def run_redemux(fcid):
    # Define the base command and parameters
    root_log_dir = alpha + "/logs/"
    
    base_cmd = 'sbatch'
    output_cmd = f"--output={root_log_dir}/{fcid}/pipeline/slurm-%j-redemux.out"
    error_cmd = f"--error={root_log_dir}/{fcid}/pipeline/slurm-%j-redemux.err"
    job_name_cmd = f"--job-name=GENEFLOW_MANAGER_({fcid})_redemux"
    
    # Find the run directory
    find_cmd = f"find -maxdepth 2 -type d -name '*{fcid}'"
    rundir_cmd = f"rundir=$({find_cmd}); echo rundir = $rundir;"
    
    # Construct the launch command with the found run directory
    # TODO: make path to launch.sh a param in the config?
    launch_cmd = f"/home/gencore/SCRIPTS/GENEFLOW/launch.sh /scratch/gencore/sequencers/$rundir {fcid} demux"
    
    # Combine the commands
    full_cmd = f"{rundir_cmd} {base_cmd} {output_cmd} {error_cmd} {job_name_cmd} {launch_cmd}"
    
    # Execute the command and get the output
    demux_output = subprocess.getoutput(full_cmd)
    
    # Format and print the output
    formatted_output = f"\n{full_cmd}\n{demux_output}\n"
    print(f"Auto-Redemuxing\n{formatted_output}")
    
    
def set_pool_index_revcom(pool_id, i):
    url = f"{tw_api_root}librarypool/{pool_id}/setindexisrevcom/{i}/"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(url, params=params)
    data = response.json()
    return data


def set_run_index_revcom(fcid, status):
    url = f"{tw_api_root}flowcell/illumina/{fcid}/setindexisrevcom/{status}"
    params = {"username": tw_user, "api_key": tw_api_key}
    response = requests.get(url, params=params)
    data = response.json()
    return data


def flip_index2_revcom(fcid):
    """
    Turn on/off revcom for index 2 for the run and all pools.
    Right now we just turn off, but we can set this to be a parameter
    or to toggle it based on the current state.
    """
    run = get_run_info(fcid)
    lanes = get_lanes(run['id'])
    for lane in lanes:
        pool = get_pool(lane['pool']['id'])
        # turn off revcom for index 2 for all pools
        set_pool_index_revcom(pool, 0)
    # turn off revcom for index 2 for the run
    set_run_index_revcom(fcid, 0)
