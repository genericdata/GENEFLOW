import requests
import json
from xml.dom import minidom
from config import tw_api_root, tw_user, tw_api_key


def get_run_info(fcid):
    run_url = f"{tw_api_root}flowcell/illumina/{fcid}/"
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


def check_do_merge(fcid):
    run = get_run_info(fcid)
    run_type = run["run_type_name"]
    return (
        "NextSeq" in run["sequencer"]["name"] or "NovaSeq" in run["sequencer"]["name"]
    ) and not run_type.startswith("XP")


# Check all lanes for the value of demux for the pool in that lane
# either they will all have the same value (miseq, nextseq, some hiseqs)
# or occasionally a hiseq will have lanes with different specs for demux
# report (return) whether lanes are demux, no demux, or inconsistent
def check_demux(run_id):
    lanes = get_lanes(run_id)
    demux_values = set()

    for lane in lanes["lanes"]:
        pool = get_pool(lane["id"])
        if pool is None:
            continue

        demux_value = pool.get("no_demux")
        if demux_value is None:
            continue

        demux_values.add(str(demux_value))

        if len(demux_values) > 1:
            return "error: inconsistent"

    return demux_values.pop() if len(demux_values) == 1 else "error: inconsistent"
