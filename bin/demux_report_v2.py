#!/usr/bin/env python3

import os
import re
import sys
import json
import glob
from collections import defaultdict
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary
from slime import check_do_merge, check_demux, get_run_dir, get_num_lanes, get_run_info
from config import alpha


def compute_lane_metrics_from_sav(summary):
    phix_aligned_per_lane = {}
    total_reads_per_lane = {}
    total_pf_reads_per_lane = {}

    lane_count = summary.at(0).lane_count()

    for lane_number in range(lane_count):
        lane_index = lane_number + 1

        total_reads_per_lane[lane_index] = summary.at(0).at(lane_number).reads()
        total_pf_reads_per_lane[lane_index] = summary.at(0).at(lane_number).reads_pf()

        num_reads, aligned_total = calculate_phix_alignment_for_lane_from_sav(
            summary, lane_number
        )
        phix_aligned_per_lane[lane_index] = (
            aligned_total / num_reads if num_reads != 0 else 0
        )

    return {
        "phix_aligned_dict": phix_aligned_per_lane,
        "total_num_reads": total_reads_per_lane,
        "total_pf_reads": total_pf_reads_per_lane,
    }


def calculate_phix_alignment_for_lane_from_sav(summary, lane_number):
    num_reads = 0
    aligned_total = 0

    for read_index in range(summary.size()):
        read_summary = summary.at(read_index)

        # Skip index reads
        if read_summary.read().is_index():
            continue

        num_reads += 1
        aligned_total += read_summary.at(lane_number).percent_aligned().mean()

    return num_reads, aligned_total


def parse_sav(run_folder):
    ## Magic Code
    run_metrics = py_interop_run_metrics.run_metrics()
    valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
    py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
    run_folder = run_metrics.read(run_folder, valid_to_load)
    summary = py_interop_summary.run_summary()
    py_interop_summary.summarize_run_metrics(run_metrics, summary)
    metrics = compute_lane_metrics_from_sav(summary)
    '''
    metrics = {}
    metrics['phix_aligned_dict']={1: 0, 2:0, 3:0, 4:0}
    metrics['total_num_reads']={1: 0, 2:0, 3:0, 4:0}
    metrics['total_pf_reads']={1: 0, 2:0, 3:0, 4:0}
    '''
    return metrics


def parse_aviti_run_stats(run_dir):
    """
    Read run_dir/AvitiRunStats.json and return a metrics dict keyed by lane.
    For each lane it extracts:
      - phix_aligned_dict: average PhiXAlignmentRate across all reads (R1, R2, etc.)
      - total_num_reads:   PolonyCount
      - total_pf_reads:    PFCount
    Supports any number of lanes.
    """
    stats_file = os.path.join(run_dir, "AvitiRunStats.json")
    with open(stats_file, "r") as f:
        data = json.load(f)

    lane_stats = data.get("LaneStats", [])
    metrics = {
        "phix_aligned_dict": {},
        "total_num_reads":   {},
        "total_pf_reads":    {},
    }

    for i, rec in enumerate(lane_stats):
        # JSON uses zero‐based lane indices; convert to 1-based
        lane = i + 1

        # per‐lane counts
        polony = rec.get("PolonyCount", 0)
        pf     = rec.get("PFCount",    0)

        # average PhiX alignment over all read entries
        reads = rec.get("Reads", [])
        rates = [
            r.get("PhiXAlignmentRate", 0)
            for r in reads
            if "PhiXAlignmentRate" in r
        ]
        phi_avg = (sum(rates) / len(rates)) if rates else 0.0

        # populate dicts
        metrics["total_num_reads"][lane]   = polony
        metrics["total_pf_reads"][lane]    = pf
        metrics["phix_aligned_dict"][lane] = phi_avg

    return metrics

def get_metrics(fcid):
    run = get_run_info(fcid)
    run_dir = get_run_dir(fcid)['run_dir']
    manufacturer = run['sequencer']['manufacturer']
    if manufacturer == "Element":
        metrics = parse_aviti_run_stats(run_dir)
    else:
        metrics = parse_sav(run_dir)
    
    return metrics


def generate_reports(fcid):
    num_lanes = get_num_lanes(fcid) 
    metrics = get_metrics(fcid)
    print(f"metrics: {metrics}")

    lanes_dict = {}
    for lane in range(1, num_lanes + 1):
        lanes_dict[lane] = {}
        no_demux = check_demux(fcid, lane)
        if not no_demux:
            pheniqs_out = get_pheniqs_output(fcid, lane)
            lanes_dict[lane].update(parse_pheniqs_output(pheniqs_out))

    do_merge = check_do_merge(fcid)
    prepare_demux_reports(lanes_dict, do_merge, fcid)
    # prepare_summary_reports() needs lanes_dict because it has the undetermined reads per lane
    prepare_summary_reports(lanes_dict, do_merge, fcid, metrics)


def create_demux_report(libs, num_lanes):
    report = "\n".join(
        [
            "# id: 'demux_report'",
            "# plot_type: 'table'",
            "# section_name: 'Demultiplexing Report'",
            '# description: "<br><strong>Total Read Count:</strong> Total number of PF (Passing Filter) reads in this library.<br><strong>Portion:</strong> The proportion of reads that represent the individual library in the entire Library Pool."',
            "Library\tTotal Read Count\tPortion (%)",
            *(
                f"{lib}\t{libs[lib]['yield']:,d}\t{libs[lib]['portion']/num_lanes:.1f}"
                for lib in libs
            ),
        ]
    )
    return report


def create_summary_report(
    lane_header,
    lane_col,
    read_count_total_unfiltered,
    read_count_total_filtered,
    undetermined,
    phix_aligned,
    no_demux
):
    undetermined_header = "" if no_demux else "\t% Undetermined"
    undetermined_value = "" if no_demux else "\t{:.2f}".format(undetermined)
    
    summary = "\n".join(
        [
            "# id: 'run_stats'",
            "# plot_type: 'table'",
            "# section_name: 'Run Statistics'",
            f"{lane_header}\tTotal # of Single-End Reads\tTotal # PF Reads{undetermined_header}\t% PhiX Aligned",
            f"{lane_col}\t{int(read_count_total_unfiltered):,d}\t{int(read_count_total_filtered):,d}{undetermined_value}\t{phix_aligned:.2f}",
        ]
    )
    return summary


def prepare_summary_reports(lanes_dict, do_merge, fcid, metrics):
    if do_merge:
        path = "merged/{}_merged_summary_mqc.txt".format(fcid)
        no_demux = check_demux(fcid, 1)        
        phix_aligned = sum(metrics["phix_aligned_dict"].values()) / len(metrics["phix_aligned_dict"])
        undetermined = (
            "" if no_demux
            else round(sum(l["undetermined"] for l in lanes_dict.values()) / len(lanes_dict), 2)
        )

        report = create_summary_report(
            "Number of Lanes",
            len(lanes_dict),
            sum(metrics["total_num_reads"].values()),
            sum(metrics["total_pf_reads"].values()),
            undetermined,
            phix_aligned,
            no_demux
        )
        write_report(report, path)

    else:
        for l in sorted(lanes_dict):
            path = "{}/{}_{}_summary_mqc.txt".format(l, fcid, l)
            no_demux = check_demux(fcid, l)
            report = create_summary_report(
                "Lane",
                l,
                metrics["total_num_reads"][l],
                metrics["total_pf_reads"][l],
                "" if no_demux else lanes_dict[l]["undetermined"],
                metrics["phix_aligned_dict"][l],
                no_demux
            )
            write_report(report, path)


def prepare_demux_reports(lanes_dict, do_merge, fcid):
    if do_merge:
        no_demux = check_demux(fcid, 1)
        if no_demux:
            return
        path = "merged/{}_merged_demux_report_mqc.txt".format(fcid)
        combined = defaultdict(dict)
        for l in sorted(lanes_dict):
            for lib, lib_data in lanes_dict[l]["libs"].items():
                combined[lib]["yield"] = (
                    combined[lib].get("yield", 0) + lib_data["yield"]
                )
                combined[lib]["portion"] = (
                    combined[lib].get("portion", 0) + lib_data["portion"]
                )
        report = create_demux_report(combined, len(lanes_dict))
        write_report(report, path)
    else:
        for l in sorted(lanes_dict):
            no_demux = check_demux(fcid, l)
            if no_demux:
                continue
            path = "{}/{}_{}_demux_report_mqc.txt".format(l, fcid, l)
            report = create_demux_report(lanes_dict[l]['libs'], 1)
            write_report(report, path)


def write_report(report, path):
    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)

    out = open(path, "w")
    out.write(report)
    out.close()

    print("Wrote report: {}".format(path))


def parse_pheniqs_output(report):
    with open(report, "r") as f:
        data = f.read()
    json_data = json.loads(data)
    lane_dict = {"libs": {}}

    for r in json_data["sample"]["classified"]:
        lib_name = str(r["LB"])
        lane_dict["libs"][lib_name] = {
            "yield": r["count"],
            "portion": round(r["pooled fraction"] * 100, 3),
        }

    lane_dict['libs']['undetermined'] = {
        "yield": json_data["sample"]["unclassified"]['count'],
        "portion": round(json_data["sample"]["unclassified"]["pooled fraction"] * 100, 3)
    }

    lane_dict["lane_read_count_total_filtered"] = json_data["sample"]["count"]
    lane_dict["undetermined"] = (
        round(json_data["sample"]["unclassified"]["pooled fraction"] * 100,3)
    )

    return lane_dict


def get_pheniqs_output(fcid, lane):
    file_path = os.path.join(alpha, "pheniqs_out", fcid, str(lane), 'demux.out')

    if os.path.isfile(file_path):
        return file_path
    else:
        print(f'demux_report.py: get_pheniqs_output: {file_path} does not exist.')


def main():    
    fcid = sys.argv[1]
    generate_reports(fcid)


if __name__ == "__main__":
    main()
