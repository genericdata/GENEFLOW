#!/usr/bin/env python3

import os
import re
import sys
import json
import glob
from collections import defaultdict
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary
from slime import check_do_merge, check_demux, get_run_dir, get_num_lanes
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
    return metrics


def generate_reports(fcid):
    run_dir = get_run_dir(fcid)['run_dir']
    num_lanes = get_num_lanes(fcid)
    sav = parse_sav(run_dir)
    print("sav: {}".format(sav))
    phix_aligned_per_lane = sav["phix_aligned_dict"]
    total_num_reads = sav["total_num_reads"]
    total_pf_reads = sav["total_pf_reads"]

    lanes_dict = {}
    for lane in range(1, num_lanes + 1):
        lanes_dict[lane] = {"lane_read_count_total_unfiltered": total_num_reads[lane]}
        no_demux = check_demux(fcid, lane)
        if not no_demux:
            pheniqs_out = get_pheniqs_output(fcid, lane)
            lanes_dict[lane].update(parse_pheniqs_output(pheniqs_out))

    do_merge = check_do_merge(fcid)
    prepare_demux_reports(lanes_dict, do_merge, fcid)
    prepare_summary_reports(lanes_dict, do_merge, fcid, phix_aligned_per_lane, total_pf_reads)


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


def prepare_summary_reports(lanes_dict, do_merge, fcid, phix_aligned_per_lane,total_pf_reads):
    if do_merge:
        path = "merged/{}_merged_summary_mqc.txt".format(fcid)
        fc_read_count = get_flowcell_read_count(lanes_dict)
        phix_aligned = sum(phix_aligned_per_lane.values()) / len(phix_aligned_per_lane)
        no_demux = check_demux(fcid, 1)
        report = create_summary_report(
            "Number of Lanes",
            len(lanes_dict),
            fc_read_count["total_unfiltered"],
            sum(total_pf_reads.values()) if no_demux else fc_read_count["total_filtered"],
            "" if no_demux else fc_read_count["undetermined"],
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
                lanes_dict[l]["lane_read_count_total_unfiltered"],
                total_pf_reads[l] if no_demux else lanes_dict[l]["lane_read_count_total_filtered"],
                "" if no_demux else lanes_dict[l]["undetermined"],
                phix_aligned_per_lane[l],
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


def get_flowcell_read_count(lanes):
    fc_read_count_total_unfiltered = 0
    fc_read_count_total_filtered = 0
    fc_undetermined = 0
    for l in sorted(lanes):
        fc_read_count_total_unfiltered += lanes[l]["lane_read_count_total_unfiltered"]
        fc_read_count_total_filtered += lanes[l]["lane_read_count_total_filtered"]
        fc_undetermined += lanes[l]["undetermined"]
    fc_undetermined = round((fc_undetermined / len(lanes)), 2)

    return {
        "total_unfiltered": fc_read_count_total_unfiltered,
        "total_filtered": fc_read_count_total_filtered,
        "undetermined": fc_undetermined,
    }


def parse_pheniqs_output(report):
    with open(report, "r") as f:
        data = f.read()
    json_data = json.loads(data)

    lane_dict = {"libs": {}}
    read_count_total = 0
    portion_total = 0

    for r in json_data["demultiplex output report"]["read group quality reports"]:
        lib_name = r["LB"]
        lane_dict["libs"][lib_name] = {
            "yield": r["pf count"],
            "portion": round(r["pf pooled fraction"] * 100, 3),
        }
        read_count_total += r["pf count"]
        portion_total += r["pf pooled fraction"] * 100

    lane_dict["lane_read_count_total_filtered"] = json_data[
        "demultiplex output report"
    ]["pf count"]
    lane_dict["undetermined"] = (
        round((1 - json_data["demultiplex output report"]["multiplex fraction"])*100,3)
    )
    lane_dict["portion_total"] = round(portion_total, 0)

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
