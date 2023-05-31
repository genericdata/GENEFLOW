import sys
import gzip
from collections import defaultdict


def get_barcodes(fastq_path):
    barcodes = defaultdict(int)
    with gzip.open(fastq_path, "rt") as fastq:
        for line in fastq:
            if line.startswith("@"):
                barcode = line.rstrip().split(":")[-1]
                barcodes[barcode] += 1
    return barcodes


def print_barcodes(barcodes):
    total = sum(barcodes.values())
    for barcode, count in sorted(
        barcodes.items(), key=lambda item: item[1], reverse=True
    ):
        print(f"{barcode}\t{count}\t{round(count/total*100, 2)}")


def main():
    barcodes = get_barcodes(sys.argv[1])
    print_barcodes(barcodes)


if __name__ == "__main__":
    main()
