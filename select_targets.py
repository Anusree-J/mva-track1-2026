"""Select carried CDS or canonical two-base splice overlaps without ranking."""

import argparse
import csv
import os
from pathlib import Path


def main():
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.input.open() as source, args.output.open("x", newline="") as dest:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(dest, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()
        count = 0
        for row in reader:
            if row["allele_carriage"] == "carried" and set(row["location_flags"].split(";")) & {
                "CDS",
                "splice_acceptor_2bp",
                "splice_donor_2bp",
            }:
                writer.writerow(row)
                count += 1
    print(f"PASS: selected {count} carried coding/splice rows")


if __name__ == "__main__":
    main()
