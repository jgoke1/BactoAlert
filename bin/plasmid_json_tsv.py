import json
import sys
import csv

def plasmidfinder_json_to_tsv(json_file, tsv_file):
    with open(json_file) as f:
        data = json.load(f)

    seq_regions = data.get("seq_regions", {})

    headers = [
        "plasmid",
        "identity",
        "coverage",
        "alignment_length",
        "ref_gene_length",
        "depth",
        "ref_id",
        "ref_acc",
        "database"
    ]

    with open(tsv_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(headers)

        # If no plasmids found, write a dummy row
        if not seq_regions:
            writer.writerow(["NA"] * (len(headers) - 1) + ["PlasmidFinder"])
            return

        for _, entry in seq_regions.items():
            writer.writerow([
                entry.get("name", "NA"),
                entry.get("identity", "NA"),
                entry.get("coverage", "NA"),
                entry.get("alignment_length", "NA"),
                entry.get("ref_gene_lenght", "NA"),  # typo is real in PlasmidFinder
                entry.get("depth", "NA"),
                entry.get("ref_id", "NA"),
                entry.get("ref_acc", "NA"),
                "PlasmidFinder"
            ])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python plasmidfinder_json_to_tsv.py <input.json> <output.tsv>")

    plasmidfinder_json_to_tsv(sys.argv[1], sys.argv[2])

