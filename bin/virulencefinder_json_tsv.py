import json
import sys
import csv

def virulencefinder_json_to_tsv(json_file, tsv_file):
    with open(json_file) as f:
        data = json.load(f)

    seq_regions = data.get("seq_regions", {})
    phenotypes = data.get("phenotypes", {})

    headers = [
        "gene",
        "function",
        "identity",
        "coverage",
        "alignment_length",
        "ref_seq_length",
        "depth",
        "grade",
        "ref_id",
        "ref_acc",
        "database"
    ]

    with open(tsv_file, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(headers)

        # Handle case where no virulence genes are found
        if not seq_regions:
            writer.writerow(
                ["NA"] * (len(headers) - 1) + ["VirulenceFinder"]
            )
            return

        for _, entry in seq_regions.items():
            gene = entry.get("name", "NA")

            # Look up gene function from phenotypes block
            function = phenotypes.get(gene, {}).get("function", "NA")

            writer.writerow([
                gene,
                function,
                entry.get("identity", "NA"),
                entry.get("coverage", "NA"),
                entry.get("alignment_length", "NA"),
                entry.get("ref_seq_length", "NA"),
                entry.get("depth", "NA"),
                entry.get("grade", "NA"),
                entry.get("ref_id", "NA"),
                entry.get("ref_acc", "NA"),
                "VirulenceFinder"
            ])


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python virulencefinder_json_to_tsv.py <input.json> <output.tsv>"
        )

    virulencefinder_json_to_tsv(sys.argv[1], sys.argv[2])

