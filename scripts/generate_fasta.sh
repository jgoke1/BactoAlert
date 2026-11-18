if [[ $# -lt 2 ]]; then usage >&2; exit 1; fi

input_dir=$1
output_csv=$2
shift 2
strip_suffix=''
force=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --strip-suffix)
            [[ $# -ge 2 && -n $2 ]] || die '--strip-suffix needs a value.'
            strip_suffix=$2
            shift 2
            ;;
        --force) force=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done

[[ -d "$input_dir" ]] || die "Input directory not found: $input_dir"
[[ -r "$input_dir" && -x "$input_dir" ]] || die "Cannot read directory: $input_dir"
input_dir=$(cd -- "$input_dir" && pwd -P)
[[ "$output_csv" == *.csv ]] || die 'The output filename must end in .csv.'
[[ ! -d "$output_csv" ]] || die "Output is a directory: $output_csv"
if [[ -e "$output_csv" || -L "$output_csv" ]]; then
    [[ $force -eq 1 ]] || die "Output exists: $output_csv. Use --force to replace it."
fi

output_dir=$(dirname -- "$output_csv")
mkdir -p -- "$output_dir"
temporary_csv=$(mktemp "$output_dir/.fasta_samplesheet.XXXXXX")
trap 'rm -f -- "$temporary_csv"' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
printf 'sample_id,fasta\n' > "$temporary_csv"

count=0
seen_ids='|'
shopt -s nullglob

for fasta in "${input_dir%/}"/*; do
    [[ ! -d "$fasta" ]] || continue
    case "$fasta" in
        *.fa|*.fna|*.fasta|*.fa.gz|*.fna.gz|*.fasta.gz) ;;
        *) continue ;;
    esac
    [[ -f "$fasta" && -r "$fasta" && -s "$fasta" ]] || die "FASTA is missing, unreadable or empty: $fasta"
    [[ ! "$fasta" =~ [[:space:]] ]] || die "Your pipeline requires FASTA paths without whitespace: $fasta"

    sample_id=${fasta##*/}
    sample_id=${sample_id%.gz}
    sample_id=${sample_id%.*}
    if [[ -n "$strip_suffix" ]]; then
        sample_id=${sample_id%"$strip_suffix"}
    fi
    [[ "$sample_id" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || die "Invalid sample ID '$sample_id' from: $fasta"
    [[ "$seen_ids" != *"|${sample_id}|"* ]] || die "Duplicate sample ID '$sample_id'; keep one FASTA per sample or use distinct filenames."
    seen_ids="${seen_ids}${sample_id}|"

    # Escape embedded quotes so paths containing commas/quotes remain valid CSV.
    csv_path=${fasta//\"/\"\"}
    printf '%s,"%s"\n' "$sample_id" "$csv_path" >> "$temporary_csv"
    count=$((count + 1))
done

[[ $count -gt 0 ]] || die "No supported FASTA files found directly inside: $input_dir"
mv -f -- "$temporary_csv" "$output_csv"
trap - EXIT INT TERM
printf 'Wrote %s samples to %s\n' "$count" "$output_csv"
