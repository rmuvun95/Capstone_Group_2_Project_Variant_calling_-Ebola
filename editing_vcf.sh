#!/bin/bash

# Process only original VCF files (not *_filtered.vcf.gz)
for file in results/vcf/*.vcf.gz; do
    # Skip files that are already filtered
    if [[ "$file" == *_filtered.vcf.gz ]]; then
        continue
    fi

    # Get sample ID from original filename
    sample_id=$(basename "$file" .vcf.gz)

    # Define filtered filename
    filtered_file="results/vcf/${sample_id}_filtered.vcf.gz"

    # Apply filters
    bcftools filter -e 'QUAL<30 || DP<10 || AC<0.05' "$file" -Oz -o "$filtered_file"

    # Index the filtered file
    tabix -p vcf "$filtered_file"
done

# Extract info from filtered files and write per-sample TSVs
for filtered_file in results/vcf/*_filtered.vcf.gz; do
    # Get sample ID by stripping _filtered.vcf.gz
    sample_id=$(basename "$filtered_file" _filtered.vcf.gz)

    output_tsv="results/vcf/${sample_id}_hq.tsv"

    (
        echo -e "CHROM\tPOS\tREF\tALT\tQUAL\tDP\tAF\tSAMPLE_ID"
        zgrep -v "^##" "$filtered_file" | awk -v sid="$sample_id" -F'\t' 'BEGIN{OFS="\t"}
            NR==1 { next }
            {
                split($8, info, ";");
                dp=ac=an="";
                for (i in info) {
                    if (info[i] ~ /^DP=/) dp=substr(info[i], 4);
                    if (info[i] ~ /^AC=/) ac=substr(info[i], 4);
                    if (info[i] ~ /^AN=/) an=substr(info[i], 4);
                }
                if (dp && ac && an && an != 0) {
                    af = ac / an;
                    print $1, $2, $4, $5, $6, dp, af, sid;
                }
            }
        '
    ) > "$output_tsv"
done


# Concatenate all tsv without headers

 awk 'FNR==1 && NR!=1 {next} {print}' results/vcf/*_hq.tsv > results/vcf/hq_allvariants.tsv
