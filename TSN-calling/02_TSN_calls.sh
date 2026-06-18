# Exact TSN as based on PRO-seq 5 pnt coverage

# 3 steps
  # 1: defining TSS search regions (R script)
    # 01_TSS_regions.R
  # 2: calling 5pnt (bash)
    # 02_TSN_calls.sh
  # 3: filter and update TSS/TSN annotation (R script)
    # 03_TSN_filter_update.R


####################################################################################################
# 02_TSN_calls.sh
#
# Language: Bash/Unix
#
# Aim:
#   Call the position with maximum 5' end coverage within each transcription start region (see 01_TSS_regions.R)
#
# Input:
#   tss_region_plus.bed
#   tss_region_minus.bed
#   merged_rep1-all_5pnt_nonsorted.bed
#   merged_rep2-all_5pnt_nonsorted.bed
#
# Output:
#   TSS_max_cov_and_pos_perregion_merged_rep1-all.txt
#   TSS_max_cov_and_pos_perregion_merged_rep2-all.txt
####################################################################################################

#!/bin/bash -l
#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -N call_5pnt_TSS
#$ -e ./call_5pnt_TSS.err.txt
#$ -pe smp 8
#$ -q all.q@narrativum.umcn.nl


# Directories

BASE_DIR="/home/femkevh/PROseq/Aag2"
TSS_DIR="${BASE_DIR}/tss"
PARALLEL_DIR="${TSS_DIR}/parallel"
BED_DIR="${BASE_DIR}/bed/mergedreplicates"

cd "${PARALLEL_DIR}"


# Input files

TSS_PLUS="${TSS_DIR}/tss_region_plus.bed"
TSS_MINUS="${TSS_DIR}/tss_region_minus.bed"

REPLICATES="
merged_rep1-all
merged_rep2-all
"


# Sort the TSS region files 
sort -k1,1 -k2,2n "${TSS_PLUS}" > sorted_tss_regions_plus.bed
sort -k1,1 -k2,2n "${TSS_MINUS}" > sorted_tss_regions_minus.bed


# Sort the PRO-seq 5pnt coverage files per replicate
for REP in ${REPLICATES}
do
  echo "Preparing ${REP}"

  sort -k1,1 -k2,2n "${BED_DIR}/${REP}_5pnt_nonsorted.bed" > "${REP}_5pnt.bed"

  awk '$6 == "+"' "${REP}_5pnt.bed" > "${REP}_5pnt_plus.bed"
  awk '$6 == "-"' "${REP}_5pnt.bed" > "${REP}_5pnt_minus.bed"
done




###################################################
# Call TSN for plus-strand genes
###################################################

for REP in ${REPLICATES}
do

  OUTPUT_FILE="plus_TSS_max_cov_and_pos_perregion_${REP}.txt"

  echo -e "chr\tstartsite-TSS\tendsite-TSS\tmaxcoverage_position\tmaxcoverage_coverage\tgenename-isoform" > "${OUTPUT_FILE}"

  FIVEPNT_FILE="${REP}_5pnt_plus.bed"

  while read chrom start end name
  do
    tss_coverage=$(awk -v c="${chrom}" -v s="${start}" -v e="${end}" '$1==c && $2>=s && $3<=e' "${FIVEPNT_FILE}")

    maxpos=$(echo "${tss_coverage}" | cut -f2 | uniq -c | sort -rn | head -n1 | awk '{print $2}')
    maxcov=$(echo "${tss_coverage}" | cut -f2 | uniq -c | sort -rn | head -n1 | awk '{print $1}')

    echo -e "${chrom}\t${start}\t${end}\t${maxpos}\t${maxcov}\t${name}" >> "${OUTPUT_FILE}"

  done < sorted_tss_regions_plus.bed
done



###################################################
# Call TSN for minus-strand genes
###################################################
for REP in ${REPLICATES}
do

  OUTPUT_FILE="minus_TSS_max_cov_and_pos_perregion_${REP}.txt"

  echo -e "chr\tstartsite-TSS\tendsite-TSS\tmaxcoverage_position\tmaxcoverage_coverage\tgenename-isoform" > "${OUTPUT_FILE}"

  FIVEPNT_FILE="${REP}_5pnt_minus.bed"

  while read chrom start end name
  do
    tss_coverage=$(awk -v c="${chrom}" -v s="${start}" -v e="${end}" '$1==c && $2>=s && $3<=e' "${FIVEPNT_FILE}")

    maxpos=$(echo "${tss_coverage}" | cut -f2 | uniq -c | sort -rn | head -n1 | awk '{print $2}')
    maxcov=$(echo "${tss_coverage}" | cut -f2 | uniq -c | sort -rn | head -n1 | awk '{print $1}')

    echo -e "${chrom}\t${start}\t${end}\t${maxpos}\t${maxcov}\t${name}" >> "${OUTPUT_FILE}"

  done < sorted_tss_regions_minus.bed
done

###################################################
# Pool strands per replicate
###################################################

for REP in ${REPLICATES}
do
  FINAL_OUTPUT="TSS_max_cov_and_pos_perregion_${REP}.txt"

  echo -e "chr\tstartsite-TSS\tendsite-TSS\tmaxcoverage_position\tmaxcoverage_coverage\tgenename-isoform" > "${FINAL_OUTPUT}"

  awk FNR-1 "plus_TSS_max_cov_and_pos_perregion_${REP}.txt" >> "${FINAL_OUTPUT}"
  awk FNR-1 "minus_TSS_max_cov_and_pos_perregion_${REP}.txt" >> "${FINAL_OUTPUT}"

done
