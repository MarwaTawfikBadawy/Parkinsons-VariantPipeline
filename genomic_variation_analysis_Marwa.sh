#!/bin/bash

SCRIPT="Genomic Variation Analysis - v1.0 (Dec 21, 2023)"

## Command line options
## -----------------------------------------------------------------------------
while getopts i:g:o:m:r:t: OPTION
do
case "${OPTION}"
in

# Input directory
i) INPUT=${OPTARG};;
# Reference genome
g) REF=${OPTARG};;
#Output directory
o) OUTPUT=${OPTARG};;
#Enable Merge VCF
m) MERGE=${OPTARG};;
# Rams
r) RAM=${OPTARG};;
# Threads
t) THREADS=${OPTARG};;

esac
done
## -----------------------------------------------------------------------------

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Input Directory  : ${INPUT}"
echo -e "[     INFO    ] Reference Genome  : ${REF}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Enable Merge VCFs : ${MERGE}"
echo -e "[     INFO    ] RAMs = ${RAM}"
echo -e "[     INFO    ] Threads = ${THREADS}"
## -----------------------------------------------------------------------------

## User confirmation
## -----------------------------------------------------------------------------
read -p "Continue (y/n)?" CHOICE
case "$CHOICE" in 
	y|Y ) 

## Main function
## -----------------------------------------------------------------------------                     
main () {

## Print pipeline info
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] ${SCRIPT}"
echo -e "[     INFO    ] ${AUTHOR}"
## -----------------------------------------------------------------------------  

## Print user args
## -----------------------------------------------------------------------------                     
echo -e "[     INFO    ] Input Directory  : ${INPUT}"
echo -e "[     INFO    ] Sample Name  : ${SAMPLE_NAME}"
echo -e "[     INFO    ] Reference Genome  : ${REF}"
echo -e "[     INFO    ] Output Directory : ${OUTPUT}"
echo -e "[     INFO    ] Enable Merge VCFs : ${MERGE}"
echo -e "[     INFO    ] RAMs = ${RAM}"
echo -e "[     INFO    ] Threads = ${THREADS}\n"
## -----------------------------------------------------------------------------

## Print start date/time
## -----------------------------------------------------------------------------                     
echo -e "[    START    ] $(date)\n"
## -----------------------------------------------------------------------------  


## Create output directory
## -----------------------------------------------------------------------------                     
echo -e "[   PROCESS   ] Creating output directory..."

# Check if output directory is exist
if [ -d "${OUTPUT}" ] 
then
    echo -e "[    ERROR    ] The output directory already exists.\n"
    exit 0
else
    mkdir -p ${OUTPUT}
fi

echo -e "[      OK     ] Output directory is ready on ${OUTPUT}\n"
## ----------------------------------------------------------------------------- 

## Zipping and Indexing Reference Genome (Run for first time)
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Zipping and Indexing Reference Genome..."

# Check if reference genome is already indexed
if [ -f "${REF}.gz" ] 
then
    echo -e "[    INFO    ] The reference genome is already indexed.\n"
else
    bgzip -c ${REF} > ${REF}.gz
    gatk CreateSequenceDictionary -R ${REF}.gz
    samtools faidx ${REF}.gz
fi

echo -e "[      OK     ] Indexing is done.\n"
## -----------------------------------------------------------------------------

## Convert GVCF to VCF using GATK
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Converting GVCF to VCF using GATK..."

VCF_DIR="${OUTPUT}/01_VCFs"
mkdir -p $VCF_DIR

# Get list of GVCFs
GVCFS=$(ls ${INPUT}/*.g.vcf.gz )
for GVCF in ${GVCFS}
do
    # Extract sample name
    SAMPLE_ID=$(basename "${GVCF}" | cut -d. -f1)

    echo -e "[   SUBPROCESS   ] Converting ${SAMPLE_ID} GVCF to VCF using GATK..."

    gatk IndexFeatureFile \
    -I ${GVCF}

    gatk --java-options "-Xmx${RAM}g" GenotypeGVCFs \
        -R ${REF}.gz \
        -V ${GVCF} \
        -O ${VCF_DIR}/${SAMPLE_ID}.vcf.gz \
        --allow-old-rms-mapping-quality-annotation-data

done

echo -e "[      OK     ] Conversion is done.\n"
## -----------------------------------------------------------------------------

## Apply variants quality control using BCFTools
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Applying variants quality control using BCFTools..."

FVCF_DIR="${OUTPUT}/02_Filtered_VCFs"
mkdir -p $FVCF_DIR

STATSQC_DIR="${OUTPUT}/03_QC_Stats"
mkdir -p $STATSQC_DIR

# Get list of VCFs
VCFS=$(ls ${VCF_DIR}/*.vcf.gz )
for VCF in ${VCFS}
do
    # Extract sample name
    SAMPLE_ID=$(basename "${VCF}" | cut -d. -f1)

    echo -e "[   SUBPROCESS   ] Apply variants QC on ${SAMPLE_ID} VCF using BCFTools..."

    bcftools filter \
        -i 'QUAL >= 30 || INFO/DP > 20' ${VCF} \
        -o ${FVCF_DIR}/${SAMPLE_ID}.filtered.vcf.gz \
        --output-type z \
        --threads ${THREADS}

    { echo "Before QC:"; bcftools query -f '%CHROM\t%TYPE\n' ${VCF} | awk '{count[$1"\t"$2]++} END {for (variant in count) print variant, count[variant]}' | sort -k1,1 -k3,3nr; echo ""; echo "After QC:"; bcftools query -f '%CHROM\t%TYPE\n' ${FVCF_DIR}/${SAMPLE_ID}.filtered.vcf.gz | awk '{count[$1"\t"$2]++} END {for (variant in count) print variant, count[variant]}' | sort -k1,1 -k3,3nr; } > ${STATSQC_DIR}/${SAMPLE_ID}.stats

done

echo -e "[      OK     ] QC is done.\n"
## -----------------------------------------------------------------------------


## Merge VCF files in one VCF file using BCFTools
## -----------------------------------------------------------------------------
# Check if merge enabled
if [[ ${MERGE} == "YES" || ${MERGE} == "yes" ]] 
then
    echo -e "[   PROCESS   ] Merging VCF files in one VCF file using BCFTools..."

    MVCF_DIR="${OUTPUT}/04_Merged_VCF"
    mkdir -p $MVCF_DIR

    # Merge command
    MERGE_CMD="bcftools merge --output-type z --threads ${THREADS}"

    #Get list of VCFs
    VCFS=$(ls ${FVCF_DIR}/*.filtered.vcf.gz )
    for VCF in ${VCFS}
    do
        
        # Index VCF file
        bcftools index -f --threads ${THREADS} ${VCF}

        # Add VCF path to merge command
        MERGE_CMD="${MERGE_CMD} ${VCF}"

    done

    # Add output vcf file
    MERGE_CMD="${MERGE_CMD} > ${MVCF_DIR}/all_samples.vcf.gz"

    # Run merge command
    eval " ${MERGE_CMD}"

    # Extract stats for merged VCF
    bcftools query -f '%CHROM\t%TYPE\n' ${MVCF_DIR}/all_samples.vcf.gz | awk '{count[$1"\t"$2]++} END {for (variant in count) print variant, count[variant]}' | sort -k1,1 -k3,3nr > ${MVCF_DIR}/all_samples.stats

    echo -e "[      OK     ] Merging is done.\n"
fi
## -----------------------------------------------------------------------------

## Perform variant annotation using VEP
## -----------------------------------------------------------------------------
echo -e "[   PROCESS   ] Performing variant annotation using VEP..."

ANN_DIR="${OUTPUT}/05_Annotated_Variants"
mkdir -p $ANN_DIR

ANNS_DIR="${OUTPUT}/06_Annotation_Stats"
mkdir -p $ANNS_DIR

# Get list of VCFs
VCFS=$(ls ${FVCF_DIR}/*.filtered.vcf.gz )
for VCF in ${VCFS}
do
    # Extract sample name
    SAMPLE_ID=$(basename "${VCF}" | cut -d. -f1)

    echo -e "[   SUBPROCESS   ] Performing variant annotation on ${SAMPLE_ID} VCF using VEP..."

    # vep --cache --everything \
    # -i ${VCF} \
    # -o ${ANN_DIR}/${SAMPLE_ID}.txt --tab \
    # --stats_file ${ANN_DIR}/${SAMPLE_ID}_summary.html\
    # --verbose \
    # --warning_file ${ANN_DIR}/${SAMPLE_ID}_warning.txt \
    # --skipped_variants_file ${ANN_DIR}/${SAMPLE_ID}_skipped_variants.txt \
    # --nearest symbol \
    # --allele_number \
    # --show_ref_allele \
    # --vcf_info_field ANN \
    # --fork ${THREADS}

    vep --cache --symbol --af_gnomadg --variant_class --gene_phenotype \
        -i ${VCF} \
        -o ${ANN_DIR}/${SAMPLE_ID}.tsv --tab \
        --stats_file ${ANN_DIR}/${SAMPLE_ID}_summary.html \
        --warning_file ${ANN_DIR}/${SAMPLE_ID}_warning.txt \
        --skipped_variants_file ${ANN_DIR}/${SAMPLE_ID}_skipped_variants.txt \
        --show_ref_allele \
        --verbose \
        --fork ${THREADS}

    # Count Genes
    echo -e "Count\tGenes" > ${ANNS_DIR}/${SAMPLE_ID}_Genes_Count.tsv
    cat ${ANN_DIR}/${SAMPLE_ID}.tsv | grep -v "^#" | cut -f20 | sort | uniq -c | sed 's/ \+/\t/g' | sed 's/^\t//g' >> ${ANNS_DIR}/${SAMPLE_ID}_Genes_Count.tsv

done

echo -e "[      OK     ] QC is done.\n"
## -----------------------------------------------------------------------------

## Print end date/time
## -----------------------------------------------------------------------------                     
echo -e "[     END     ] $(date)\n"
## -----------------------------------------------------------------------------  

} ## End of main function
## -----------------------------------------------------------------------------                     

## Prepare output log file
## -----------------------------------------------------------------------------                     
LOG=$( { time main > "$(dirname "${OUTPUT}")"/output.log 2>&1; } 2>&1 )
echo -e "Duration:${LOG}" >> "$(dirname "${OUTPUT}")"/output.log 2>&1
## -----------------------------------------------------------------------------                     

exit 0

;;

	n|N ) echo -e "[      OK     ] Process stopped.";;
	* ) echo -e   "[     ERROR   ] invalid";;
esac
