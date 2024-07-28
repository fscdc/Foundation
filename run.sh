#!/bin/bash

# output_dir
OUTPUT_DIR="/home/foundation/program/foundation-new/record/txt/"
mkdir -p "${OUTPUT_DIR}"

export CUDA_VISIBLE_DEVICES="1"

echo "CUDA_VISIBLE_DEVICES set to: $CUDA_VISIBLE_DEVICES"

GENE_LISTS=("mt_genes" "all_genes")  
EMBEDDING_METHODS=("scVI")  

declare -A DATASET_BATCHES
# DATASET_BATCHES["GSE206785"]="Patient,Tissue,Platform"
# DATASET_BATCHES["GSE206785_tumor"]="Patient,Tissue,Platform"
#DATASET_BATCHES["GSE261157"]="Sample"
#DATASET_BATCHES["kidney"]="tissue_general"
DATASET_BATCHES["pancreas"]="tissue_general"
# DATASET_BATCHES["Brain-visual-cortex"]="donor_id,Source"
# DATASET_BATCHES["Brain-excitatory-neurons"]="donor_id"
# DATASET_BATCHES["Brain-entorhinal-maturation"]="donor_id"
# DATASET_BATCHES["Brain-snRNA-seq"]="batch,donor_id"

CLUSTER_METHODS=("louvain")  
RESOLUTIONS=(0.3) 

# iterate over all combinations of parameters
for GENE_LIST in "${GENE_LISTS[@]}"; do
  for EMBEDDING_METHOD in "${EMBEDDING_METHODS[@]}"; do
    for DATASET in "${!DATASET_BATCHES[@]}"; do
      BATCH=${DATASET_BATCHES[$DATASET]}
      for CLUSTER_METHOD in "${CLUSTER_METHODS[@]}"; do
        for RESOLUTION in "${RESOLUTIONS[@]}"; do
          # output file name
          OUTPUT_FILE="${OUTPUT_DIR}/${GENE_LIST}-${EMBEDDING_METHOD}-${DATASET}-${CLUSTER_METHOD}-${RESOLUTION}.txt"

          # Start time
          START_TIME=$(date +%s)


          # run the pipeline       
          python -u main.py \
            --gene_list "$GENE_LIST" \
            --embedding_method "$EMBEDDING_METHOD" \
            --dataset "$DATASET" \
            --batch_list "$BATCH" \
            --cluster_method "$CLUSTER_METHOD" \
            --resolution "$RESOLUTION" > "$OUTPUT_FILE"


          # End time and duration
          END_TIME=$(date +%s)
          DURATION=$((END_TIME - START_TIME))
          echo "Execution time : $DURATION seconds"
        done
      done
    done
  done
done

echo "Finish!---by Sicheng Feng"
