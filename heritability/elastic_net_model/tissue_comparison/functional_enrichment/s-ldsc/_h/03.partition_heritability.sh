#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=partition_heritability  # Job name
#SBATCH --output=logs/partition_heritability/output_%j.log  # Standard output log
#SBATCH --error=logs/partition_heritability/error_%j.log    # Standard error log

# Log function
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

python /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/ldsc/ldsc.py \
  --h2 /projects/b1213/resources/public_data/gwas/alz/AD_sumstats_Jansenetal_2019sept.txt.gz \
  --ref-ld-chr baselineLD_v2.2/baselineLD.,brain_enhancer_annot/ \
  --w-ld-chr weights_hm3_no_hla/weights. \
  --overlap-annot \
  --frqfile-chr 1000G.EUR.QC. \
  --out results/brain_enhancer

log_message "**** Job ends ****"