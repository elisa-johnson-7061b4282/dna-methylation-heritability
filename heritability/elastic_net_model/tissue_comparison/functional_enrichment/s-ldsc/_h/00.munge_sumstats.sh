#!/bin/bash
#SBATCH --account=p32505        # Replace with your allocation
#SBATCH --partition=short       # Partition (queue) name
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=1     # Number of cores (CPU)
#SBATCH --mem=25G                # Memory limit
#SBATCH --job-name=munge_sumstats  # Job name
#SBATCH --output=logs/munge_sumstats/output_%j.log  # Standard output log
#SBATCH --error=logs/munge_sumstats/error_%j.log    # Standard error log

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

python munge_sumstats.py \
    --sumstats /projects/b1213/resources/public_data/gwas/alz/AD_sumstats_Jansenetal_2019sept.txt.gz \
    --merge-alleles /projects/p32505/users/elisa/dna-methylation-heritability/heritability/elastic_net_model/tissue_comparison/functional_enrichment/s-ldsc/w_hm3.snplist.gz \
    --out alz \
    --a1-inc

log_message "**** Job ends ****"