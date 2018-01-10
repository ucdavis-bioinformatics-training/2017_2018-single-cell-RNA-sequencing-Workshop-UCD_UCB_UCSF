#!/bin/bash
## SLURM parameters
#SBATCH --time=10-0  # days-minutes
#SBATCH --job-name=counts # Job name
#SBATCH --array=1-6  # lines in samples.txt
#SBATCH --nodes=1
#SBATCH --ntasks=32 # Number of cores
#SBATCH --mem=64000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=gc128 # Partition to submit to
#SBATCH --output=countArrayJob_%A_%a.out # File to which STDOUT will be written
#SBATCH --error=countArrayJob_%A_%a.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jnfass@ucdavis.edu # Email to which notifications will be sent

## Record the start time
start=`date +%s`
## Record the host being run on
hostname

## Where cellranger and bcl2fastq executables are located
## a) by loading a module
module load cellranger/2.0.1
module load bcl2fastq/2.17.1.14
## b) or, by placing the location of the executables on the path (edit to your location)
# export PATH=/share/biocore/software/bin/:$PATH
## c) or if they are already on the path, do nothing

## Set the number of theads and max memory,
##  either determined from Slurm or set yourself if running on command line
THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)
echo "Slurm job id: " $SLURM_JOB_ID
echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

## Set the parameters for the run
transcriptome="/share/genomes/cellranger_genomes/refdata-cellranger-mm10-1.2.0"
basepath="."
fastqs="${basepath}/cellranger-fastqs"
## in a slurm task array the index is set with the SLURM_ARRAY_TASK_ID variable
## if running on the command line, need to set to the row # of the sample to be run 
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt`

## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
## Create the call
call="cellranger count \
  --id=${sample} \
  --sample=${sample} \
  --transcriptome=${transcriptome} \
  --fastqs=${fastqs} \
  --localcores=${THREADS} \
  --localmem=${MEM}"

## Echo the call
echo $call
## Evaluate the call
#eval $call

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
