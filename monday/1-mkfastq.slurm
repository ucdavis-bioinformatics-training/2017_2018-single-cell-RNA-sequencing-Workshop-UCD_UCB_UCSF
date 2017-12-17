#!/bin/bash
## SLURM parameters
#SBATCH --time=5-0  # days, hours
#SBATCH --job-name=mkfastq  # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=16  # Number of cores
#SBATCH --mem=64000  # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=gc128  # Partition to submit to
#SBATCH --output=mkfastq.out  # File to which STDOUT will be written
#SBATCH --error=mkfastq.err  # File to which STDERR will be written
#SBATCH --mail-type=END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jnfass@ucdavis.edu  # Email to which notifications will be sent

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
basepath="."
bclpath="${basepath}/170721_J00113_0284_BHG27NBBXX_run470"
## https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
samplesheet="${basepath}/samplesheet.csv"
outfolder="${basepath}/cellranger-fastqs"

## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
## Create the call
call="cellranger mkfastq \
  --run=${bclpath} \
  --csv=${samplesheet} \
  --output-dir=${outfolder} \
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
