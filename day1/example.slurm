#!/bin/bash 
#SBATCH --partition=gc # partition to submit to
#SBATCH --reservation=workshop # this is the workshop reservation
#SBATCH --job-name=“myjob” # Job name
#SBATCH --nodes=1 # single node, anything more than 1 will not run
#SBATCH --ntasks=1 # equivalent to cpus, stick to around 20 max on gc64, or gc128 nodes
#SBATCH --mem=1000 # in MB, memory pool all cores, default is 2GB per cpu
#SBATCH --time=0-01:30:00  # expected time of completion in hours, minutes, seconds, default 1-day
#SBATCH --output=my_output.out # STDOUT
#SBATCH --error=my_output.err # STDERR
#SBATCH --mail-user=you@gmail.com
#SBATCH --mail-type=ALL

# This will be run once for a single process 
/bin/hostname
start=`date +%s`

module load blat/v.35-fastq

echo “HELLO”

end=`date +%s`
runtime=$((end-start))
echo $runtime


