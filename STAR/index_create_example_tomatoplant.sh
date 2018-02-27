#!/usr/bin/env bash
# Leave only one comment symbol on selected options
# Those with two commets will be ignored:
# The name to show in queue lists for this job:
##SBATCH -J GENERATED_SCRIPT_FILE

# Number of desired cpus (can be in any node):
#SBATCH --ntasks=1

# Number of desired cpus (all in same node):
##SBATCH --cpus=48

# Amount of RAM needed for this job:
#SBATCH --mem=500gb

# The time the job will be running:
#SBATCH --time=168:00:00

# To use GPUs you have to request them:
##SBATCH --gres=gpu:1

# If you need nodes with special features uncomment the desired constraint line:
# * to request only the machines with 80 cores and 2TB of RAM
##SBATCH --constraint=bigmem
# * to request only machines with 16 cores and 64GB with InfiniBand network
##SBATCH --constraint=cal
# * to request only machines with 24 cores and Gigabit network
##SBATCH --constraint=slim

# Set output and error files
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

# MAKE AN ARRAY JOB, SLURM_ARRAYID will take values from 1 to 100
##SARRAY --range=1-100

# To load some software (you can show the list with 'module avail'):

module load star/2.5.1b


# the program to execute with its parameters:

time STAR --runMode genomeGenerate --runThreadN 48 --sjdbOverhang 74  --genomeDir ./  --genomeFastaFiles ./Solanum_lycopersicum.SL2.50.31.dna.genome.fa --sjdbGTFfile ./Solanum_lycopersicum.SL2.50.31.gtf $


