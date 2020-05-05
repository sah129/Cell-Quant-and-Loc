#!/bin/sh
#SBATCH --time=00:05:00
#SBATCH --mem-per-cpu=1024
#SBATCH --job-name=TestODonnell
#SBATCH --error=TestODonnell.stdout
#SBATCH --output=TestODonnell.stderr
 
module load r
R CMD BATCH PipelineTest.R