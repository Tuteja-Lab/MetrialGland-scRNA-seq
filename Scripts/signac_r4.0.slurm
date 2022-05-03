#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename

#SBATCH --time=48:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --mail-user=hhvu@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=128G
#SBATCH --output="./reports/deg2timepoints-%j.out" # job standard output file (%j replaced by job id)
#SBATCH --error="./reports/deg2timepoints-%j.errors" # job standard error file (%j replaced by job id)



module load gcc/7.3.0-xegsmw4
module load r/4.0.2-py3-icvulwq
module load gsl/2.5-fpqcpxf
module load udunits/2.2.24-yldmp4h
module load gdal/2.4.4-nw2drgf
module load geos/3.8.1-2m7gav4

#Rscript 3_cluster.R
#Rscript 3_cluster_19.5.R
#Rscript 4_DEG.R
#Rscript 4_DEG_19.5.R
Rscript 5_DEG-2timepoints.R
#Rscript 3_cluster-15.5_rmChrMT.R
#Rscript 3_cluster-19.5_rmChrMT.R

