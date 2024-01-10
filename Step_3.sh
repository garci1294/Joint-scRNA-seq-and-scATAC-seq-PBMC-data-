#!/bin/bash -l
#SBATCH --time=04:00:00
#SBATCH --ntasks=4
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=garci624@umn.edu

cd /scratch.global/garci624/Single_Cell_Multiome/data

# sequencing data (FASTQ - 95.6 GB)
curl -O https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_fastqs.tar

# untar pbmc_granulocyte_sorted_10k_fastqs
tar -xvf pbmc_granulocyte_sorted_10k_fastqs.tar

# library (CSV - 220 B	)
curl -O https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_library.csv


