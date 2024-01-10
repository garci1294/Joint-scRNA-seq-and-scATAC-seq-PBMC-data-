#!/bin/bash -l
#SBATCH --time=16:00:00
#SBATCH --ntasks=16
#SBATCH --mem=64g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=garci624@umn.edu

cd /scratch.global/garci624/Single_Cell_Multiome

export PATH=/home/aventeic/garci624/local/src/lsoft/10xGen/cellranger-arc-1.0.1:$PATH

cellranger-arc count --id=10xMultiomic \
--reference=/home/aventeic/garci624/local/src/ref/refdata-cellranger-arc-GRCh38-2020-A-1.0.0 \
--libraries=/scratch.global/garci624/Single_Cell_Multiome/data/pbmc_granulocyte_sorted_10k_library.csv \
--localcores=16 \
--localmem=64