#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=garci624@umn.edu

cd /home/aventeic/garci624/local/src/ref

# download and install cellranger arc GRCh38-2020 reference from 10x genomics
curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A.tar.gz

tar -xzvf refdata-cellranger-arc-GRCh38-2020-A.tar.gz




