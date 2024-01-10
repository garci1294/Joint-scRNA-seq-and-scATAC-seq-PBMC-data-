#!/bin/bash -l
#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --mem=16g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=garci624@umn.edu

cd /home/aventeic/garci624/local/src/lsoft/10xGen

# download and install cellranger arc 1.0.1 from 10x genomics
curl -o cellranger-arc-1.0.1.tar.gz .... [KEY]

tar -xzvf cellranger-arc-1.0.1.tar.gz






