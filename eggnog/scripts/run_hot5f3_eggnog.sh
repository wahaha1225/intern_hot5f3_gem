#!/bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V

# load miniconda
module load miniconda
# activate eggnog environment 
conda activate /projectnb/hfsp/marine_models_summer25/Yu_intern/EggNOG_CarveMe/envs/eggnog_env

emapper.py -i ../refs/NDPXLJ_6_hot5f3.faa \
    --itype proteins \
    --output /projectnb/hfsp/marine_models_summer25/Yu_intern/EggNOG_CarveMe/results/eggnog_output \
    --cpu 16 \
    --override