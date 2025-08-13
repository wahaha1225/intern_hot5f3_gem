#!/bin/bash
#$ -pe smp 16
#$ -cwd
#$ -V

# load miniconda
module load miniconda

conda activate /projectnb/hfsp/marine_models_summer25/Yu_intern/EggNOG_CarveMe/envs/eggnog_env

emapper.py -i ../refs/Marinovum_algicola_DG898_assembly.faa \
--itype proteins \
--output /projectnb/hfsp/marine_models_summer25/Yu_intern/EggNOG_CarveMe/results/eggnog_algicola_output \
--cpu 16