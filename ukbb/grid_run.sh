#!/bin/bash
#$ -N ddn_ukbb
#$ -cwd
#$ -l h_vmem=32G
#$ -l h_rt=2:00:00

module load Anaconda/3.7
python generate_disease_disease_edges.py

