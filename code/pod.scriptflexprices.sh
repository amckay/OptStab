#!/bin/bash

#PBS -q B30
#PBS -N OptStabFlex
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=02:00:00
#PBS -j oe
#PBS -o logs/
#PBS -m e
#PBS -M alisdair.mckay@gmail.com

module load anaconda/5.0.1/python2.7
cd /home/amckay/OptStab/code/global
python OptStab_FlexiblePrices.py
