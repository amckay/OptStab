#!/bin/bash
# Script to compute all the different versions


for version in 1 2 3 4 10 12 13 15 16 22 32
do
  source ComputeVersion.sh $version
done


python OptStab_FlexiblePrices.py
source ComputeCyclicalBenefit.sh


python Table2.py  # writes Table 2 to data/results/MacroStabTable.tex
jupyter-lab Figures123.ipynb # plots Figures 1, 2, and 3 in Jupyter notebook
