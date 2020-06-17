#!/bin/bash
# Script to compute all the different versions

# ----- Main results -----------
source ComputeVersion.sh 1
jupyter-lab Figures123.ipynb

# ------- Make Table 2 --------
for version in 2 3 4 10 12 13 15 16 22 32
do
  source ComputeVersion.sh $version
done
python OptStab_FlexiblePrices.py

python Table2.py  # writes Table 2 to data/results/MacroStabTable.tex

# -------- Robustness extensions ----------
python Figure5.py
source ComputeCyclicalBenefit.sh  # results with cyclical benefits

python Section_6_4_1.py  # report results with wage rule that responds to benefits and with zeta twice as high


# ------- Figure 4 (with savings) ----------
cd ../extension_with_saving/FullModel
matlab Main.m
cd ../NoSavings
matlab Main.m
cd ..
matlab PlotResults.m
