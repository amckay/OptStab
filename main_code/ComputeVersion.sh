#!/bin/bash
# script to compute results for one version

V=$1

if [ $V -eq 16 ]
then # if version 16 we use the grid from version 1
  cp ../data/results/GridBounds_v1.npz ../data/results/GridBounds_v16.npz
  JOB_Optim=$(qsub pod.scriptoptimize.sh  -F "ver $V")
  qsub -l depend=afterok:$JOB_Optim pod.scriptunpackserial.sh -F "ver $V param b mode point"

else # not version 16
  # find the grid
  python OptStab_find_grid.py ver $V

  # grid search
  python OptStab_search.py ver $V GridSearch NumProc 8

  # optimize over policy
  python OptStab_search.py ver $V

  # apply the propositions to unpack results  (used for Table 2)
  python OptStab_Unpack.py ver $V param b mode point

  # if version 1, we apply propositions to make the figures
  if [ $V -eq 1 ]
  then
    python OptStab_Unpack.py ver $V param b mode fig
    python OptStab_Unpack.py ver $V param tau mode fig
  fi



fi
