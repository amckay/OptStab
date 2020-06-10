# "Optimal Automatic Stabilizers" Replication Material

Codes to replicate results in "Optimal Automatic Stabilizers" by Alisdair McKay and Ricardo Reis.

## Setting up the computing environment

Most computations were conducted using Python 2.7.  The extension with savings was conducted with Matlab.  To set up the Python environment using the Anaconda virtual environment manager:
```
conda create -n py27 python=2.7 numpy=1.15.4 matplotlib ipykernel ad scipy jinja2 pandas fredapi statsmodels
conda activate py27 
python -m ipykernel install --user
```

## Model versions

In the paper, we consider several variants of the model corresponding to different parameter values and model specifications.  In the code, these versions are referenced by the variable `VERSION`, which is set by a command-line argument when the scripts are called.


## The main results

The main results in our paper are the effect of business cycles on the optimal social insurance system reported in sections 6.2 and 6.3.  These results correspond to `VERSION = 1`.

First we need to find the grid on which we solve the model:
```
cd main_code
python OptStab_find_grid.py ver 1
```
Next we evaluate the objective function on a grid over the policy parameters, which is used for plotting purposes:
```
python OptStab_search.py ver 1 GridSearch NumProc 8
```
The option ```NumProc 8`` tells the program to use 8 processor cores in parallel to do the grid search. Finally we optimize over policy parameters:
```
python OptStab_search.py ver 1
```
Finally, to see the results and generate figures 1, 2, and 3:
```
jupyter-lab Figures123.ipynb
```

## All results (excluding Figure 4)

Run `ComputeAll.sh`

## Extension with savings (Figure 4)

```
cd extension_with_saving/FullModel
matlab Main.m
cd extension_with_saving/NoSavings
matlab Main.m
cd extension_with_saving
matlab PlotResults.m
```
