# "Optimal Automatic Stabilizers" Replication Material

Codes to replicate results in "Optimal Automatic Stabilizers" by Alisdair McKay and Ricardo Reis.

## Setting up the computing environment

Most computations were conducted using Python 2.7 with plotting in JupyterLab.  The extension with savings was conducted with Matlab.  To set up the Python environment using the Anaconda virtual environment manager:
```
conda config --append channels conda-forge
conda config --append channels terhorst
conda create -n py27 python=2.7 numpy=1.15.4 matplotlib ipykernel ad  jinja2 pandas==0.20.3 fredapi scipy==0.19.1 statsmodels==0.8.0 jupyterlab
conda activate py27
python -m ipykernel install --user
```

## Model versions

In the paper, we consider several variants of the model corresponding to different parameter values and model specifications.  In the code, these versions are referenced by the variable `VERSION`, which is set by a command-line argument when the scripts are called.

## Data and empirical moments

Some of the calibration targets are empirical moments computed from data downloaded from the Federal Reserve Bank of St Louis FRED database.<sup>1</sup> The script
`main_code/EmpiricalMoments.py` reads the raw data series from csv files in `data/FRED` and saves the relevant moments to
`data/calibration/Moments.txt`.  The output file already exists in the replication package so you can skip this
step if you are not interested.


## To compute all results

The script `main_code/ComputeAll.sh` performs all calculations, produces figures and tables, and writes other results
to the screen.  You should run this script from within the main_code directory. Running this script takes many hours primarily because it solves the model with many different
parameter combinations to construct Table 2. Rather than running the script as it is written, you may prefer to look at the comments in the script to find the part that you are interested in.

## Explanation of steps to compute the main results

The main results in our paper are the effect of business cycles on the optimal social insurance system reported in sections 6.2 and 6.3.  These results correspond to `VERSION = 1`.  The steps below are executed by `ComputeAll.sh`, but we list them here to give some more explanation.

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
To apply the results from the propositions to unpack the numerical results:
```
python OptStab_Unpack.py ver 1 param b mode fig
python OptStab_Unpack.py ver 1 param tau mode fig
```
Finally, to see the results reported in the text and generate figures 1, 2, and 3:
```
jupyter-lab Figures123.ipynb
```

## References

<sup>1</sup> See [https://fred.stlouisfed.org/](https://fred.stlouisfed.org/). The following data series were accessed on Sep. 3, 2015:  
[GDP](https://fred.stlouisfed.org/series/GDP),
[GCE](https://fred.stlouisfed.org/series/GCE),
[unrate](https://fred.stlouisfed.org/series/unrate),
[UEMPLT5](https://fred.stlouisfed.org/series/UEMPLT5),
[UNEMPLOY](https://fred.stlouisfed.org/series/UNEMPLOY),
[PRS85006023](https://fred.stlouisfed.org/series/PRS85006023),
[EMRATIO](https://fred.stlouisfed.org/series/EMRATIO),
[JCXFE](https://fred.stlouisfed.org/series/JCXFE),
[FEDFUNDS](https://fred.stlouisfed.org/series/FEDFUNDS),
[CNP16OV](https://fred.stlouisfed.org/series/CNP16OV),
[DPI](https://fred.stlouisfed.org/series/DPI),
[W825RC1Q027SBEA](https://fred.stlouisfed.org/series/W825RC1Q027SBEA),
[CNP16OV](https://fred.stlouisfed.org/series/CNP16OV).
