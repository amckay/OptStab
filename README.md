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

## Figure 4

```
cd extension_with_saving/FullModel
matlab Main.m
cd extension_with_saving/NoSavings
matlab Main.m
cd extension_with_saving
matlab PlotResults.m
```
