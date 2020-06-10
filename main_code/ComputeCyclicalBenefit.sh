#!/bin/bash

# find the grid
python OptStab_find_grid.py ver 1 CyclicalPolicy 2

# grid search
python OptStab_search.py ver 1 GridSearch NumProc 8 CyclicalPolicy 2

# optimize over policy
python OptStab_search.py ver 1 CyclicalPolicy 2
