#!/usr/bin/env python
# coding: utf-8

'''This script creates Table 2.  Results are written to data/results/MacroStabTable.tex '''


import jinja2
import numpy as np
import os
from OptStab import beta

OutDir = os.path.join(os.path.pardir,'data','results')


def FormatNumber(x,precision=3):
    if type(x) == str:
        return x
    elif type(x) == np.ndarray:
        return np.array([FormatNumber(xi) for xi in x])
    else:
        formatstring = '{x:.' + str(precision) + 'f}'
        return formatstring.format(x = x)



env = jinja2.Environment(
        "%<", ">%",
        "<<", ">>",
        "[§", "§]",
        loader=jinja2.FileSystemLoader(".")
        )


RowTemplate =   '(<<rownum>>) & <<rowname>>   & <<b>> & <<tau>> & <<std_dWdx>> & <<std_dxdb>> & <<cov_dWdx_dxdb>> \\\\'

template = env.from_string(RowTemplate)

def RepRate(b,tau):
    return 2 * b **(1.0/(1-tau)) -1

def GetRow(Version, rownum,rowname,Shocks=True):

    with np.load(os.path.join(OutDir,'SearchResults_v'+str(Version)+'.npz')) as X:
        if Shocks:
            b = X['OptCyc'][0]
            tau = X['OptCyc'][1]
        else:
            b = X['OptSS'][0]
            tau = X['OptSS'][1]

    if Shocks:
        with np.load(os.path.join(OutDir,'Unpack_v'+str(Version)+'_param_b_mode_point.npz')) as X:
            std_dWdx = X['std_dWdx'][0] *(1-beta)
            std_dxdb = X['std_dxdparam'][0] *(1-beta)
#             cov_dWdx_dxdb = X['Cov'][0].sum() *(1-beta)
            cov_dWdx_dxdb = X['Cov'][0,2:].sum() *(1-beta)
    else:
        std_dWdx = std_dxdb = cov_dWdx_dxdb = '--'


    return template.render(rownum = rownum,
                    rowname = rowname,
                    b = FormatNumber(b,3),
                    tau = FormatNumber(tau,3),
                    std_dWdx = FormatNumber(std_dWdx,3),
                    std_dxdb = FormatNumber(std_dxdb,3),
                    cov_dWdx_dxdb = FormatNumber(cov_dWdx_dxdb,3) )

def GetRowFlex(Version, rownum,rowname):

    with np.load(os.path.join(OutDir,'SearchResults_v'+str(Version)+'.npz')) as X:
            b = X['OptCyc'][0]
            tau = X['OptCyc'][1]


    return template.render(rownum = rownum,
                    rowname = rowname,
                    b = FormatNumber(b,3),
                    tau = FormatNumber(tau,3),
                    std_dWdx = '--',
                    std_dxdb = 'See note.',
                    cov_dWdx_dxdb = '--' )





fname = os.path.join(OutDir,'MacroStabTable.tex')
f = open(fname, "w")

f.write(os.linesep)
f.write(GetRow(1,'i','No aggr. shocks',Shocks=False))
f.write(os.linesep)
f.write(GetRow(1,'ii','Baseline'))
f.write(os.linesep)
f.write(GetRowFlex(8,'iii','Flexible prices'))
f.write(os.linesep)
f.write(GetRow(3,'iv','Aggressive monetary policy'))
f.write(os.linesep)
f.write(GetRow(22,'v','Smaller mon. shock'))
f.write(os.linesep)
f.write(GetRow(12,'vi','No $Q^u$'))
f.write(os.linesep)
f.write(GetRow(16,'vii','No skill risk'))
f.write(os.linesep)
f.write(GetRow(10,'viii','Acyclical G'))
f.write(os.linesep)
f.write(GetRow(32,'ix','Acyclical G, Higher G/Y ratio'))
f.write(os.linesep)
f.write('% Note: with flexible prices, $x_t$ is pinned down by the firm\'s first order condition and we find $|dx/db| < 3\times10^5$ across the state space.')


f.close()
