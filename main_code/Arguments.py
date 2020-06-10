# argument syntax
# key-values
# ver  [int], default 1
# CyclicalPolicy [int], default 0
# NumProc [int], default -1\
# param [str], default 'b'
# mode [str], default 'point'
#
# flags
# GridSearch, default false

import sys
A = [arg.lower() for arg in sys.argv]
n = len(A)


def find(key):
    place = [i for i in range(n) if A[i] == key]
    if len(place) == 0:
        return -1
    elif len(place) > 1:
        raise RuntimeError("multiple keys found in argument list")
    else:
        return place[0]

def identity(x):
    return x

def parse(key,default,transform = identity):
    i = find(key.lower())
    if i == -1:
        return default
    else:
        return transform(A[i+1])

def checkflag(flag):
    return find(flag.lower()) >0





class ArgObject:
    def __init__(self,ver, GridSearch, CyclicalPolicy, NumProc, param, mode):
        self.ver = ver
        self.GridSearch = GridSearch
        self.CyclicalPolicy = CyclicalPolicy
        self.NumProc = NumProc
        self.param = param
        self.mode = mode

    def __str__(self):
        return "------\nCommand-line arguments\nver = {0},\nGridSearch = {1},\nCyclicalPolicy = {2},\nNumProc = {3},\nparam = {4},\nmode = {5}\n------\n".format(self.ver, self.GridSearch, self.CyclicalPolicy, self.NumProc, self.param, self.mode)


CLArgs = ArgObject(
                    ver = parse("ver",1,lambda x: int(x)),
                    GridSearch = checkflag("GridSearch"),
                    CyclicalPolicy = parse("CyclicalPolicy",0,lambda x: int(x)),
                    NumProc = parse("NumProc",-1,lambda x: int(x)),
                    param = parse("param",'b'),
                    mode = parse("mode",'point')
                    )


print(CLArgs)
