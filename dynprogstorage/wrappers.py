# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Sat May 30 15:46:28 2020
#
# @author: robin.girard
# """
#
from dynprogstorage.Wrapper_dynprogstorage import Pycplfunctionvec
#from dynprogstorage import linearfunc
import numpy
import math


def pmin(x,y):
    n=len(x)
    TMP= numpy.concatenate((numpy.array(x),numpy.array(y))).reshape(2,n)
    return(numpy.amin(TMP,axis=0).tolist())

def pmax(x,y):
    n=len(x)
    TMP= numpy.concatenate((numpy.array(x),numpy.array(y))).reshape(2,n)
    return(numpy.amax(TMP,axis=0).tolist())

def GenCostFunctionFromMarketPrices(Prices,r_in=1,r_out=1,valueAtZero=0):
    #import dynprogstorage
    Prices=numpy.array(Prices)
    n=len(Prices)

    if ((type(valueAtZero)==int)|(type(valueAtZero)==float)) : valueAtZero=[valueAtZero]*n

    if (type(valueAtZero)==list) : valueAtZero = numpy.array(valueAtZero)

    if (r_in*r_out==1):
        Res= Pycplfunctionvec(Prices.tolist(),[-math.inf]*n,valueAtZero.tolist())
    else:
        Prices_in = Prices * r_out
        Prices_out = Prices / r_in
        Res= Pycplfunctionvec(Prices_in.tolist(),[-math.inf]*n,valueAtZero.tolist(),Prices_out.tolist(),[0]*n)

    return Res;

def GenCostFunctionFromMarketPrices_dict(Prices,r_in=1,r_out=1,valueAtZero=0):
    #import dynprogstorage
    Prices=numpy.array(Prices)
    n=len(Prices)

    if ((type(valueAtZero)==int)|(type(valueAtZero)==float)) : valueAtZero=[valueAtZero]*n

    if (type(valueAtZero)==list) : valueAtZero = numpy.array(valueAtZero)


    if (r_in*r_out==1):
        Res= {"S1" : Prices.tolist(),"B1" : [-math.inf]*n,"f0" :valueAtZero.tolist()}
    else:
        Prices_in = Prices * r_out
        Prices_out = Prices / r_in
        Res= {"S1" : Prices_in.tolist(),"B1" : [-math.inf]*n, "f0" : valueAtZero.tolist(), "S2" : Prices_out.tolist(),"B2" : [0]*n}

    return Res;