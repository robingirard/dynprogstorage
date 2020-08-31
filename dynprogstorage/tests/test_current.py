#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 23 14:49:39 2020

@author: robin.girard
"""
import dynprogstorage
#from dynprogstorage import linearfunc
#cplfunction_obj = dynprogstorage.Pycplfunction(x0, y0, x1, y1)


#### test Pycplfunction
A=dynprogstorage.Pycplfunction([-1.,1.],[-10.,0.],0.)
A.getBreakPoints()
A.evalf(1)
A.evalf(0)
A.AddSimple(1.,2.,3.,4.)
A.getBreakPoints()

import math
B=dynprogstorage.Pycplfunction([-1,1],[-math.inf,0],0.)
B.getBreakPoints()

#### test Pycplfunctionvec
C= dynprogstorage.Pycplfunctionvec([-1.,1.],[-10.,0.],[0,0])
C.size()
C.Maxf_1Breaks_withO([-1.,3.],[-12.,2.],[0,0])
E= dynprogstorage.Pycplfunctionvec([-2.,2.],[-12.,2.],[0,0])
E=E.Max_(D,C)

import dynprogstorage
from numpy import random
import math


#### cas avec S1, B1, f0
n=100
S1=random.uniform(low=1,high=10,size=n).tolist()
E= dynprogstorage.Pycplfunctionvec(S1,[-math.inf]*n,[0]*n)
PP=1; CC=5*PP;
xEtoile=E.OptimMargInt([-PP]*n,[PP]*n,[0]*n,[CC]*n)
xEtoile

#### cas avec 2 breakpoints S1, B1, f0 , S2, B2
n=100
S1=random.uniform(low=1,high=10,size=n)
S2=S1+2
E= dynprogstorage.Pycplfunctionvec(S1.tolist(),[-math.inf]*n,[0]*n,#S1, B1, f0
                                   S2.tolist(),[0]*n) #S2, B2
PP=1; CC=5*PP;
xEtoile=E.OptimMargInt([-PP]*n,[PP]*n,[0]*n,[CC]*n)

from dynprogstorage import wrappers
from numpy import random

### exemples simples d'utilisation de l'outil de programmation dynamique
###### storage operation example
nbTime=250
Prices=random.uniform(1, 1000, nbTime)
p_max=1.;  c_max=10.*p_max;

CostFunction= wrappers.GenCostFunctionFromMarketPrices(Prices)
## x_i>0 : on stocke (on consomme du réseau)
## x_i<0 : on produit
### --> phi_i(x_i) est donc un coût Achat - vente que l'on veut minimiser
res=CostFunction.OptimMargInt([-p_max]*nbTime,[p_max]*nbTime,[0]*nbTime,[c_max]*nbTime)
## min sum_i phi_i(x_i)
## -p_max <= x_i <=  p_max forall i
## 0 <= sum_j=1^ix_j <= C_max  forall i


#plt.plot(res)
#plt.show()