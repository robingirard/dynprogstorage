from dynprogstorage.wrappers import GenCostFunctionFromMarketPrices
from numpy import random
import time

def test_dynprog():
    ### exemples simples d'utilisation de l'outil de programmation dynamique
    ###### storage operation example
    nbTime=250
    Prices=random.uniform(1, 1000, nbTime)
    p_max=1.
    c_max=10.*p_max

    ## x_i>0 : on stocke (on consomme du réseau)
    ## x_i<0 : on produit
    ### --> phi_i(x_i) est donc un coût Achat - vente que l'on veut minimiser
    cpl_func = GenCostFunctionFromMarketPrices(Prices.tolist())
    res = cpl_func.OptimMargInt([-p_max]*nbTime,[p_max]*nbTime,[0]*nbTime,[c_max]*nbTime)
    print(res)



def test_dyn_prog_perf():
    p_max = 1.
    c_max = 10. * p_max
    for nbTime in [1000, 10000, 100000, 10000000]:
        Prices = random.uniform(1, 1000, nbTime)
        start = time.time()
        cpl_func = GenCostFunctionFromMarketPrices(Prices.tolist())
        res = cpl_func.OptimMargInt([-p_max] * nbTime, [p_max] * nbTime, [0] * nbTime, [c_max] * nbTime)
        print('Elapsed time for ' + str(nbTime) + ' steps : ' + str(time.time() - start))

