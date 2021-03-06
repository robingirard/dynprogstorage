Dynamic programming optimisation tool
=======================================

This project implements the dynamic programming tool proposed in
[this paper](https://hal-mines-paristech.archives-ouvertes.fr/hal-01110689) (published [here](https://ieeexplore.ieee.org/document/6863551/))and
available in the R Software in [the ConConPiWiFun package](https://cran.r-project.org/web/packages/ConConPiWiFun/index.html)
---------------

### Table of Contents

* [1. Installation](#1.Installation)
* [2. Documentation](#2.documentation)
    * [2.1. Optimisation problem](#2.1.optimisation)
* [3. Examples](#3.examples)
    * [3.1. Simplest example](#3.1.simple)


## 1. Installation <a class="anchor" id="1.Installation"></a>

If you want to install the package from source, you can do

    pip install git+https://github.com/robingirard/dynprogstorage#egg=dynprogstorage

## Documentation <a class="anchor" id="2.documentation"></a> 

### optimisation problem <a class="anchor" id="2.1.optimisation"></a> 

This tool allows you to solve problems with the form 

    ## min_x  sum_i phi_i(x_i)   phi_i : convex piecewise linear function
    ## P_i^-<= x_i <=P_i^+
    ## C_i^-<= x_0 + sum_j=0^i x_j <= C_i^+ 



## Examples <a class="anchor" id="3.examples"></a> 
Let us give a few examples of use 
### storage operation example   <a class="anchor" id="3.1.simple"></a> 

While participating in the market with a 100% efficiency storage you want to maximize the profit 

    ## min_x  sum_i Pi_i x_i  ( phi_i linear function)
    ## -p_max <= x_i <=p_max
    ## 0<= x_0 + sum_j=0^i x_j <= c_max
    ## x_i>0 : consumption from the network
    ## x_i<0 : producing (injection to network)
    ### --> phi_i(x_i) is a buying cost we want to minimize

The code you need to use is  : 

    ## Definition of values    
    x_0=0
    nbTime=250
    Prices=random.uniform(1, 1000, nbTime)
    p_max=1.
    c_max=10.*p_max
    
    ## Generation of a vector of cost functions 
    cpl_func = GenCostFunctionFromMarketPrices(Prices.tolist())
    cpl_func.vec_get(0).getBreakPoints() ## what does the first cost function look like
    ## now solve the optimisation problem
    res = cpl_func.OptimMargInt([-p_max]*nbTime,[p_max]*nbTime,[-x_0]*nbTime,[c_max-x_0]*nbTime)
    print(res)
    
    ## Visualisation of results (power) with prices            
    period=100
    plt.plot(res[:100])
    plt.plot(-(Prices[:100]-Prices.mean())/Prices.max())
    plt.ylabel("Puissance (MW)")
    plt.xlabel("Index")
    plt.show()
    
    ## Visualisation of Energy evolution 
    energie=np.cumsum(res)
    plt.plot(energie[:100], color='g')
    plt.plot([0]*100, color='b')
    plt.plot([c_max]*100, color='b')
    plt.ylabel("Energie (MWh)")
    plt.xlabel("Index")
    plt.show()