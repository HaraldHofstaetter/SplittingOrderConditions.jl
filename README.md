#SplittingOrderConditions.jl

A Julia package for generating order conditions for splitting methods.

For more information, see the examples below and the papers:

>[W. Auzinger](http://www.asc.tuwien.ac.at/~winfried), [W. Herfort](http://www.asc.tuwien.ac.at/~herfort/), [H. Hofstätter](http://www.harald-hofstaetter.at), [O. Koch](http://othmar-koch.org), [Setup of Order Conditions for Splitting Methods](http://arxiv.org/pdf/1605.00445.pdf), to appear in [Proceedings of CASC 2016](http://www.casc.cs.uni-bonn.de/2016/).

>[W. Auzinger](http://www.asc.tuwien.ac.at/~winfried), [H. Hofstätter](http://www.harald-hofstaetter.at), [D. Ketcheson](http://www.kaust.edu.sa/faculty/ketcheson.html), [O. Koch](http://othmar-koch.org), [Practical Splitting Methods for the Adaptive Integration of Nonlinear Evolution Equations. Part I: Construction of Optimized Schemes and Pairs of Schemes](http://link.springer.com/content/pdf/10.1007%2Fs10543-016-0626-9.pdf), [ BIT Numer. Math. (2016), doi:10.1007/s10543-016-0626-9](http://link.springer.com/article/10.1007/s10543-016-0626-9).


##Algorithm
The function [`generate_equations(q,s)`](https://github.com/HaraldHofstaetter/SplittingOrderConditions.jl/blob/master/src/SplittingOrderConditions.jl#L64) of this package implements the following algorithm:
>![](https://raw.githubusercontent.com/HaraldHofstaetter/SplittingOrderConditions.jl/master/generate_equations1.png)

##Installation
In a Julia notebook type
```julia
Pkg.clone("https://github.com/HaraldHofstaetter/SplittingOrderConditions.jl")
```
##Examples
To get easy access to the examples, copy them into the home directory:
```julia
cp(joinpath(homedir(), ".julia/v0.4/SplittingOrderConditions/examples/"), joinpath(homedir(), "SplittingOrderConditions_examples"), remove_destination=true)
```
Then 'SplittingOrderConditions_examples' will be listed in the JuliaBox home screen. The examples contain among others
+ [TestSplittingOrderConditions.ipynb](https://github.com/HaraldHofstaetter/SplittingOrderConditions.jl/blob/master/examples/TestSplittingOrderConditions.ipynb)
