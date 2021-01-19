# Symplectic Tao's Method

A direct implementation of a symplectic integration method for generic non-separable hamiltonian 
systems of any degrees of freedom, as proposed by Molei Tao.

# Usage and tests

The files `ode_sys_solv.cpp` and `ode_sys_solv.hpp` is all you need to include into your code in 
order to use the method, with no extra compilation flags. An minimal example usage is given in 
`example/` folder. 

A series of performance tests were made for a 2D lattice hamiltonian model and its code for reproduction 
along with results can be found at `test/` folder.

# Credits

The method implemented here was originally developed by Molei Tao in: 
 
`M. Tao, "Explicit symplectic approximation of nonseparable Hamiltonians: algorithm and long time 
performance", ArXiv, 7 Sep 2016.`

and was implemented here without any modification besides small optmization flags.
