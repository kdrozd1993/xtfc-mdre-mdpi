%-------------------------------------------------------------------------%
function [tau, w] = CG_NodesWeights(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 10/19/2021

Description:

This functions returns Chebyshev-Gauss points and arranged in 
ascending order. It also returns the Chebyshev-Gauss weights. See 
page 108 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods: 
Algorithms, Analysis and Applications, Springer Series in Compuational 
Mathematics, 41, Springer, 2011. 


Inputs:

n - Number of desired points.

Outputs:
tau - Chebyshev-Gauss points.
w   - Chebyshev-Gauss weights.

Version:

10/19/2021 - Initial completion.

%}

j = linspace(0, n, n)';
tau = -cos((2*j+1)*pi/(2*n+2));
w = pi/(n+1);

end