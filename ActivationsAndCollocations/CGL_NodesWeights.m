%-------------------------------------------------------------------------%
function [tau, w] = CGL_NodesWeights(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 10/18/2021

Description:

This functions returns Chebyshev-Gauss-Lobatoo points and arranged in 
ascending order. It also returns the Chebyshev-Gauss-Lobatto weights. See 
page 108 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods: 
Algorithms, Analysis and Applications, Springer Series in Compuational 
Mathematics, 41, Springer, 2011. 


Inputs:

n - Number of desired points.

Outputs:
tau - Chebyshev-Gauss points.
w   - Chebyshev-Gauss weights.

Version:

10/18/2021 - Initial completion.

%}

j = linspace(0, n, n)';
tau = -cos(pi*j/n);
w = [pi/(2*n); ones(n - 2, 1)*pi/n; pi/(2*n)];

end