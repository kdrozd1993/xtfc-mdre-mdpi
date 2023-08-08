%-------------------------------------------------------------------------%
function [tau, w] = CGR_NodesWeights(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 10/18/2021

Description:

This functions returns Chebyshev-Gauss-Radau points with -1 at the start 
and arranged in ascending order. It also returns the Chebyshev-Gauss-Radau 
weights. See page 108 of the book: J. Shen, T. Tang and L. Wang, Spectral 
Methods: Algorithms, Analysis and Applications, Springer Series in 
Compuational Mathematics, 41, Springer, 2011. 


Inputs:

n - Number of desired points.

Outputs:
tau - Chebyshev-Gauss points.
w   - Chebyshev-Gauss weights.

Version:

10/18/2021 - Initial completion.

%}

j = linspace(0, n, n)';
tau = -cos((2*pi*j)/(2*n+1));
w = [pi/(2*n+1); ones(n - 1, 1)*2*pi/(2*n+1)];

end