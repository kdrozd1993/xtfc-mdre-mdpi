%-------------------------------------------------------------------------%
function [tau, w] = LGR_NodesWeights_NewtonMethod(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 09/30/2021

Description:

This functions returns Legendre-Gauss-Radau points with -1 at the start and 
arranged in ascending order. It also returns the Legnedre-Gauss-Radau 
weights. The method used to calculate the points is the Newton iteration 
method. See Page 99 of the book: J. Shen, T. Tang and L. Wang, Spectral 
Methods: Algorithms, Analysis and Applications, Springer Series in 
Compuational Mathematics, 41, Springer, 2011. Note that this function uses 
LegendrePoly().


Inputs:

n - Number of desired points.

Outputs:
tau - Legendre-Gauss-Radau points.
w   - Legendre-Gauss-Radau weights.

Version:

09/30/2021 - Initial completion.

%}

% Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
nn = n-1;
ze = -cos(2*pi*(0:nn)/(2*nn+1))';
ep = eps*10;                      % error tolerance for stopping iteration
ze1 = ze+ep+1;

% Newton's iteration procedure
i = 0;
while max(abs(ze1-ze)) >= ep            
    ze1 = ze;
    [dy1,y1] = LegendrePoly(nn,ze);
    [dy2,y2] = LegendrePoly(n,ze);
    ze = ze - (y1+y2)./(dy1+dy2);  % see Page 99 of the book
end               

tau = [ze];

[~,y] = LegendrePoly(nn, tau);

% Use the weight expression (3.180) to compute the weights
w = (1-tau)./(n^2*y.^2); 