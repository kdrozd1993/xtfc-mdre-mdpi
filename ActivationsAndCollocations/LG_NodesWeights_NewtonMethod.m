%-------------------------------------------------------------------------%
function [tau, w] = LG_NodesWeights_NewtonMethod(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 09/15/2021

Description:

This functions returns Legendre-Gauss points arranged in ascending order.
It also returns the Legnedre-Gauss weights. The method used to calculate
the points is the Newton iteration method. See Page 99 of the book: J. 
Shen, T. Tang and L. Wang, Spectral Methods: Algorithms, Analysis and 
Applications, Springer Series in Compuational Mathematics, 41, Springer, 
2011. Note that this function uses LegendrePoly().


Inputs:

n - Number of desired points.

Outputs:
tau - Legendre-Gauss points.
w   - Legendre-Gauss weights.

Version:

09/15/2021 - Initial completion.

%}

% Compute the initial guess of the interior LGL points
thetak = (4*(1:n)-1)*pi/(4*n+2);
ze=-(1-(n-1)/(8*n^3)-(39-28./sin(thetak).^2)/(384*n^4)).*cos(thetak);
ep = eps*10;                      % error tolerance for stopping iteration
ze1 = ze+ep+1;

% Newton's iteration procedure
while max(abs(ze1-ze)) >= ep            
    ze1 = ze;
    [dy,y] = LegendrePoly(n,ze);
    ze=ze-y./dy;  % see Page 99 of the book
end               % around 6 iterations are required for n=100

tau = ze'; % column vector

% Use the weight expression (3.178) to compute the weights
w = (2./((1-ze.^2).*dy.^2))';

end