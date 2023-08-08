%-------------------------------------------------------------------------%
function [tau, w] = LGL_NodesWeights_NewtonMethod(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 09/24/2021

Description:

This functions returns Legendre-Gauss-Lobatto points arranged in ascending 
order. It also returns the Legnedre-Gauss-Lobatto weights. The method used 
to calculate the points is the Newton iteration method. See Page 99 of the 
book: J. Shen, T. Tang and L. Wang, Spectral Methods: Algorithms, Analysis 
and Applications, Springer Series in Compuational Mathematics, 41, 
Springer, 2011. Note that this function uses LegendrePoly().


Inputs:

n - Number of desired points.

Outputs:
tau - Legendre-Gauss-Lobatto points.
w   - Legendre-Gauss-Lobatto weights.

Version:

09/24/2021 - Initial completion.

%}

% Compute the initial guess of the interior LGL points
nn = n-1;
thetak = (4*(1:nn)-1)*pi/(4*nn+2);
sigmak = -(1-(nn-1)/(8*nn^3)-(39-28./sin(thetak).^2)/(384*nn^4)).*cos(thetak);
ze = (sigmak(1:nn-1)+sigmak(2:nn))/2;
ep = eps*10;                      % error tolerance for stopping iteration
ze1 = ze+ep+1;

% Newton's iteration procedure
while max(abs(ze1-ze)) >= ep        
    ze1 = ze;
    [dy,y] = LegendrePoly(nn, ze);
    ze = ze-(1-ze.*ze).*dy./(2*ze.*dy-nn*(nn+1)*y); % see Page 99 of the book
end               % around 6 iterations are required for n=100

if exist('y', 'var') == 1
    [~,y] = LegendrePoly(nn, ze);
end

tau = [-1, ze, 1]'; % column vector

% Use the weight expression (3.180) to compute the weights
w = [2/(nn*(nn+1)),2./(nn*(nn+1)*y.^2),2/(nn*(nn+1))]';

end
