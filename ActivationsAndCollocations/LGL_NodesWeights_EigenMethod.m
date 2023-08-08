%-------------------------------------------------------------------------%
function [tau, w] = LGL_NodesWeights_EigenMethod(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 09/24/2021

Description:  

Inputs:

n - Number of desired points.

Outputs:

tau - Legendre-Gauss-Lobatto points.
w   - Legendre-Gauss-Lobatto weights.


Version:

09/24/2021 - Initial completion.

%}
%z=leglndm(n) returns n Legendre-Gauss-Lobatto quadrature nodes with
%   z(1)=-1, z(n)=1 computed by the eigen-method
%   Recall that the interior nodes are zeros of L_{N-1}'(x)=c J_{n-2}^{1,1}(x) 
%   See Page 84 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%   Algorithms, Analysis and Applications, Springer Series in Compuational
%   Mathematics, 41, Springer, 2011.     
%
%Last modified on August 30, 2011       

if n <= 1
    disp('n should be bigger than 1'); 
    tau=[]; 
    return
end

if n == 2
    tau = [-1;1]; 
    return
end

if n == 3
    tau = [-1;0;1];
    return
end

av = zeros(1,n-2); 
j = [1:n-3]'; bv=j.*(j+2)./((2*j+1).*(2*j+3));
A = diag(av)+diag(sqrt(bv),1)+diag(sqrt(bv),-1);  % form the Jacobi matrix (3.142) with alpha=beta=1
ze = sort(eig(sparse(A)))';     % find the eigenvalues

tau = [-1, ze, 1]';

nn = n-1;
[~,y] = LegendrePoly(nn, ze);

% Use the weight expression (3.180) to compute the weights
w = [2/(nn*(nn+1)), 2./(nn*(nn+1)*y.^2), 2/(nn*(nn+1))]';

 return