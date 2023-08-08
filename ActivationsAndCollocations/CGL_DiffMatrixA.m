%-------------------------------------------------------------------------%
function [D, x, w] = CGL_DiffMatrixA(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 10/19/2021

Description:

This function generates the flipped Chebyshev-Gauss-Lobatto collocation 
points and weights (x, w) as well as the the differential approximation 
matrix. See chapter 11 in the book A Graduate Introduction to Numerical 
Methods by Corless and Fillion (2013) for more info on calculating the
differentiation matrix for Lagrange orthogonal polynomials.

Inputs:

n - Number of Gauss points.

Outputs:

x    - Chebyshev-Gauss-Lobatto points.
w    - Chebyshev-Gauss-Lobatto weights.
Dg   - Chebyshev-Gauss-Lobatto differentiation matrix.

Version:

10/19/2021 - Initial completion.

%}

% Radau Pts
[x, w] = CGL_NodesWeights(n);

% Eval derivative of Lagrange polynomials
for j = 1:n
    
    for i = 1:n
        
        prod = 1;
        summ = 0;
        
        if j == i           
            for k = 1:n
                
                if k~=i   
                    summ = summ+1/(x(i)-x(k));
                end
                
            end
            
            D(i,j) = summ;            
        else           
            for k = 1:n
                
                if (k~=i)&&(k~=j)
                    
                    prod = prod * (x(i)-x(k));
                    
                end
                
            end
            
            for k = 1:n
                
                if k~=j
                    prod = prod/(x(j)-x(k));
                end
                
            end
            
            D(i,j) = prod;            
        end
        
    end
    
end
