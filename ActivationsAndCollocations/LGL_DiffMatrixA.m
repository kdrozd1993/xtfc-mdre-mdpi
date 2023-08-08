%-------------------------------------------------------------------------%
function [D, x, w] = LGL_DiffMatrixA(n)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 10/12/2021

Description:

This function generates the Legendre-Gauss-Lobatto collocation points and 
weights (x, w) as well as the the differential approximation matrix D. 
Note that these differentiation matrices are associated directly with the
Legendre psuedospectral method. This means that the differentiation matrix 
calculated for the LGL points. See chapter 11 in the book A Graduate 
Introduction to Numerical Methods by Corless and Fillion (2013) for more 
info on calculated the differentiation matrix for Lagrange orthogonal 
polynomials.

Inputs:

n - Number of LGL points.

Outputs:

x  - LGL points.
w  - LGL weights.
Dg - LGL differentiation matrix.

Version:

10/12/2021 - Initial completion.

%}

% Gauss Pts
[x, w] = LGL_NodesWeights_EigenMethod(n);

% Add initial point -1
x = [x];
x = sort(x);
n = n;

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

end