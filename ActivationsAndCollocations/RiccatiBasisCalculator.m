%-------------------------------------------------------------------------%
function [H, Hd, tau] = RiccatiBasisCalculator(varargin)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/15/2021


Description:

This function is used to solve the finite-horizon continuout-time linear
quadratic regulator problem.


Inputs:

struct - This stucture is used to store the parameters used to calculate
         basis matrix.

  .m - Amount of basis functions desired.
  .n - Amount of collocation points.
  .freeFunc - The selected free function.
  .collocType - The type of collocation points.
  .w - weight values for elm free functions.
  .b - bias values for elm free functions.


Outputs:

H  - Basis matrix where columns corresponds to the degree and rows
     corresponds to the time points.

Hd - First derivative of the basis matrix.

tau - Collocation points.

Version:

11/15/2021 - Initial completion.
12/15/2021 - Added lagrange polynomials.

%}

struct = varargin{1};

m = struct.m;
n = struct.n;
freeFunc = struct.freeFunc;
collocType = struct.collocType;
tau0 = struct.tau0;
tauf = struct.tauf;

switch collocType

    case 'chebyshev_gauss_lobatto'

        [tau, ~] = CGL_NodesWeights(n);

    case 'chebyshev_gauss_radau'

        [tau, ~] = CGR_NodesWeights(n);

    case 'chebyshev_gauss'

        [tau, ~] = CG_NodesWeights(n);

    case 'legendre_gauss_lobatto'

        [tau, ~] = LGL_NodesWeights_NewtonMethod(n);

    case 'legendre_gauss_radau'

        [tau, ~] = LGR_NodesWeights_NewtonMethod(n);

    case 'legendre_gauss'

        [tau, ~] = LG_NodesWeights_NewtonMethod(n);

    case 'equidistant'

        tau = linspace(tau0, tauf, n)';

end
    
if nargin == 2
    
    tauCo = tau;
    tau = varargin{2};
    
end

switch freeFunc
    
    case 'lagrange_orthogonal_polynomial'
        
        switch collocType
            
            case 'legendre_gauss_lobatto'
            
                if nargin == 1
                    H = LagrangePolynomials(tau, tau);
                    [Hd, ~, ~] = LGL_DiffMatrixA(n);
                    %Hd = LGL_DiffMatrixB(tau, tau);
                else
                    H = LagrangePolynomials(tau, tauCo)';
                    Hd = LGL_DiffMatrixB(tau, tauCo);
                end
                
            case 'chebyshev_gauss_lobatto'
                
                if nargin == 1
                    H = LagrangePolynomials(tau, tau);
                    [Hd, ~, ~] = CGL_DiffMatrixA(n);
                    %Hd = CGL_DiffMatrixB(tau, tau);
                else
                    H = LagrangePolynomials(tau, tauCo)';
                    Hd = CGL_DiffMatrixB(tau, tauCo);
                end
                
        end
        
    case 'chebyshev_orthogonal_polynomial'
        
        [T, Td, ~, ~, ~] = CP(tau, m+1);
        H = T(:,2:m+1);
        Hd = Td(:,2:m+1);
        
    case 'legendre_orthogonal_polynomial'
        
        [T, Td, ~, ~, ~] = LeP(tau, m+1);
        H = T(:,2:m+1);
        Hd = Td(:,2:m+1);
        
    case 'logistic_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'logistic');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'tanh_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'tanh');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'sin_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'sin');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'cos_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'cos');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'gaussian_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'gaussian');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'arctan_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'arctan');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'sinh_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'sinh');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'soft_plus_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'soft_plus');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'bent_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'bent');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'inv_sinh_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'inv_sinh');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'soft_sign_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'soft_sign');
        H = T(:,1:m);
        Hd = Td(:,1:m);
        
    case 'elu_elm'
        
        w = struct.w;
        b = struct.b;
        [T, Td, ~, ~, ~] = Activation(tau, w, b, 'soft_sign');
        H = T(:,1:m);
        Hd = Td(:,1:m);
    
end