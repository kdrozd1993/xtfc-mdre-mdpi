%-------------------------------------------------------------------------%
function [xi, CE, lossNorm, compTime] = TFC_MDLE_LS(lossFunc, setup)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/10/2021


Description:

This function performs least-squares in the XTFC framework for computing 
the Matrix Differential Lyapunov Equation.


Inputs:

lossFunc - String of a loss function.

setup - A structure of setup parameters.


Outputs:

xi - The final computed unkown values.

CE - The final computed constrained expression values.

lossNorm - The norm of the loss at each iteration.

compTime - An array structure that computes the jacobian computation and 
           inversion computation at each iteration.


Version:

11/10/2021 - Initial completion.

12/02/2021 - Adapted to handle the loss function just dealing with the 
             upper triangular elements of the riccati matrix instead of all 
             of them. The outputted constrained expressions and xi
             variables are for all of the elements though.

01/19/2022 - Added lsType.

%}

H = setup.tfcVariables.H;
S0fvutr = setup.S0fvutr;
nx = setup.tfcVariables.nx;
nCE = 1/2*nx*(nx+1);
n = setup.tfcVariables.n;
m = setup.tfcVariables.m;
lsType = setup.lsType;

[A, bb] = lossFunc(setup);
    
startInversionTime = tic;
xi = LS_Type(A, bb, lsType);
computeInversionTime = toc(startInversionTime);

compTime = struct('inversionComputationTime', computeInversionTime);

losses = reshape(A*xi - bb, n, nCE);

lossNorm = vecnorm(losses, 2, 2);

xi = reshape(xi, m, nCE);

CEut = (H - ones(n,1)*H(end,:))*xi + S0fvutr;

CE = zeros(n, nx*nx);
for j = 1:n
    V = CEut(j,:);
    ind = find(tril(ones(nx,nx), 0));
    CEn = zeros(nx);
    CEn(ind) = V;
    CEn = CEn'+tril(CEn, -1);
    CE(j,:) = reshape(CEn', 1, nx*nx);
end

xiut = xi;
xi = zeros(m, nx*nx);
for j = 1:m
    V = xiut(j,:);
    ind = find(tril(ones(nx,nx), 0));
    xin = zeros(nx);
    xin(ind) = V;
    xin = xin'+tril(xin, -1);
    xi(j,:) = reshape(xin', 1, nx*nx);
end

end