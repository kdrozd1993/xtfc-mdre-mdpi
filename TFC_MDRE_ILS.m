%-------------------------------------------------------------------------%
function [xi, iter, lossNorm, CE, compTime] = TFC_MDRE_ILS(xi, ...
    jacFunc, jacCase, maxIter, tol, setup)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/10/2021


Description:

This function performs iterative least-squares in the XTFC framework for 
computing the Riccati Differential Equation.


Inputs:

xi - A vector of orthogonal polynomial unkowns concatenated for each 
     element in the Riccati matrix.

jacFunc - A string that gives the name of the function that computes the
          jacobian.

jacCase - A string that specifies whether to calculate the jacobian with a
         'numeric' method or with the 'adigator' method.

maxIter - Maximum iterations of iterative least squares to perform.

tol - The tolerance that specifies when to stop iterative least squares.

setup - A structure of setup parameters.


Outputs:

xi - The final computed unkown values.

iter - The amount of iterations performed.

lossNorm - The norm of the loss at each iteration.

CE - The final computed constrained expression values.

compTime - An array structure that computes the jacobian computation and inversion
           computation at each iteration.


Version:

11/10/2021 - Initial completion.

12/02/2021 - Adapted to handle the loss function just dealing with the 
             upper triangular elements of the riccati matrix instead of all 
             of them. The outputted constrained expressions and xi
             variables are for all of the elements though.

01/19/2022 - Added lsType.

%}

err2 = [1e15, 1e15-1];
iter = 0;
nx = setup.tfcVariables.nx;
nCE = 1/2*nx*(nx+1);
n = setup.tfcVariables.n;
m = setup.tfcVariables.m;
H = setup.tfcVariables.H;
Sfvutr = setup.Sfvutr;

lsType = setup.lsType;

while abs(err2(2)) > tol && iter < maxIter && abs(err2(1) - err2(2)) > tol

    iter = iter + 1;

    err2(1) = err2(2);
    
    if strcmp(jacCase, 'numeric') == 1
        startJacobianTime = tic;
        J = numeric_jacobian(@TFC_MDRE_Loss_5, xi, setup);
        computeJacobianTime(iter) = toc(startJacobianTime);
        Lo = TFC_MDRE_Loss_5(xi, setup);
    elseif strcmp(jacCase, 'adigator') == 1
        startJacobianTime = tic;
        [J, Lo] = jacFunc(xi, setup);
        J = full(J);
        computeJacobianTime(iter) = toc(startJacobianTime);
    end

    startInversionTime = tic;
    [dxi] = LS_Type(J, Lo, lsType);
    computeInversionTime(iter) = toc(startInversionTime);
    
    losses = reshape(Lo, n, nCE);

    lossNorm = vecnorm(losses, 2, 2);

    dxi = reshape(dxi, m, nCE);
    
    xi = xi - dxi;

    err2(2) = norm(Lo);

end

compTime = struct('jacobianComputationTime', computeJacobianTime, ...
                  'inversionComputationTime', computeInversionTime);             
              
CEut = (H - ones(n,1)*H(end,:))*xi + Sfvutr;

CE = zeros(n, nx*nx);
for i = 1:n
    V = CEut(i,:);
    ind = find(tril(ones(nx,nx), 0));
    CEn = zeros(nx);
    CEn(ind) = V;
    CEn = CEn'+tril(CEn, -1);
    CE(i,:) = reshape(CEn', 1, nx*nx);
end

xiut = xi;
xi = zeros(m, nx*nx);
for i = 1:m
    V = xiut(i,:);
    ind = find(tril(ones(nx,nx), 0));
    xin = zeros(nx);
    xin(ind) = V;
    xin = xin'+tril(xin, -1);
    xi(i,:) = reshape(xin', 1, nx*nx);
end

end