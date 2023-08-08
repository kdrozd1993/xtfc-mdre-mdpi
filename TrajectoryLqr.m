%-------------------------------------------------------------------------%
function z = TrajectoryLqr(timeSpan, ltvMatFuncs, Sf, tol)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/10/2021


Description:

This function is used to solve the finite-horizon continuout-time linear
quadratic regulator problem.


Inputs:

timeSpan - Monotonically increasing vector of times at which the riccati 
           matrix to be calculated.

ltvMatFuncs - Structure of anonymous functions for ltv A, B, Q, and R
              matrices

Sf - Final state cost matrix.

tol - Relative and absolute tolerance for MATLAB's Runge Kutta schemes.


Outputs:

z - A matrix representing the riccati matrix values. The riccati matrix 
    values are reshaped to a row vector, such that each row of the z matrix 
    represents the time points.


Version:

11/10/2021 - Initial completion.

%}


% INPUTS:
%
% OUTPUTS:
%   Soln = struct array with solution at each point in t
%   Soln(i).t = t(i);
%   Soln(i).K = gain matrix at t(i)
%   Soln(i).S = Cost to go
%   Soln(i).E = close-loop eigen-values for system at t(i)
%
% NOTES:
%
%   J = x'Fx + Integral {x'Qx + u'Ru} dt
%
% See Also LQR

% nSoln = length(t);
% Soln(nSoln).t = 0;
% Soln(nSoln).K = zeros(nState, nInput);
% Soln(nSoln).S = zeros(nState, nState);
% Soln(nSoln).E = zeros(nState, 1);
% 
% for idx=1:nSoln
%     i = nSoln-idx+1;
%     zNow = z(:,i);
%     tNow = t(i);
%     S = reshape(zNow,nState,nState);
%     [A,B] = linSys(tNow);
%     K = R\(B'*S);
%     Soln(i).t = tNow;
%     Soln(i).K = K;
%     Soln(i).S = S;
%     Soln(i).E = eig(A-B*K);
% end

nx = size(Sf, 1);

userFun = @(t, z) rhs(t, z, ltvMatFuncs.AB, ltvMatFuncs.Q(t), ltvMatFuncs.R(t), nx);
z0 = reshape(Sf, nx*nx, 1);
timeSpanInt = flip(timeSpan);
options = odeset();
options.RelTol = tol;
options.AbsTol = tol;
%sol = ode45(userFun, timeSpanInt, z0, options);
sol = ode4(userFun, timeSpanInt, z0);
%z = deval(sol, timeSpan)';
z = flip(sol);

end

function dz = rhs(t, z, linSys, Q, R, nState)
    P = reshape(z, nState, nState);
    [AB] = linSys(t);
    A = AB{1};
    B = AB{2};
    dP = ricatti(A, B, Q, R, P);
    dz = reshape(dP, nState*nState, 1);
end

function [dP, K] = ricatti(A, B, Q, R, P)
    K = R\B'*P;
    dP = -(A'*P + P*A - P*B*K + Q);
end


