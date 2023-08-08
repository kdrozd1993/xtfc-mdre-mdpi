%-------------------------------------------------------------------------%
function [P] = KalmanEnglar(ltv3DMats, timeSpan, Pf)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/29/2021


Description:

This function performs the Kalman Englar Method for computing the matrix
differential riccati equation.


Inputs:

mat3Struct - Structure that contains the following 3D matrices.
  .A  - State matrix.
  .AT - State matrix transposed.
  .B  - Control matrix.
  .Q  - State difference matrix.
  .R  - Control effort matrix.
  .Ri - Inverse of control matrix.
  .G  - B*R*B'

timeSpan - Vector of time points.

Pf - Riccati matrix at the final time.


Outputs:

P - 3D riccati matrix where the third dimension represents the time points.

Version:

11/29/2021 - Initial completion.

%}

A3 = ltv3DMats.A;
AT3 = ltv3DMats.AT;
Q3 = ltv3DMats.Q;
G3 = ltv3DMats.G;

nT = size(A3, 3);
nx = size(A3, 1);

P = zeros(nx, nx, nT);
P(:,:,end) = Pf;

for i = 1:nT-1
    
    stepT = timeSpan(end-i+1) - timeSpan(end-i);
    
    Z = [A3(:,:,end-i+1), -G3(:,:,end-i+1); -Q3(:,:,end-i+1), -AT3(:,:,end-i+1)];
    
    transitionM = expm(-Z*stepT);
    
    transitionM11 = transitionM(1:nx,1:nx);
    transitionM12 = transitionM(1:nx,nx+1:nx*2);
    transitionM21 = transitionM(nx+1:nx*2,1:nx);
    transitionM22 = transitionM(nx+1:nx*2,nx+1:nx*2);
    
    P(:,:,end-i) = (transitionM21 + transitionM22*P(:,:,end-i+1))/ ...
        (transitionM11 + transitionM12*P(:,:,end-i+1));

    %P(:,:,end-i) = 0.5*(P(:,:,end-i) + P(:,:,end-i)');
    
end