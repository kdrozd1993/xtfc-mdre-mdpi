%-------------------------------------------------------------------------%
function [P, P0] = LyapunovApproachRiccati(ltv3DMats, timeSpan, Pf)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/29/2021


Description:

This function performs the Lyapunov approach for computing the matrix
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
B3 = ltv3DMats.B;
Q3 = ltv3DMats.Q;
R3 = ltv3DMats.R;
G3 = ltv3DMats.G;

nT = size(A3, 3);
nx = size(A3, 1);
nu = size(B3, 2);

[~, Pp, ~] = lqr(-A3(:,:,end), B3(:,:,end), Q3(:,:,end), R3(:,:,end), zeros(nx, nu));   
Pm = -Pp;
P0f = pinv(Pf - Pm);

P0 = zeros(nx, nx, nT);
P0(:,:,end) = P0f;

P = zeros(nx, nx, nT);
P(:,:,end) = Pf;

for i = 1:nT-1
    
    stepT = timeSpan(end-i+1) - timeSpan(end-i);
    
    A0 = A3(:,:,end-i+1) - G3(:,:,end-i+1)*Pm;
    
    E = lyap(A0, -G3(:,:,end-i+1));
    
    P0(:,:,end-i) = expm(-A0*stepT)*(P0(:,:,end-i+1)-E)*expm(-A0*stepT)' + E; 
    
    P0(:,:,end-i) = 1/2*(P0(:,:,end-i)+P0(:,:,end-i)');
    
    P(:,:,end-i) = Pm + pinv(P0(:,:,end-i));
    
end