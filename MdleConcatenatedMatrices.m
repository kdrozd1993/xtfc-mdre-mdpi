%-------------------------------------------------------------------------%
function [mat3Struct, vecStruct] = MdleConcatenatedMatrices( ...
    Sf, mat3Struct)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 12/06/2021


Description:

This function generates 3D matrices that are necessary to transform the 
MRDE into a matrix differential lyapunov equation. The third dimension are 
the matrices at the time point. The function also generates the vectorized 
matrices.


Inputs:

Sf - The final condition of the MDRE.

mat3Struct - Structure that contains the 3D LTV matrices.


Outputs:

matStruct - Structure that contains the following 3D matrices.
  .A0  - MDLE state matrix.
  .A0T - MDLE state matrix transposed.
  .S0f - Final condition of the matrix differential lyapunov equation.
  .SfMinusGSm - (Sf - G3*Sm)
  .Sm - antistabilzing algebraic riccati equation solution at the end.

vecStruct - Structure of vetorized matrics similar to above. Each row
            represents a vectorized matrix at a time point. Matrices at
            each time point are vectorized with reshape(MAT', 1, nx*nx).
 
Version:

12/06/2021 - Initial completion.

%}

A3 = mat3Struct.A;
B3 = mat3Struct.B;
Q3 = mat3Struct.Q;
R3 = mat3Struct.R;
G3 = mat3Struct.G;
nx = size(A3, 1);
nu = size(B3, 2);
n = size(A3, 3);

A03 = zeros(nx, nx, n);
A0T3 = zeros(nx, nx, n);
for i = 1:n
    [~, Sp, ~] = lqr(-A3(:,:,i), B3(:,:,i), Q3(:,:,i), R3(:,:,i), zeros(nx, nu));   
    Sm3(:,:,i) = -Sp;
    A03(:,:,i) = A3(:,:,i) - G3(:,:,i)*Sm3(:,:,i);
    A0T3(:,:,i) = A03(:,:,i)';
    SfMinusSm3 = Sf - Sm3(:,:,i);
    S03 = pinv(SfMinusSm3);
end

mat3Struct.A0 = A03;
mat3Struct.A0T = A0T3;
mat3Struct.SfMinusSm = SfMinusSm3;
mat3Struct.S0 = S03;
mat3Struct.Sm = Sm3;

SmVec = zeros(1, nx*nx);
A0Vec = zeros(n, nx*nx);
A0TVec = zeros(n, nx*nx);
SfMinusSmVec = zeros(1, nx*nx);
S0Vec = zeros(1, nx*nx);
for i = 1:n
    SmVec(i,:) = reshape(Sm3(:,:,i)', 1, nx*nx);
    A0Vec(i,:) = reshape(A03(:,:,i)', 1, nx*nx);
    A0TVec(i,:) = reshape(A0T3(:,:,i)', 1, nx*nx);
    SfMinusSmVec(i,:) = reshape(SfMinusSm3', 1, nx*nx);
    S0Vec(i,:) = reshape(S03', 1, nx*nx);
end

vecStruct.A0 = A0Vec;
vecStruct.A0T = A0TVec;
vecStruct.SfMinusSm = SfMinusSmVec;
vecStruct.S0 = S0Vec;
vecStruct.Sm = SmVec;

end
