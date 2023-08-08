%-------------------------------------------------------------------------%
function [mat3Struct, concatStruct] = MdreConcatenatedMatricesFromFuncs( ...
    t, nx, nu, ltvMatFuncs)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/10/2021


Description:

This function generates 3D LTV matrices where the 3rd dimension represents
the time point. The function also generates 2D concatenated LTV matrices.


Inputs:

t - a vector where each element is a moment in time.

nx - scalar representing the length of the state vector.

nu - scalar representing the length of the control vector.

ltvMatFuncs - Structure of weight mat functions with respect to the 
              independent variable for linear time varying systems (can 
              be autonomous too).

  .AB  - Autonomous function for linear state & control matrices.

  .Q  - Autonomous function for state difference matrix.

  .R  - Autonomous function for control effort matrix.


Outputs:

mat3Struct - Structure that contains the following 3D matrices.
  .A  - State matrix.
  .AT - State matrix transposed.
  .B  - Control matrix.
  .Q  - State difference matrix.
  .R  - Control effort matrix.
  .Ri - Inverse of control matrix.
  .G  - B*R*B'

concatStruct - Structure of vetorized matrics similar to above. Each row
               represents a vectorized matrix at a time point. Matrices at
               each time point are vectorized with reshape(MAT', 1, nx*nx).
 
Version:

11/10/2021 - Initial completion.
12/02/2021 - Changed concatStruct to have vectorized matrices.

%}

n = length(t);

A3 = zeros(nx, nx, length(t));
AT3 = zeros(nx, nx, length(t));
B3 = zeros(nx, nu, length(t));
Q3 = zeros(nx, nx, length(t));
R3 = zeros(nu, nu, length(t));
Ri3 = zeros(nu, nu, length(t));
G3 = zeros(nx, nx, length(t));
for i = 1:n
    AB = ltvMatFuncs.AB(t(i)); 
    A3(:,:,i) = AB{1};
    B3(:,:,i) = AB{2};
    AT3(:,:,i) = A3(:,:,i)';
    Q3(:,:,i) = ltvMatFuncs.Q(t(i)); 
    R3(:,:,i) = ltvMatFuncs.R(t(i));  
    Ri3(:,:,i) = inv(R3(:,:,i));
    G3(:,:,i) = B3(:,:,i)*Ri3(:,:,i)*B3(:,:,i)';
end

mat3Struct.A = A3;
mat3Struct.AT = AT3;
mat3Struct.B = B3;
mat3Struct.Q = Q3;
mat3Struct.R = R3;
mat3Struct.Ri = Ri3;
mat3Struct.G = G3;

ACon = zeros(length(t), nx*nx);
ATCon = zeros(length(t), nx*nx);
BCon = zeros(length(t), nx*nu);
QCon = zeros(length(t), nx*nx);
RCon = zeros(length(t), nu*nu);
RiCon = zeros(length(t), nu*nu);
GCon = zeros(length(t), nx*nx); 

for i = 1:n
    ACon(i,:) = reshape(A3(:,:,i)', 1, nx*nx);
    ATCon(i,:) = reshape(AT3(:,:,i)', 1, nx*nx);
    BCon(i,:) = reshape(B3(:,:,i)', 1, nx*nu);
    QCon(i,:) = reshape(Q3(:,:,i)', 1, nx*nx);
    RCon(i,:) = reshape(R3(:,:,i)', 1, nu*nu);
    RiCon(i,:) = reshape(Ri3(:,:,i)', 1, nu*nu);
    GCon(i,:) = reshape(G3(:,:,i)', 1, nx*nx);
end

concatStruct.A = ACon;
concatStruct.AT = ATCon;
concatStruct.B = BCon;
concatStruct.Q = QCon;
concatStruct.R = RCon;
concatStruct.Ri = RiCon;
concatStruct.G = GCon;

end