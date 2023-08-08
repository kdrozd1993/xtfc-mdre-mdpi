%-------------------------------------------------------------------------%
function [loss] = TFC_MDRE_Loss_5(xi, setup)
%-------------------------------------------------------------------------%

% %
% 
% Authors: Kristofer Drozd
% Created: 12/02/2021
% 
% 
% Description:
% 
% This function generates the loss vector for the 5x5 Matrix Riccati 
% Differential Equation for the PINN/XTFC/TFC methodology.
% 
% 
% Inputs:
% 
% xi - A vector of orthogonal polynomial unkowns where columns represent 
%      the element in the Riccati matrix and rows are polynomial order.
% 
% setup - A structure representing inputs.
% 
% 
% Outputs:
% 
% loss - Loss vector.
% 
% 
% Version:
% 
% 12/02/2021 - Initial completion.
% 
% %

Q = setup.mdreVecMats.Q;
nx = sqrt(size(Q, 2));
A = setup.mdreVecMats.A;
AT = setup.mdreVecMats.AT;
G = setup.mdreVecMats.G;

n = setup.tfcVariables.n;
m = setup.tfcVariables.m;
t0 = setup.t0;
tf = setup.tf;
tau = setup.tau;

Sfvutr = setup.Sfvutr;

H = setup.tfcVariables.H;
Hd = setup.tfcVariables.Hd;

c = (tau(end) - tau(1)) / (tf - t0);

CE = (H - ones(n,1)*H(end,:))*xi + Sfvutr;
CEd = c*Hd*xi;

CeColInds = zeros(nx,nx);
CeColInds(1,:) = [1, 2, 3, 4, 5];
CeColInds(2,:) = [2, 6, 7, 8, 9];
CeColInds(3,:) = [3, 7, 10, 11, 12];
CeColInds(4,:) = [4, 8, 11, 13, 14];
CeColInds(5,:) = [5, 9, 12, 14, 15];

CeRowInds = zeros(nx,nx);
CeRowInds(1,:) = [1, 2, 3, 4, 5];
CeRowInds(2,:) = [2, 6, 7, 8, 9];
CeRowInds(3,:) = [3, 7, 10, 11, 12];
CeRowInds(4,:) = [4, 8, 11, 13, 14];
CeRowInds(5,:) = [5, 9, 12, 14, 15];

colInds = zeros(nx,nx);
colInds(1,:) = 1:1:5;
colInds(2,:) = 6:1:10;
colInds(3,:) = 11:1:15;
colInds(4,:) = 16:1:20;
colInds(5,:) = 21:1:25;

rowInds = zeros(nx,nx);
rowInds(1,:) = 1:5:25;
rowInds(2,:) = 2:5:25;
rowInds(3,:) = 3:5:25;
rowInds(4,:) = 4:5:25;
rowInds(5,:) = 5:5:25;

S11 = zeros(n, nx);
S11(:,1) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(1,1));
S11(:,2) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(1,2));
S11(:,3) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(1,3));
S11(:,4) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(1,4));
S11(:,5) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(1,5));

S12 = zeros(n, nx);
S12(:,1) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(2,1));
S12(:,2) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(2,2));
S12(:,3) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(2,3));
S12(:,4) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(2,4));
S12(:,5) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(2,5));

S13 = zeros(n, nx);
S13(:,1) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(3,1));
S13(:,2) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(3,2));
S13(:,3) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(3,3));
S13(:,4) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(3,4));
S13(:,5) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(3,5));

S14 = zeros(n, nx);
S14(:,1) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(4,1));
S14(:,2) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(4,2));
S14(:,3) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(4,3));
S14(:,4) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(4,4));
S14(:,5) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(4,5));

S15 = zeros(n, nx);
S15(:,1) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(5,1));
S15(:,2) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(5,2));
S15(:,3) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(5,3));
S15(:,4) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(5,4));
S15(:,5) = dot(CE(:,CeColInds(1,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(5,5));

% S21 = zeros(n, nx);
% S21(:,1) = dot(CE(:,colInds(2,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(1,1));
% S21(:,2) = dot(CE(:,colInds(2,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(1,2));
% S21(:,3) = dot(CE(:,colInds(2,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(1,3));
% S21(:,4) = dot(CE(:,colInds(2,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(1,4));
% S21(:,5) = dot(CE(:,colInds(2,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(1,5));

S22 = zeros(n, nx);
S22(:,1) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(2,1));
S22(:,2) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(2,2));
S22(:,3) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(2,3));
S22(:,4) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(2,4));
S22(:,5) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(2,5));

S23 = zeros(n, nx);
S23(:,1) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(3,1));
S23(:,2) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(3,2));
S23(:,3) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(3,3));
S23(:,4) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(3,4));
S23(:,5) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(3,5));

S24 = zeros(n, nx);
S24(:,1) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(4,1));
S24(:,2) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(4,2));
S24(:,3) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(4,3));
S24(:,4) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(4,4));
S24(:,5) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(4,5));

S25 = zeros(n, nx);
S25(:,1) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(5,1));
S25(:,2) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(5,2));
S25(:,3) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(5,3));
S25(:,4) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(5,4));
S25(:,5) = dot(CE(:,CeColInds(2,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(5,5));

% S31 = zeros(n, nx);
% S31(:,1) = dot(CE(:,colInds(3,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(1,1));
% S31(:,2) = dot(CE(:,colInds(3,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(1,2));
% S31(:,3) = dot(CE(:,colInds(3,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(1,3));
% S31(:,4) = dot(CE(:,colInds(3,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(1,4));
% S31(:,5) = dot(CE(:,colInds(3,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(1,5));
% 
% S32 = zeros(n, nx);
% S32(:,1) = dot(CE(:,colInds(3,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(2,1));
% S32(:,2) = dot(CE(:,colInds(3,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(2,2));
% S32(:,3) = dot(CE(:,colInds(3,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(2,3));
% S32(:,4) = dot(CE(:,colInds(3,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(2,4));
% S32(:,5) = dot(CE(:,colInds(3,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(2,5));

S33 = zeros(n, nx);
S33(:,1) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(3,1));
S33(:,2) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(3,2));
S33(:,3) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(3,3));
S33(:,4) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(3,4));
S33(:,5) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(3,5));

S34 = zeros(n, nx);
S34(:,1) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(4,1));
S34(:,2) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(4,2));
S34(:,3) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(4,3));
S34(:,4) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(4,4));
S34(:,5) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(4,5));

S35 = zeros(n, nx);
S35(:,1) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(5,1));
S35(:,2) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(5,2));
S35(:,3) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(5,3));
S35(:,4) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(5,4));
S35(:,5) = dot(CE(:,CeColInds(3,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(5,5));

% S41 = zeros(n, nx);
% S41(:,1) = dot(CE(:,colInds(4,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(1,1));
% S41(:,2) = dot(CE(:,colInds(4,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(1,2));
% S41(:,3) = dot(CE(:,colInds(4,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(1,3));
% S41(:,4) = dot(CE(:,colInds(4,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(1,4));
% S41(:,5) = dot(CE(:,colInds(4,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(1,5));
% 
% S42 = zeros(n, nx);
% S42(:,1) = dot(CE(:,colInds(4,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(2,1));
% S42(:,2) = dot(CE(:,colInds(4,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(2,2));
% S42(:,3) = dot(CE(:,colInds(4,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(2,3));
% S42(:,4) = dot(CE(:,colInds(4,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(2,4));
% S42(:,5) = dot(CE(:,colInds(4,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(2,5));
% 
% S43 = zeros(n, nx);
% S43(:,1) = dot(CE(:,colInds(4,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(3,1));
% S43(:,2) = dot(CE(:,colInds(4,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(3,2));
% S43(:,3) = dot(CE(:,colInds(4,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(3,3));
% S43(:,4) = dot(CE(:,colInds(4,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(3,4));
% S43(:,5) = dot(CE(:,colInds(4,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(3,5));

S44 = zeros(n, nx);
S44(:,1) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(4,1));
S44(:,2) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(4,2));
S44(:,3) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(4,3));
S44(:,4) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(4,4));
S44(:,5) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(4,5));

S45 = zeros(n, nx);
S45(:,1) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(5,1));
S45(:,2) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(5,2));
S45(:,3) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(5,3));
S45(:,4) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(5,4));
S45(:,5) = dot(CE(:,CeColInds(4,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(5,5));

% S51 = zeros(n, nx);
% S51(:,1) = dot(CE(:,colInds(5,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(1,1));
% S51(:,2) = dot(CE(:,colInds(5,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(1,2));
% S51(:,3) = dot(CE(:,colInds(5,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(1,3));
% S51(:,4) = dot(CE(:,colInds(5,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(1,4));
% S51(:,5) = dot(CE(:,colInds(5,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(1,5));
% 
% S52 = zeros(n, nx);
% S52(:,1) = dot(CE(:,colInds(5,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(2,1));
% S52(:,2) = dot(CE(:,colInds(5,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(2,2));
% S52(:,3) = dot(CE(:,colInds(5,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(2,3));
% S52(:,4) = dot(CE(:,colInds(5,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(2,4));
% S52(:,5) = dot(CE(:,colInds(5,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(2,5));
% 
% S53 = zeros(n, nx);
% S53(:,1) = dot(CE(:,colInds(5,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(3,1));
% S53(:,2) = dot(CE(:,colInds(5,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(3,2));
% S53(:,3) = dot(CE(:,colInds(5,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(3,3));
% S53(:,4) = dot(CE(:,colInds(5,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(3,4));
% S53(:,5) = dot(CE(:,colInds(5,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(3,5));
% 
% S54 = zeros(n, nx);
% S54(:,1) = dot(CE(:,colInds(5,:)),G(:,rowInds(1,:)),2).*CE(:,rowInds(4,1));
% S54(:,2) = dot(CE(:,colInds(5,:)),G(:,rowInds(2,:)),2).*CE(:,rowInds(4,2));
% S54(:,3) = dot(CE(:,colInds(5,:)),G(:,rowInds(3,:)),2).*CE(:,rowInds(4,3));
% S54(:,4) = dot(CE(:,colInds(5,:)),G(:,rowInds(4,:)),2).*CE(:,rowInds(4,4));
% S54(:,5) = dot(CE(:,colInds(5,:)),G(:,rowInds(5,:)),2).*CE(:,rowInds(4,5));

S55 = zeros(n, nx);
S55(:,1) = dot(CE(:,CeColInds(5,:)),G(:,rowInds(1,:)),2).*CE(:,CeRowInds(5,1));
S55(:,2) = dot(CE(:,CeColInds(5,:)),G(:,rowInds(2,:)),2).*CE(:,CeRowInds(5,2));
S55(:,3) = dot(CE(:,CeColInds(5,:)),G(:,rowInds(3,:)),2).*CE(:,CeRowInds(5,3));
S55(:,4) = dot(CE(:,CeColInds(5,:)),G(:,rowInds(4,:)),2).*CE(:,CeRowInds(5,4));
S55(:,5) = dot(CE(:,CeColInds(5,:)),G(:,rowInds(5,:)),2).*CE(:,CeRowInds(5,5));

L11 = CEd(:,1) + dot(CE(:,CeColInds(1,:)),A(:,rowInds(1,:)),2) + dot(AT(:,colInds(1,:)),CE(:,CeRowInds(1,:)),2) - sum(S11,2) + Q(:,1);
L12 = CEd(:,2) + dot(CE(:,CeColInds(1,:)),A(:,rowInds(2,:)),2) + dot(AT(:,colInds(1,:)),CE(:,CeRowInds(2,:)),2) - sum(S12,2) + Q(:,2);
L13 = CEd(:,3) + dot(CE(:,CeColInds(1,:)),A(:,rowInds(3,:)),2) + dot(AT(:,colInds(1,:)),CE(:,CeRowInds(3,:)),2) - sum(S13,2) + Q(:,3);
L14 = CEd(:,4) + dot(CE(:,CeColInds(1,:)),A(:,rowInds(4,:)),2) + dot(AT(:,colInds(1,:)),CE(:,CeRowInds(4,:)),2) - sum(S14,2) + Q(:,4);
L15 = CEd(:,5) + dot(CE(:,CeColInds(1,:)),A(:,rowInds(5,:)),2) + dot(AT(:,colInds(1,:)),CE(:,CeRowInds(5,:)),2) - sum(S15,2) + Q(:,5);

%L21 = CEd(:,6) + dot(CE(:,colInds(2,:)),A(:,rowInds(1,:)),2) + dot(AT(:,colInds(2,:)),CE(:,rowInds(1,:)),2) - sum(S21,2) + Q(:,6);
L22 = CEd(:,6) + dot(CE(:,CeColInds(2,:)),A(:,rowInds(2,:)),2) + dot(AT(:,colInds(2,:)),CE(:,CeRowInds(2,:)),2) - sum(S22,2) + Q(:,7);
L23 = CEd(:,7) + dot(CE(:,CeColInds(2,:)),A(:,rowInds(3,:)),2) + dot(AT(:,colInds(2,:)),CE(:,CeRowInds(3,:)),2) - sum(S23,2) + Q(:,8);
L24 = CEd(:,8) + dot(CE(:,CeColInds(2,:)),A(:,rowInds(4,:)),2) + dot(AT(:,colInds(2,:)),CE(:,CeRowInds(4,:)),2) - sum(S24,2) + Q(:,9);
L25 = CEd(:,9) + dot(CE(:,CeColInds(2,:)),A(:,rowInds(5,:)),2) + dot(AT(:,colInds(2,:)),CE(:,CeRowInds(5,:)),2) - sum(S25,2) + Q(:,10);

%L31 = CEd(:,11) + dot(CE(:,colInds(3,:)),A(:,rowInds(1,:)),2) + dot(AT(:,colInds(3,:)),CE(:,rowInds(1,:)),2) - sum(S31,2) + Q(:,11);
%L32 = CEd(:,12) + dot(CE(:,colInds(3,:)),A(:,rowInds(2,:)),2) + dot(AT(:,colInds(3,:)),CE(:,rowInds(2,:)),2) - sum(S32,2) + Q(:,12);
L33 = CEd(:,10) + dot(CE(:,CeColInds(3,:)),A(:,rowInds(3,:)),2) + dot(AT(:,colInds(3,:)),CE(:,CeRowInds(3,:)),2) - sum(S33,2) + Q(:,13);
L34 = CEd(:,11) + dot(CE(:,CeColInds(3,:)),A(:,rowInds(4,:)),2) + dot(AT(:,colInds(3,:)),CE(:,CeRowInds(4,:)),2) - sum(S34,2) + Q(:,14);
L35 = CEd(:,12) + dot(CE(:,CeColInds(3,:)),A(:,rowInds(5,:)),2) + dot(AT(:,colInds(3,:)),CE(:,CeRowInds(5,:)),2) - sum(S35,2) + Q(:,15);

%L41 = CEd(:,16) + dot(CE(:,colInds(4,:)),A(:,rowInds(1,:)),2) + dot(AT(:,colInds(4,:)),CE(:,rowInds(1,:)),2) - sum(S41,2) + Q(:,16);
%L42 = CEd(:,17) + dot(CE(:,colInds(4,:)),A(:,rowInds(2,:)),2) + dot(AT(:,colInds(4,:)),CE(:,rowInds(2,:)),2) - sum(S42,2) + Q(:,17);
%L43 = CEd(:,18) + dot(CE(:,colInds(4,:)),A(:,rowInds(3,:)),2) + dot(AT(:,colInds(4,:)),CE(:,rowInds(3,:)),2) - sum(S43,2) + Q(:,18);
L44 = CEd(:,13) + dot(CE(:,CeColInds(4,:)),A(:,rowInds(4,:)),2) + dot(AT(:,colInds(4,:)),CE(:,CeRowInds(4,:)),2) - sum(S44,2) + Q(:,19);
L45 = CEd(:,14) + dot(CE(:,CeColInds(4,:)),A(:,rowInds(5,:)),2) + dot(AT(:,colInds(4,:)),CE(:,CeRowInds(5,:)),2) - sum(S45,2) + Q(:,20);

%L51 = CEd(:,21) + dot(CE(:,colInds(5,:)),A(:,rowInds(1,:)),2) + dot(AT(:,colInds(5,:)),CE(:,rowInds(1,:)),2) - sum(S51,2) + Q(:,21);
%L52 = CEd(:,22) + dot(CE(:,colInds(5,:)),A(:,rowInds(2,:)),2) + dot(AT(:,colInds(5,:)),CE(:,rowInds(2,:)),2) - sum(S52,2) + Q(:,22);
%L53 = CEd(:,23) + dot(CE(:,colInds(5,:)),A(:,rowInds(3,:)),2) + dot(AT(:,colInds(5,:)),CE(:,rowInds(3,:)),2) - sum(S53,2) + Q(:,23);
%L54 = CEd(:,24) + dot(CE(:,colInds(5,:)),A(:,rowInds(4,:)),2) + dot(AT(:,colInds(5,:)),CE(:,rowInds(4,:)),2) - sum(S54,2) + Q(:,24);
L55 = CEd(:,15) + dot(CE(:,CeColInds(5,:)),A(:,rowInds(5,:)),2) + dot(AT(:,colInds(5,:)),CE(:,CeRowInds(5,:)),2) - sum(S55,2) + Q(:,25);

% loss = [L11; L12; L13; L14; L15; ...
%         L21; L22; L23; L24; L25; ...
%         L31; L32; L33; L34; L35; ...
%         L41; L42; L43; L44; L45; ...
%         L51; L52; L53; L54; L55];
    
loss = [L11; L12; L13; L14; L15; ...
        L22; L23; L24; L25; ...
        L33; L34; L35; ...
        L44; L45; ...
        L55];