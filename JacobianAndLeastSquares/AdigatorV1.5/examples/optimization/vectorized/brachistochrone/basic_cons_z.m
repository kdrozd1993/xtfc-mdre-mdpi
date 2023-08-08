% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function [C,Ceq] = basic_cons_z(z,probinfo)
global ADiGator_basic_cons_z
if isempty(ADiGator_basic_cons_z); ADiGator_LoadData(); end
Gator1Data = ADiGator_basic_cons_z.basic_cons_z.Gator1Data;
% ADiGator Start Derivative Computations
X.dz = z.dz(Gator1Data.Index1);
X.f = z.f(probinfo.xind);
%User Line: X = z(probinfo.xind);
U.dz = z.dz(Gator1Data.Index2);
U.f = z.f(probinfo.uind);
%User Line: U = z(probinfo.uind);
tf.dz = z.dz(325);
tf.f = z.f(probinfo.tfind);
%User Line: tf = z(probinfo.tfind);
%User Line: % Calculate h
numintervals = probinfo.numintervals;
%User Line: numintervals = probinfo.numintervals;
h.dz = tf.dz./numintervals;
h.f = tf.f/numintervals;
%User Line: h = tf/numintervals;
%User Line: % Calculate F
cadainput2_1.dz = X.dz; cadainput2_1.f = X.f;
%User Line: cadainput2_1 = X;
cadainput2_2.dz = U.dz; cadainput2_2.f = U.f;
%User Line: cadainput2_2 = U;
cadaoutput2_1 = ADiGator_dynamics1(cadainput2_1,cadainput2_2);
% Call to function: dynamics
F.dz = cadaoutput2_1.dz; F.f = cadaoutput2_1.f;
%User Line: F = cadaoutput2_1;
%User Line: % Get Reference Indices
k = probinfo.k;
%User Line: k    = probinfo.k;
kbp1 = probinfo.kbp1;
%User Line: kbp1 = probinfo.kbp1;
kp1 = probinfo.kp1;
%User Line: kp1  = probinfo.kp1;
Xk.dz = X.dz(Gator1Data.Index3);
Xk.f = X.f(k,:);
%User Line: Xk    = X(k,:);
Xkbp1.dz = X.dz(Gator1Data.Index4);
Xkbp1.f = X.f(kbp1,:);
%User Line: Xkbp1 = X(kbp1,:);
Xkp1.dz = X.dz(Gator1Data.Index5);
Xkp1.f = X.f(kp1,:);
%User Line: Xkp1  = X(kp1,:);
Fk.dz = F.dz(Gator1Data.Index6);
Fk.f = F.f(k,:);
%User Line: Fk    = F(k,:);
Fkbp1.dz = F.dz(Gator1Data.Index7);
Fkbp1.f = F.f(kbp1,:);
%User Line: Fkbp1 = F(kbp1,:);
Fkp1.dz = F.dz(Gator1Data.Index8);
Fkp1.f = F.f(kp1,:);
%User Line: Fkp1  = F(kp1,:);
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index9) = Xkp1.dz;
cada1td1(Gator1Data.Index10) = cada1td1(Gator1Data.Index10) + Xk.dz;
cada1f1dz = cada1td1;
cada1f1 = Xkp1.f + Xk.f;
cada1f2dz = 0.5.*cada1f1dz;
cada1f2 = 0.5*cada1f1;
cada1td1 = zeros(360,1);
cada1td1(Gator1Data.Index11) = Xkbp1.dz;
cada1td1(Gator1Data.Index12) = cada1td1(Gator1Data.Index12) + -cada1f2dz;
cada1f3dz = cada1td1;
cada1f3 = Xkbp1.f - cada1f2;
cada1f4dz = h.dz./8;
cada1f4 = h.f/8;
cada1td1 = zeros(400,1);
cada1td1(Gator1Data.Index13) = Fk.dz;
cada1td1(Gator1Data.Index14) = cada1td1(Gator1Data.Index14) + -Fkp1.dz;
cada1f5dz = cada1td1;
cada1f5 = Fk.f - Fkp1.f;
cada1tempdz = cada1f4dz(Gator1Data.Index15);
cada1td1 = zeros(520,1);
cada1td1(Gator1Data.Index16) = cada1f5(:).*cada1tempdz;
cada1td1(Gator1Data.Index17) = cada1td1(Gator1Data.Index17) + cada1f4.*cada1f5dz;
cada1f6dz = cada1td1;
cada1f6 = cada1f4*cada1f5;
cada1td1 = zeros(880,1);
cada1td1(Gator1Data.Index18) = cada1f3dz;
cada1td1(Gator1Data.Index19) = cada1td1(Gator1Data.Index19) + -cada1f6dz;
C1.dz = cada1td1;
C1.f = cada1f3 - cada1f6;
%User Line: C1 = Xkbp1 - 1/2.*(Xkp1+Xk) - h/8.*(Fk-Fkp1);
cada1td1 = zeros(240,1);
cada1td1(Gator1Data.Index20) = Xkp1.dz;
cada1td1(Gator1Data.Index21) = cada1td1(Gator1Data.Index21) + -Xk.dz;
cada1f1dz = cada1td1;
cada1f1 = Xkp1.f - Xk.f;
cada1f2dz = h.dz./6;
cada1f2 = h.f/6;
cada1f3dz = 4.*Fkbp1.dz;
cada1f3 = 4*Fkbp1.f;
cada1td1 = zeros(400,1);
cada1td1(Gator1Data.Index22) = Fkp1.dz;
cada1td1(Gator1Data.Index23) = cada1td1(Gator1Data.Index23) + cada1f3dz;
cada1f4dz = cada1td1;
cada1f4 = Fkp1.f + cada1f3;
cada1td1 = zeros(600,1);
cada1td1(Gator1Data.Index24) = cada1f4dz;
cada1td1(Gator1Data.Index25) = cada1td1(Gator1Data.Index25) + Fk.dz;
cada1f5dz = cada1td1;
cada1f5 = cada1f4 + Fk.f;
cada1tempdz = cada1f2dz(Gator1Data.Index26);
cada1td1 = zeros(720,1);
cada1td1(Gator1Data.Index27) = cada1f5(:).*cada1tempdz;
cada1td1(Gator1Data.Index28) = cada1td1(Gator1Data.Index28) + cada1f2.*cada1f5dz;
cada1f6dz = cada1td1;
cada1f6 = cada1f2*cada1f5;
cada1td1 = zeros(960,1);
cada1td1(Gator1Data.Index29) = cada1f1dz;
cada1td1(Gator1Data.Index30) = cada1td1(Gator1Data.Index30) + -cada1f6dz;
C2.dz = cada1td1;
C2.f = cada1f1 - cada1f6;
%User Line: C2 = Xkp1 - Xk - h/6.*(Fkp1 + 4.*Fkbp1 + Fk);
cada1f1dz = C1.dz;
cada1f1 = C1.f(:);
cada1f2dz = C2.dz;
cada1f2 = C2.f(:);
cada1td1 = zeros(1840,1);
cada1td1(Gator1Data.Index31) = cada1f1dz;
cada1td1(Gator1Data.Index32) = cada1f2dz;
Ceq.dz = cada1td1;
Ceq.f = [cada1f1;cada1f2];
%User Line: Ceq = [C1(:);C2(:)];
C.f = [];
%User Line: C =[];
Ceq.dz_size = [240,325];
Ceq.dz_location = Gator1Data.Index33;
end
function Xdot = ADiGator_dynamics1(X,U)
global ADiGator_basic_cons_z
Gator1Data = ADiGator_basic_cons_z.ADiGator_dynamics1.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Brachistochrone Dynamics
cada1f1 = size(X.f);
Xdot.f = zeros(cada1f1);
%User Line: Xdot = zeros(size(X));
X3.dz = X.dz(Gator1Data.Index1);
X3.f = X.f(:,3);
%User Line: X3 = X(:,3);
cada1f1dz = cos(U.f(:)).*U.dz;
cada1f1 = sin(U.f);
cada1td1 = zeros(162,1);
cada1td1(Gator1Data.Index2) = cada1f1(:).*X3.dz;
cada1td1(Gator1Data.Index3) = cada1td1(Gator1Data.Index3) + X3.f(:).*cada1f1dz;
cada1f2dz = cada1td1;
cada1f2 = X3.f.*cada1f1;
Xdot.dz = cada1f2dz;
Xdot.f(:,1) = cada1f2;
%User Line: Xdot(:,1) = X3.*sin(U);
cada1f1dz = -X3.dz;
cada1f1 = uminus(X3.f);
cada1f2dz = -sin(U.f(:)).*U.dz;
cada1f2 = cos(U.f);
cada1td1 = zeros(162,1);
cada1td1(Gator1Data.Index4) = cada1f2(:).*cada1f1dz;
cada1td1(Gator1Data.Index5) = cada1td1(Gator1Data.Index5) + cada1f1(:).*cada1f2dz;
cada1f3dz = cada1td1;
cada1f3 = cada1f1.*cada1f2;
cada1td1 = zeros(324,1);
cada1td1(Gator1Data.Index6) = cada1f3dz;
cada1td1(Gator1Data.Index7) = Xdot.dz(Gator1Data.Index8);
Xdot.dz = cada1td1;
Xdot.f(:,2) = cada1f3;
%User Line: Xdot(:,2) = -X3.*cos(U);
cada1f1dz = -sin(U.f(:)).*U.dz;
cada1f1 = cos(U.f);
cada1f2dz = 9.81.*cada1f1dz;
cada1f2 = 9.81*cada1f1;
cada1td1 = zeros(405,1);
cada1td1(Gator1Data.Index9) = cada1f2dz;
cada1td1(Gator1Data.Index10) = Xdot.dz(Gator1Data.Index11);
Xdot.dz = cada1td1;
Xdot.f(:,3) = cada1f2;
%User Line: Xdot(:,3) = 9.81.*cos(U);
end


function ADiGator_LoadData()
global ADiGator_basic_cons_z
ADiGator_basic_cons_z = load('basic_cons_z.mat');
return
end