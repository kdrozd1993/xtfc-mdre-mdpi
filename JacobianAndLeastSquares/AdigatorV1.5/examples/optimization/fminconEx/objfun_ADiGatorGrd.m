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

function f = objfun_ADiGatorGrd(x)
global ADiGator_objfun_ADiGatorGrd
if isempty(ADiGator_objfun_ADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_objfun_ADiGatorGrd.objfun_ADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Objective Function File
cada1f1dx = x.dx(1);
cada1f1 = x.f(1);
cada1f2dx = exp(cada1f1).*cada1f1dx;
cada1f2 = exp(cada1f1);
cada1f3dx = x.dx(1);
cada1f3 = x.f(1);
cada1f4dx = 2.*cada1f3.^(2-1).*cada1f3dx;
cada1f4 = cada1f3^2;
cada1f5dx = 4.*cada1f4dx;
cada1f5 = 4*cada1f4;
cada1f6dx = x.dx(2);
cada1f6 = x.f(2);
cada1f7dx = 2.*cada1f6.^(2-1).*cada1f6dx;
cada1f7 = cada1f6^2;
cada1f8dx = 2.*cada1f7dx;
cada1f8 = 2*cada1f7;
cada1td1 = zeros(2,1);
cada1td1(1) = cada1f5dx;
cada1td1(2) = cada1td1(2) + cada1f8dx;
cada1f9dx = cada1td1;
cada1f9 = cada1f5 + cada1f8;
cada1f10dx = x.dx(1);
cada1f10 = x.f(1);
cada1f11dx = 4.*cada1f10dx;
cada1f11 = 4*cada1f10;
cada1f12dx = x.dx(2);
cada1f12 = x.f(2);
cada1td1 = zeros(2,1);
cada1td1(1) = cada1f12.*cada1f11dx;
cada1td1(2) = cada1td1(2) + cada1f11.*cada1f12dx;
cada1f13dx = cada1td1;
cada1f13 = cada1f11*cada1f12;
cada1td1 = cada1f9dx;
cada1td1 = cada1td1 + cada1f13dx;
cada1f14dx = cada1td1;
cada1f14 = cada1f9 + cada1f13;
cada1f15dx = x.dx(2);
cada1f15 = x.f(2);
cada1f16dx = 2.*cada1f15dx;
cada1f16 = 2*cada1f15;
cada1td1 = cada1f14dx;
cada1td1(2) = cada1td1(2) + cada1f16dx;
cada1f17dx = cada1td1;
cada1f17 = cada1f14 + cada1f16;
cada1f18dx = cada1f17dx;
cada1f18 = cada1f17 + 1;
cada1td1 = zeros(2,1);
cada1td1(1) = cada1f18.*cada1f2dx;
cada1td1 = cada1td1 + cada1f2.*cada1f18dx;
f.dx = cada1td1;
f.f = cada1f2*cada1f18;
%User Line: f = exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);
f.dx_size = 2;
f.dx_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_objfun_ADiGatorGrd
ADiGator_objfun_ADiGatorGrd = load('objfun_ADiGatorGrd.mat');
return
end