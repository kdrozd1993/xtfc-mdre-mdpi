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

function f = objfun_ADiGatorHes(x)
global ADiGator_objfun_ADiGatorHes
if isempty(ADiGator_objfun_ADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_objfun_ADiGatorHes.objfun_ADiGatorHes.Gator1Data;
Gator2Data = ADiGator_objfun_ADiGatorHes.objfun_ADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: % Objective Function File
cada1f1dx = x.dx(1);
cada1f1 = x.f(1);
cada2f1dx = exp(cada1f1).*cada1f1dx;
cada2f1 = exp(cada1f1);
cada1f2dxdx = cada1f1dx.*cada2f1dx;
cada1f2dx = cada2f1*cada1f1dx;
cada1f2 = exp(cada1f1);
cada1f3dx = x.dx(1);
cada1f3 = x.f(1);
cada2f1dx = 1.*cada1f3.^(1-1).*cada1f3dx;
cada2f1 = cada1f3^1;
cada2f2dx = 2.*cada2f1dx;
cada2f2 = 2*cada2f1;
cada1f4dxdx = cada1f3dx.*cada2f2dx;
cada1f4dx = cada2f2*cada1f3dx;
cada1f4 = cada1f3^2;
cada1f5dxdx = 4.*cada1f4dxdx;
cada1f5dx = 4*cada1f4dx;
cada1f5 = 4*cada1f4;
cada1f6dx = x.dx(2);
cada1f6 = x.f(2);
cada2f1dx = 1.*cada1f6.^(1-1).*cada1f6dx;
cada2f1 = cada1f6^1;
cada2f2dx = 2.*cada2f1dx;
cada2f2 = 2*cada2f1;
cada1f7dxdx = cada1f6dx.*cada2f2dx;
cada1f7dx = cada2f2*cada1f6dx;
cada1f7 = cada1f6^2;
cada1f8dxdx = 2.*cada1f7dxdx;
cada1f8dx = 2*cada1f7dx;
cada1f8 = 2*cada1f7;
cada1td1 =  zeros(2,1);
cada1td1dx = cada1f5dxdx;
cada1td1(1) = cada1f5dx;
cada2f1 = cada1td1(2);
cada2f2dx = cada1f8dxdx;
cada2f2 = cada2f1 + cada1f8dx;
cada2td1 = zeros(2,1);
cada2td1(2) = cada2f2dx;
cada2td1(1) = cada1td1dx(1);
cada1td1dx = cada2td1;
cada1td1(2) = cada2f2;
cada1f9dxdx = cada1td1dx; cada1f9dx = cada1td1;
cada1f9 = cada1f5 + cada1f8;
cada1f10dx = x.dx(1);
cada1f10 = x.f(1);
cada1f11dx = 4*cada1f10dx;
cada1f11 = 4*cada1f10;
cada1f12dx = x.dx(2);
cada1f12 = x.f(2);
cada1td1 =  zeros(2,1);
cada2f1dx = cada1f11dx.*cada1f12dx;
cada2f1 = cada1f12*cada1f11dx;
cada1td1dx = cada2f1dx;
cada1td1(1) = cada2f1;
cada2f1 = cada1td1(2);
cada2f2dx = cada1f12dx.*cada1f11dx;
cada2f2 = cada1f11*cada1f12dx;
cada2f3dx = cada2f2dx;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(2,1);
cada2td1(1) = cada2f3dx;
cada2td1(2) = cada1td1dx(1);
cada1td1dx = cada2td1;
cada1td1(2) = cada2f3;
cada1f13dxdx = cada1td1dx; cada1f13dx = cada1td1;
cada1f13 = cada1f11*cada1f12;
cada1td1dx = cada1f9dxdx; cada1td1 = cada1f9dx;
cada2td1 = zeros(4,1);
cada2td1(Gator2Data.Index1) = cada1td1dx;
cada2td1(Gator2Data.Index2) = cada2td1(Gator2Data.Index2) + cada1f13dxdx;
cada1td1dx = cada2td1;
cada1td1 = cada1td1 + cada1f13dx;
cada1f14dxdx = cada1td1dx; cada1f14dx = cada1td1;
cada1f14 = cada1f9 + cada1f13;
cada1f15dx = x.dx(2);
cada1f15 = x.f(2);
cada1f16dx = 2*cada1f15dx;
cada1f16 = 2*cada1f15;
cada1td1dx = cada1f14dxdx; cada1td1 = cada1f14dx;
cada2f1dx = cada1td1dx(Gator2Data.Index3);
cada2f1 = cada1td1(2);
cada2f2dx = cada2f1dx;
cada2f2 = cada2f1 + cada1f16dx;
cada2td1 = zeros(4,1);
cada2td1(Gator2Data.Index4) = cada2f2dx;
cada2td1(Gator2Data.Index5) = cada1td1dx(Gator2Data.Index6);
cada1td1dx = cada2td1;
cada1td1(2) = cada2f2;
cada1f17dxdx = cada1td1dx; cada1f17dx = cada1td1;
cada1f17 = cada1f14 + cada1f16;
cada1f18dxdx = cada1f17dxdx; cada1f18dx = cada1f17dx;
cada1f18 = cada1f17 + 1;
cada1td1 =  zeros(2,1);
cada2td1 = cada1f2dx.*cada1f18dx;
cada2td1(1) = cada2td1(1) + cada1f18.*cada1f2dxdx;
cada2f1dx = cada2td1;
cada2f1 = cada1f18*cada1f2dx;
cada1td1dx = cada2f1dx;
cada1td1(1) = cada2f1;
cada2tempdx = cada1f2dx(Gator2Data.Index7);
cada2td1 = zeros(4,1);
cada2td1(Gator2Data.Index8) = cada1f18dx(:).*cada2tempdx;
cada2td1 = cada2td1 + cada1f2.*cada1f18dxdx;
cada2f1dx = cada2td1;
cada2f1 = cada1f2*cada1f18dx;
cada2td1 = zeros(4,1);
cada2td1(Gator2Data.Index9) = cada1td1dx;
cada2td1 = cada2td1 + cada2f1dx;
cada1td1dx = cada2td1;
cada1td1 = cada1td1 + cada2f1;
f.dxdx = cada1td1dx; f.dx = cada1td1;
f.f = cada1f2*cada1f18;
%User Line: f = exp(x(1))*(4*x(1)^2+2*x(2)^2+4*x(1)*x(2)+2*x(2)+1);
f.dx_size = 2;
f.dx_location = Gator1Data.Index1;
f.dxdx_size = [f.dx_size,2];
f.dxdx_location = [f.dx_location(Gator2Data.Index10,:), Gator2Data.Index11];
end


function ADiGator_LoadData()
global ADiGator_objfun_ADiGatorHes
ADiGator_objfun_ADiGatorHes = load('objfun_ADiGatorHes.mat');
return
end