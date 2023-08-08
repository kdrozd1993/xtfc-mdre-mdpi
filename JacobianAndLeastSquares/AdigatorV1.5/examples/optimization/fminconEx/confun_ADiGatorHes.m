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

function [c,ceq] = confun_ADiGatorHes(x)
global ADiGator_confun_ADiGatorHes
if isempty(ADiGator_confun_ADiGatorHes); ADiGator_LoadData(); end
Gator1Data = ADiGator_confun_ADiGatorHes.confun_ADiGatorHes.Gator1Data;
Gator2Data = ADiGator_confun_ADiGatorHes.confun_ADiGatorHes.Gator2Data;
% ADiGator Start Derivative Computations
%User Line: % Constraints function
cada1f1dx = x.dx(1);
cada1f1 = x.f(1);
cada1f2dx = x.dx(2);
cada1f2 = x.f(2);
cada1td1 =  zeros(2,1);
cada2f1dx = cada1f1dx.*cada1f2dx;
cada2f1 = cada1f2*cada1f1dx;
cada1td1dx = cada2f1dx;
cada1td1(1) = cada2f1;
cada2f1 = cada1td1(2);
cada2f2dx = cada1f2dx.*cada1f1dx;
cada2f2 = cada1f1*cada1f2dx;
cada2f3dx = cada2f2dx;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(2,1);
cada2td1(1) = cada2f3dx;
cada2td1(2) = cada1td1dx(1);
cada1td1dx = cada2td1;
cada1td1(2) = cada2f3;
cada1f3dxdx = cada1td1dx; cada1f3dx = cada1td1;
cada1f3 = cada1f1*cada1f2;
cada1f4dxdx = cada1f3dxdx; cada1f4dx = cada1f3dx;
cada1f4 = 1.5 + cada1f3;
cada1f5dx = x.dx(1);
cada1f5 = x.f(1);
cada1td1dx = cada1f4dxdx; cada1td1 = cada1f4dx;
cada2f1dx = cada1td1dx(2);
cada2f1 = cada1td1(1);
cada2f2 = uminus(cada1f5dx);
cada2f3dx = cada2f1dx;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(2,1);
cada2td1(2) = cada2f3dx;
cada2td1(1) = cada1td1dx(1);
cada1td1dx = cada2td1;
cada1td1(1) = cada2f3;
cada1f6dxdx = cada1td1dx; cada1f6dx = cada1td1;
cada1f6 = cada1f4 - cada1f5;
cada1f7dx = x.dx(2);
cada1f7 = x.f(2);
cada1td1dx = cada1f6dxdx; cada1td1 = cada1f6dx;
cada2f1dx = cada1td1dx(1);
cada2f1 = cada1td1(2);
cada2f2 = uminus(cada1f7dx);
cada2f3dx = cada2f1dx;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(2,1);
cada2td1(1) = cada2f3dx;
cada2td1(2) = cada1td1dx(2);
cada1td1dx = cada2td1;
cada1td1(2) = cada2f3;
cada1f8dxdx = cada1td1dx; cada1f8dx = cada1td1;
cada1f8 = cada1f6 - cada1f7;
cada1temp1 = Gator1Data.Data1;
c.dxdx = cada1f8dxdx; c.dx = cada1f8dx;
c.f = cada1temp1;
c.f(1) = cada1f8;
%User Line: c(1) = 1.5 + x(1) * x(2) - x(1) - x(2);
cada1f1dx = x.dx(1);
cada1f1 = x.f(1);
cada1f2dx = uminus(cada1f1dx);
cada1f2 = uminus(cada1f1);
cada1f3dx = x.dx(2);
cada1f3 = x.f(2);
cada1td1 =  zeros(2,1);
cada2f1dx = cada1f2dx.*cada1f3dx;
cada2f1 = cada1f3*cada1f2dx;
cada1td1dx = cada2f1dx;
cada1td1(1) = cada2f1;
cada2f1 = cada1td1(2);
cada2f2dx = cada1f3dx.*cada1f2dx;
cada2f2 = cada1f2*cada1f3dx;
cada2f3dx = cada2f2dx;
cada2f3 = cada2f1 + cada2f2;
cada2td1 = zeros(2,1);
cada2td1(1) = cada2f3dx;
cada2td1(2) = cada1td1dx(1);
cada1td1dx = cada2td1;
cada1td1(2) = cada2f3;
cada1f4dxdx = cada1td1dx; cada1f4dx = cada1td1;
cada1f4 = cada1f2*cada1f3;
cada1f5dxdx = cada1f4dxdx; cada1f5dx = cada1f4dx;
cada1f5 = cada1f4 - 10;
cada1td1 =  zeros(4,1);
cada1td1dx = cada1f5dxdx;
cada1td1(Gator1Data.Index1) = cada1f5dx;
cada2f1dx = c.dxdx(Gator2Data.Index1);
cada2f1 = c.dx(Gator1Data.Index3);
cada2td1 = zeros(4,1);
cada2td1(Gator2Data.Index2) = cada2f1dx;
cada2td1(Gator2Data.Index3) = cada1td1dx(Gator2Data.Index4);
cada1td1dx = cada2td1;
cada1td1(Gator1Data.Index2) = cada2f1;
c.dxdx = cada1td1dx; c.dx = cada1td1;
c.f(2) = cada1f5;
%User Line: c(2) = -x(1) * x(2)-10;
%User Line: % No nonlinear equality constraints
ceq.f = [];
%User Line: ceq=[];
c.dx_size = [2 2];
c.dx_location = Gator1Data.Index4;
c.dxdx_size = [c.dx_size,2];
c.dxdx_location = [c.dx_location(Gator2Data.Index5,:), Gator2Data.Index6];
end


function ADiGator_LoadData()
global ADiGator_confun_ADiGatorHes
ADiGator_confun_ADiGatorHes = load('confun_ADiGatorHes.mat');
return
end