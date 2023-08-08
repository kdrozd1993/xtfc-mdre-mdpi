%{ 

Author: Kristofer Drozd
Created: 12/06/2021

Description:

A script for testing the XTFC for solving a transformed differential 
ricatti equation. The problem tested comes from the paper Solving the 
Matrix Differential Ricatti Equation: A Lyapunov Equation Approach.

%}

baseFolder = 'C:\Users\kdrozd\Documents\MATLAB\xtfc-mdre-mdpi';
addpath([baseFolder, '\ActivationsAndCollocations\']);
addpath([baseFolder, '\JacobianAndLeastSquares\']);
addpath(genpath([baseFolder, '\JacobianAndLeastSquares\']));
addpath([baseFolder, '\Benchmarks\']);
addpath([baseFolder, '\OdeSolvers\']);
addpath(baseFolder);

clear all
clc
close all





%% Setting up the Problem

% How should the jacobian be calculatted?
jacCase = 'adigator';
%jacCase = 'numeric';

% Least squares type:
lsType = 7;
tol = 1e-13;

% Loss function
lossFunc = str2func('TFC_MDLE_Loss_5');

% Setting up the tim intervals
t00 = 0; % Initial time
T = .001; % IVP time horizon
tff = .5; % Final time
nT = tff/T; % Number of propagations
timeStore = t00:T:tff;
Sff = diag([0.05, 0.05, 0.01, 0.01, 0.01]); % Ricatti cost matrix at final time

n = 20; % # of collocation points for each IVP
m = 15; % # of basis functions for each IVP
nu = 2; % # of control variables
nx = 5; % # of state variables

% Anonymous functions for the LTV matrices A, B, Q, and R.
ltvMatFuncs.Q   = @(t) eye(nx);
ltvMatFuncs.R   = @(t) eye(nu);
ltvMatFuncs.AB  = @(t) ABMats(t);

% Setting up structure for computing riccati basis matrices.
inputStruct.n = n;
inputStruct.m = m;
inputStruct.freeFunc = 'gaussian_elm';
inputStruct.collocType = 'chebyshev_gauss_lobatto';
inputStruct.tau0 = -1;
inputStruct.tauf = 1;

p = haltonset(1, 'Skip', 1, 'Leap', 0);
%p = scramble(p,'MatousekAffineOwen');
inputStruct.w = (net(p,m)-0.5)*2;
inputStruct.b = net(p,m)*1; %ones(m,1);

% inputStruct.w = unifrnd(-1, 1, m, 1);
% inputStruct.b = unifrnd(0, 1, m, 1);

[H, Hd, tau] = RiccatiBasisCalculator(inputStruct);

c = (tau(end)-tau(1))/(T);
t = t00 + (1/c)*(tau - tau(1)); % Collocation points in the time domain

% Obtaining necessary matrices.
tic
[mdre3DMats, mdreVecMats] = MdreConcatenatedMatricesFromFuncs(t, nx, nu, ltvMatFuncs);
[mdle3DMats, mdleVecMats] = MdleConcatenatedMatrices(Sff, mdre3DMats);
S0ff = mdle3DMats.S0(:,:,end);
toc

% Setting up a structure of inputs for computing the loss.
tfcVariables.n = n;
tfcVariables.tau = tau;
tfcVariables.m = m;
tfcVariables.H = H;
tfcVariables.Hd = Hd;

S0fft = S0ff.';
ss  = tril(true(size(S0fft)), 0);
S0fvut  = S0fft(ss).';
S0fvutr = repmat(S0fvut, n, 1);

setup.mdreVecMats = mdreVecMats;
setup.mdleVecMats = mdleVecMats;
setup.S0f = S0ff;
setup.S0fvutr = S0fvutr;
setup.t0 = t00;
setup.tf = tff;
setup.tfcVariables = tfcVariables;
setup.tfcVariables.nx = nx;
setup.lsType = lsType;
setup.tau = tau;





%% Computing the MDLE solution by splitting it up into multiple IVPs.

Sf = Sff;
S0f = S0ff;
S0fv = reshape(S0ff', 1, nx*nx);

xiStore = zeros(length(timeStore) - 1, m*nx*nx);
lossNormStore = zeros(n+(n-1)*(length(timeStore)-2), 1);
lossNormEndStore = zeros(length(timeStore) - 1, 1);
invCompStore = zeros(length(timeStore) - 1, 1);
S0fStore = zeros(length(timeStore) - 1, nx*nx);
ceStore = zeros(n+(n-1)*(length(timeStore)-2), nx*nx);
ceStore(end,:) = reshape(S0f', 1, nx*nx);
ricStore = zeros(n+(n-1)*(length(timeStore)-2), nx*nx);
ricStore(end,:) = reshape(Sf', 1, nx*nx);
ceEndStore = zeros(length(timeStore), nx*nx);
ceEndStore(end,:) = reshape(S0f', 1, nx*nx);
ricEndStore = zeros(length(timeStore), nx*nx);
ricEndStore(end,:) = reshape(Sf', 1, nx*nx);
timeInterpStore = zeros((length(timeStore)-2)*(n-1)+n, 1);
timeInterpStore(end) = tff;

tic
for i = 1:(length(timeStore)-1)

    t0 = timeStore(end-i);
    tf = timeStore(end-i+1);

    c = (tau(end)-tau(1))/(tf - t0);
    t = t0 + (1/c)*(tau - tau(1));
    
    setup.mdreVecMats = mdreVecMats;
    setup.mdleVecMats = mdleVecMats;
    setup.S0f = S0f;
    setup.t0 = t0;
    setup.tf = tf;
    
    S0ft = S0f';
    ss  = tril(true(size(S0ft)), 0);
    S0fvut  = S0ft(ss).';
    S0fvutr = repmat(S0fvut, n, 1);
    setup.S0fvutr = S0fvutr;

    S0fStore(end-(i-1),:) = S0fv;
    
    %tic
    [xi, CE, lossNorm, compTime] = TFC_MDLE_LS(lossFunc, setup);
    %toc

    lossNormEnd = sqrt(sum(lossNorm.^2));

    CE3 = reshape(CE', nx, nx, n);
    S0f = CE3(:,:,1);
    S0fv = reshape(S0f', 1, nx*nx);
    Sf = mdle3DMats.Sm(:,:,1) + inv(S0f);
    
    ric = zeros(nx, nx, n);
    ricv = zeros(size(CE, 1), nx*nx);
    for j = 1:size(CE, 1)
        ric(:,:,j) = (mdle3DMats.Sm(:,:,j) + inv(CE3(:,:,j)));
        ricv(j,:) = reshape(ric(:,:,j)', 1, nx*nx);
    end
    
    ricStore((end-i*(n-1)):(end-(i-1)*(n-1)-1),:) =  ricv(1:end-1,:);
    ricEndStore(end-i,:) =  ricv(1,:);
    ceStore((end-i*(n-1)):(end-(i-1)*(n-1)-1),:) = CE(1:end-1,:);
    ceEndStore(end-i,:) = CE(1,:);
    xiStore(end-(i-1),:) = reshape(xi, 1, nx*nx*m);
    invCompStore(end-(i-1),1) = sum(compTime.inversionComputationTime);
    timeInterpStore(n+(n-1)*(nT-i-1):n+(n-1)*(nT-i)-1) = t(1:end-1,1);
    lossNormEndStore(end-i+1,:) = lossNormEnd;
    
    if i == 1
        lossNormStore(end-n+1:end,1) = ...
            lossNorm;
    else
        lossNormStore((end-i*(n-1)):(end-(i-1)*(n-1)-1),1) = ...
            lossNorm(1:end-1,:);
    end
    
end
toc





%% Intperoloation Test

% Interpolution time points
timeInterp = unifrnd(t00, tff, 50, 1);

% Constrained expression results at interpolation points
[ceInterp] = TFC_MDRE_Interp(timeInterp, timeStore, xiStore, ...
    S0fStore, inputStruct);





%% Benchmarks

timeSpanKE = timeInterpStore;
[mdre3DMats, ~] = MdreConcatenatedMatricesFromFuncs(timeSpanKE, nx, nu, ltvMatFuncs);
SKE = KalmanEnglar(mdre3DMats, timeSpanKE, Sff);

SKEV = zeros(length(timeSpanKE), nx*nx);
for j = 1:length(timeSpanKE)
    SKEV(j,:) = reshape(SKE(:,:,j)', 1, nx*nx);
end

resNorm_PinnKE = vecnorm((SKEV-ricStore), 2, 2);

[SL] = LyapunovApproachRiccati(mdre3DMats, timeSpanKE, Sff);

SLV = zeros(length(timeSpanKE), nx*nx);
for j = 1:length(timeSpanKE)
    SLV(j,:) = reshape(SL(:,:,j)', 1, nx*nx);
end

resNorm_LyapKE = vecnorm((SKEV-SLV), 2, 2);

resNorm_PinnLyap = vecnorm((SLV-ricStore), 2, 2);

SRK4 = TrajectoryLqr(timeSpanKE, ltvMatFuncs, Sff, tol);

resNorm_PinnRK4 = vecnorm((SRK4-ricStore), 2, 2);

resNorm_RK4KE = vecnorm((SRK4-SKEV), 2, 2);

f = @() KalmanEnglar(mdre3DMats, timeSpanKE, Sff);
timeit(f)

f = @() LyapunovApproachRiccati(mdre3DMats, timeSpanKE, Sff);
timeit(f)

f = @() TrajectoryLqr(timeSpanKE, ltvMatFuncs, Sff, tol);
timeit(f)

f = @() XtfcMdleInterpTimeTest(timeStore, xiStore, S0fStore, inputStruct, mdle3DMats);
timeit(f)

f = @() interpn(timeInterpStore, SKEV(:,:), 0.4711);
timeit(f)

f = @() interpn(timeInterpStore, SLV(:,:), 0.4711);
timeit(f)

f = @() interpn(timeInterpStore, SRK4(:,:), 0.4711);
timeit(f)





%% Saving Results

results.ce = ceStore;
results.ceEnd = ceEndStore;
results.xi = xiStore;
results.lossNorm.interp = lossNormStore;
results.lossNorm.end = lossNormEndStore;
results.time.interp = timeInterpStore;
results.time.end = (t00:T:tff)';
results.compTime.inv = invCompStore;
results.residualNorm.PinnLKE = resNorm_PinnKE;
results.residualNorm.LyapKE = resNorm_LyapKE;
results.residualNorm.PinnLLyap = resNorm_PinnLyap;
results.residualNorm.PinnLRK4 = resNorm_PinnRK4;
results.residualNorm.RK4KE = resNorm_RK4KE;

options.amountCollocPts = m;
options.amountBasisFuncs = n;
options.timeHorizon = T;
options.freeFuncType = inputStruct.freeFunc;
options.collocType = inputStruct.collocType;
options.intialTau = inputStruct.tau0;
options.finalTau = inputStruct.tauf;

save('mdle_xtfc_mc.mat', 'results', 'options');





%% Plots

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%set(groot,'defaultLineLineWidth',2)
set(0,'DefaultAxesFontSize', 10)
%set(0,'DefaultLineMarkerSize',10);

f1 = figure(1);
plot(timeInterpStore, ceStore(:,1), 'r-'); hold on
plot(timeInterpStore, ceStore(:,7), 'b-'); hold on
plot(timeInterpStore, ceStore(:,13), 'c-'); hold on
plot(timeInterpStore, ceStore(:,19), 'm-'); hold on
plot(timeInterpStore, ceStore(:,25), 'g-'); hold on
ylabel('Lyapunov Matrix Values')
xlabel('Times (sec)')
legend('$E_{11}$', '$E_{22}$', '$E_{33}$', '$E_{44}$', '$E_{55}$', 'Location', 'NorthWest');
grid
% set(f1, 'Position',  [100, 100, 800, 500])
set(f1, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f2 = figure(2);
plot(timeInterpStore, lossNormStore, 'rx');
ylabel('Norm of Loss')
xlabel('Times (sec)')
grid
% set(f2, 'Position',  [100, 100, 800, 500])
set(f2, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f3 = figure(3);
plot(timeStore(2:end), lossNormEndStore, 'rx');
ylabel('Norm of Loss')
xlabel('Times (sec)')
grid
% set(f3, 'Position',  [100, 100, 800, 500])
set(f3, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f4 = figure(4);
plot(timeInterpStore, ceStore(:,1), 'r-'); hold on
plot(timeInterpStore, ceStore(:,7), 'b-'); hold on
plot(timeInterpStore, ceStore(:,13), 'c-'); hold on
plot(timeInterpStore, ceStore(:,19), 'm-'); hold on
plot(timeInterpStore, ceStore(:,25), 'g-'); hold on
plot(timeInterp, ceInterp(:,1), 'kx'); hold on
plot(timeInterp, ceInterp(:,7), 'kx'); hold on
plot(timeInterp, ceInterp(:,13), 'kx'); hold on
plot(timeInterp, ceInterp(:,19), 'kx'); hold on
plot(timeInterp, ceInterp(:,25), 'kx');
ylabel('Lyapunov Matrix Values')
xlabel('Times (s)')
legend('$P_{0_{11}}$', '$P_{0_{22}}$', '$P_{0_{33}}$', '$P_{0_{44}}$', '$P_{0_{55}}$', 'Location', 'NorthWest');
grid
% set(f4, 'Position',  [100, 100, 800, 500])
set(f4, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f5 = figure(5);
plot(timeSpanKE, resNorm_PinnKE, 'b-');
ylabel('Normwise PINN-KE Error')
xlabel('Times (sec)')
grid
% set(f5, 'Position',  [100, 100, 800, 500])
set(f5, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f6 = figure(6);
plot(timeSpanKE, resNorm_LyapKE, 'b-');
ylabel('Normwise KE-LA Error')
xlabel('Times (sec)')
grid
% set(f6, 'Position',  [100, 100, 800, 500])
set(f6, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f7 = figure(7);
plot(timeSpanKE, resNorm_PinnLyap, 'b-');
ylabel('Normwise PINN-LA Error')
xlabel('Times (sec)')
grid
% set(f7, 'Position',  [100, 100, 800, 500])
set(f7, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f8 = figure(8);
semilogy(timeSpanKE, resNorm_PinnKE, 'b-'); hold on
semilogy(timeSpanKE, resNorm_LyapKE, 'r-.'); hold on
semilogy(timeSpanKE, resNorm_RK4KE, 'k--'); hold on
legend('PINN', 'LA', 'RK4', 'Location', 'Best');
ylabel('Normwise Error wrt. KE')
xlabel('Time (sec)')
grid
% set(f8, 'Position',  [100, 100, 800, 500])
set(f8, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');

f9 = figure(9);
plot(timeInterpStore, ricStore(:,1), 'r-'); hold on
plot(timeInterpStore, ricStore(:,7), 'b-'); hold on
plot(timeInterpStore, ricStore(:,13), 'c-'); hold on
plot(timeInterpStore, ricStore(:,19), 'm-'); hold on
plot(timeInterpStore, ricStore(:,25), 'g-'); hold on
ylabel('Riccati Matrix Values')
xlabel('Times (sec)')
legend('$K_{11}$', '$K_{22}$', '$K_{33}$', '$K_{44}$', '$K_{55}$', 'Location', 'NorthWest');
grid
% set(f9, 'Position',  [100, 100, 800, 500])
set(f9, 'Units', 'inches', 'Position',  [1, 1, 5, 3])
set(gcf, 'PaperPositionMode', 'auto');





%% Nested functions

% A B matrix
function [AB] = ABMats(t)

A = [-16, -0.39, 27.20, 0, 0; ...
     0.01, -16.99, 0, 0, 12.47; ...
     15.11, 0, -53.60, -16.57, 71.78; ...
     -53.36, 0, 0, -107.20, 232.11; ...
     2.27, 69.10, 0, 2.273, -102.99];
 
 B = [11.12, -12.60; ...
     -3.61, 3.36; ...
     -21.91, 0; ...
     -53.60, 0; ...
     69.10, 0];
 
%  A = eye(5);
%  
%  B = zeros(5,2);
 
 AB = {A; B};

end

% If A or B is nonlinear then use this to linearize
function [AB] = LinearizedABMats(xd, ud, td, t, auxdata)

x = zeros(nx, 1);
for i = 1:nx
    x(i) = interp1(td, xd(:,i), t, 'spline');
end

u = zeros(nu, 1);
for i = 1:nu
    u(i) = interp1(td, ud(:,i), t, 'spline');
end

%A =  <Put linearized A func here> (x, u, auxdata);
%B = <Put linearized B func here> (x, u, auxdata);

AB = {A; B};

end