%{ 

Author: Kristofer Drozd
Created: 03/15/2022

Description:

A script for testing the XTFC for solving the differential ricatti 
equation. A Monte Carlo test is performed where the weights and bias of the
ELM are varied

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

% Number of Monte Carlo simulations.
mcn = 100;





%% Setting up the Problem

%How should the jacobian be calculatted?
jacCase = 'adigator';
%jacCase = 'numeric';

% Least squares type and tolerance:
lsType = 7;
tol = 1e-13;

% Setting up the tim intervals
t00 = 0; % Initial time
T = 0.01; % IVP time horizon
tff = 0.5; % Final time
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

for k1 = 1:mcn

    inputStruct.w = unifrnd(-1, 1, m, 1);
    inputStruct.b = unifrnd(0, 1, m, 1);

    [H, Hd, tau] = RiccatiBasisCalculator(inputStruct);

    c = (tau(end)-tau(1))/(T);
    t = t00 + (1/c)*(tau - tau(1)); % Collocation points in the time domain

    % Obtaining matrices in a concatenated format necessary for adigator to
    % work.
    [~, mdreVecMats] = MdreConcatenatedMatricesFromFuncs(t, nx, nu, ...
        ltvMatFuncs);

    % Setting up a structure of inputs for computing the loss.
    tfcVariables.n = n;
    tfcVariables.tau = tau;
    tfcVariables.m = m;
    tfcVariables.H = H;
    tfcVariables.Hd = Hd;

    Sfft = Sff.';
    ss  = tril(true(size(Sfft)), 0);
    Sfvut  = Sfft(ss).';
    Sfvutr = repmat(Sfvut, n, 1);

    setup.mdreVecMats = mdreVecMats;
    setup.Sf = Sff;
    setup.Sfvutr = Sfvutr;
    setup.t0 = t00;
    setup.tf = tff;
    setup.tfcVariables = tfcVariables;
    setup.tfcVariables.nx = nx;
    setup.lsType = lsType;
    setup.tau = tau;





    %% Adigator jacobian generation

    % This is used to generate the jacobian function with adigator. Itc an be 
    % commented out if the code already exists or jacoiban is calculated 
    % numerically.

    options = adigatorOptions('overwrite', 1, 'comments', 1, 'auxdata', 0, ...
        'unroll', 0);
    xi = adigatorCreateDerivInput([m, 1/2*nx*(nx+1)], 'xi');
    adigatorGenJacFile('TFC_MDRE_Loss_5',{xi, setup}, options);

    % Don't comment this out (Gives the name of the adigator generated jacobian
    % function.
    jacFunc = str2func(['TFC_MDRE_Loss_5', '_Jac']);





    %% Computing the DRE solution by splitting it up into multiple IVPs.

    Sf = Sff;
    Sfv = reshape(Sff', 1, nx*nx);

    xiStore = zeros(length(timeStore) - 1, m*nx*nx);
    lossNormStore = zeros(n+(n-1)*(length(timeStore)-2), 1);
    lossNormEndStore = zeros(length(timeStore) - 1, 1);
    jacCompStore = zeros(length(timeStore) - 1, 1);
    invCompStore = zeros(length(timeStore) - 1, 1);
    iterStore = zeros(length(timeStore) - 1, 1);
    SfStore = zeros(length(timeStore) - 1, nx*nx);
    ceStore = zeros(n+(n-1)*(length(timeStore)-2), nx*nx);
    ceStore(end,:) = reshape(Sf', 1, nx*nx);
    ricStore = zeros(n+(n-1)*(length(timeStore)-2), nx*nx);
    ricStore(end,:) = reshape(Sf', 1, nx*nx);
    ceEndStore = zeros(length(timeStore), nx*nx);
    ceEndStore(end,:) = reshape(Sf', 1, nx*nx);
    ricEndStore = zeros(length(timeStore), nx*nx);
    ricEndStore(end,:) = reshape(Sf', 1, nx*nx);
    timeInterpStore = zeros((length(timeStore)-2)*(n-1)+n, 1);
    timeInterpStore(end) = tff;

    for i = 1:(length(timeStore)-1)

        t0 = timeStore(end-i);
        tf = timeStore(end-i+1);

        c = (tau(end)-tau(1))/(tf - t0);
        t = t0 + (1/c)*(tau - tau(1));

        [~, mdreVecMats] = MdreConcatenatedMatricesFromFuncs(t, ...
            nx, nu, ltvMatFuncs);

        setup.mdreVecMats = mdreVecMats;
        setup.Sf = Sf;
        setup.t0 = t0;
        setup.tf = tf;

        Sft = Sf';
        ss  = tril(true(size(Sft)), 0);
        Sfvut  = Sft(ss).';
        Sfvutr = repmat(Sfvut, n, 1);
        setup.Sfvutr = Sfvutr;

        SfStore(end-(i-1),:) = Sfv;
        CEig = Sfvutr;

        [xi] = TFC_MDRE_InitialGuess(CEig, setup);

        tic
        [xi, iter, lossNorm, CE, compTime] = TFC_MDRE_ILS(xi, jacFunc, ...
            jacCase, 50, tol, setup);
        toc

        lossNormEnd = sqrt(sum(lossNorm.^2));

        CESF = reshape(CE', nx, nx, n);
        Sf = CESF(:,:,1);
        Sfv = reshape(Sf', 1, nx*nx);

        ceStore((end-i*(n-1)):(end-(i-1)*(n-1)-1),:) = CE(1:end-1,:);
        ceEndStore(end-i,:) = CE(1,:);
        xiStore(end-(i-1),:) = reshape(xi, 1, nx*nx*m);
        lossNormEndStore(end-(i-1),1) = lossNormEnd;
        iterStore(end-(i-1),1) = iter;
        jacCompStore(end-(i-1),1) = sum(compTime.jacobianComputationTime);
        invCompStore(end-(i-1),1) = sum(compTime.inversionComputationTime);
        timeInterpStore(n+(n-1)*(nT-i-1):n+(n-1)*(nT-i)-1) = t(1:end-1,1);

        if i == 1
            lossNormStore(end-n+1:end,1) = ...
                lossNorm;
        else
            lossNormStore((end-i*(n-1)):(end-(i-1)*(n-1)-1),1) = ...
                lossNorm(1:end-1,:);
        end

    end
    
    timeSpanKE = timeInterpStore;
    [mdre3DMats, ~] = MdreConcatenatedMatricesFromFuncs(timeSpanKE, nx, nu, ...
        ltvMatFuncs);
    SKE = KalmanEnglar(mdre3DMats, timeSpanKE, Sff);

    SKEV = zeros(length(timeSpanKE), nx*nx);
    for j = 1:length(timeSpanKE)
        SKEV(j,:) = reshape(SKE(:,:,j)', 1, nx*nx);
    end

    resNorm_PinnKE = vecnorm((SKEV-ceStore), 2, 2);

    [SL] = LyapunovApproachRiccati(mdre3DMats, timeSpanKE, Sff);

    SLV = zeros(length(timeSpanKE), nx*nx);
    for j = 1:length(timeSpanKE)
        SLV(j,:) = reshape(SL(:,:,j)', 1, nx*nx);
    end

    resNorm_LyapKE = vecnorm((SKEV-SLV), 2, 2);

    resNorm_PinnLyap = vecnorm((SLV-ceStore), 2, 2);

    SRK4 = TrajectoryLqr(timeSpanKE, ltvMatFuncs, Sff, tol);

    resNorm_PinnRK4 = vecnorm((SRK4-ceStore), 2, 2);

    resNorm_RK4KE = vecnorm((SRK4-SKEV), 2, 2);
    
    results(k1).ce = ceStore;
    results(k1).ceEnd = ceEndStore;
    results(k1).xi = xiStore;
    results(k1).lossNorm.interp = lossNormStore;
    results(k1).lossNorm.end = lossNormEndStore;
    results(k1).time.interp = timeInterpStore;
    results(k1).time.end = (t00:T:tff)';
    results(k1).compTime.inv = invCompStore;
    results(k1).compTime.jac = jacCompStore;
    results(k1).compTime.iter= iterStore;
    results(k1).residualNorm.PinnKE = resNorm_PinnKE;
    results(k1).residualNorm.LyapKE = resNorm_LyapKE;
    results(k1).residualNorm.PinnLyap = resNorm_PinnLyap;
    results(k1).residualNorm.PinnRK4 = resNorm_PinnRK4;
    results(k1).residualNorm.RK4KE = resNorm_RK4KE;
    
end

options.amountCollocPts = m;
options.amountBasisFuncs = n;
options.timeHorizon = T;
options.freeFuncType = inputStruct.freeFunc;
options.collocType = inputStruct.collocType;
options.intialTau = inputStruct.tau0;
options.finalTau = inputStruct.tauf;
options.amountMcSimulations = mcn;

save('mdre_xtfc_mc.mat', 'results', 'mcn', 'options');





%% Intperoloation Test

% % Interpolution time points
% timeInterp = unifrnd(t00, tff, 50, 1);
% 
% % Constrained expression results at interpolation points
% [ceInterp] = TFC_MDRE_Interp(timeInterp, timeStore, xiStore, ...
%     SfStore, inputStruct);





%% Plots

set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');
%set(groot,'defaultLineLineWidth',2)
set(0,'DefaultAxesFontSize', 20)
%set(0,'DefaultLineMarkerSize',10);

int = floor(linspace(1, length(timeSpanKE), 30));

f1 = figure(1);
p2 = semilogy(timeSpanKE, resNorm_LyapKE, 'r-.'); hold on
for i = 1:mcn
    res(:,i) = results(i).residualNorm.PinnKE;
    %semilogy(timeSpanKE, res(:,i), 'b-'); hold on
end
[b, tf] = rmoutliers(res(1,:), 'percentiles', [5, 95]);
res(:, tf) = [];
meanRes = mean(res, 2);
stdRes = std(res, 1, 2);
ci95 = 1.96*stdRes(int)./sqrt(mcn);
p1 = semilogy(timeSpanKE, meanRes, 'b--'); hold on
e = errorbar(timeSpanKE(int), ...
             meanRes(int), ...
             abs(min(res(int,:)*1, [], 2) - meanRes(int)), abs(max(res(int,:), [], 2) - meanRes(int)), 'b--'); hold on

e.LineStyle = 'none';
p3 = semilogy(timeSpanKE, resNorm_RK4KE, 'k-'); hold off
legend([p1, p2, p3],'PINN', 'LA', 'RK4', 'Location', 'NorthWest');
ylabel('Normwise Error wrt. KE')
xlabel('Time (sec)')
grid
set(f1, 'Position',  [100, 100, 800, 500])






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