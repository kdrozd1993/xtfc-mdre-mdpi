%-------------------------------------------------------------------------%
function [ce] = TFC_MDRE_Interp(t, timeIntervals, xi, Sf, inputStruct)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/10/2021


Description:

This function interpolates the constrained expressions from computed unkown
values at split interval collocation point.


Inputs:

t - The time points to be interpolated.

timeIntervals - A vector of time interval boundaries.

xi - The computed unkown values for each interval. Each row represents the
     interval while the columns are the rconcatenated unkowns for each 
     ricatti element.

Sf - The cost matrix for each interval reshaped as vector such that this a
     matrix where each row represents the time interval.

inputStruct - Structure of tfc parameters for running 
              RiccatiBasisCalculator.


Outputs:

ce - The computed constrained expression values at the interpolated time 
     points.


Version:

11/10/2021 - Initial completion.
12/02/2021 - Adapted to handle the upper triangular elements of the new 
             the new constrained expression formular that is vectorized.

%}

nx = sqrt(size(Sf,2));
ce = zeros(length(t),nx*nx);
m = size(xi,2)/(nx*nx);
tau0 = inputStruct.tau0;
tauf = inputStruct.tauf;

for i = 1:length(t)

    greaterCheck = find(timeIntervals>=t(i));
    lessCheck = find(timeIntervals<=t(i));

    if ~any(timeIntervals==t(i))
        tt = greaterCheck(1);
        tb = lessCheck(end);
    else
        if length(greaterCheck) == 1
            tt = greaterCheck(1);
            tb = lessCheck(end-1);
        else
            tt = greaterCheck(2);
            tb = lessCheck(end);
        end
    end
    
    c = (tauf-tau0)/(timeIntervals(tt) - timeIntervals(tb));
    %tLoop = timeIntervals(tb) + (1/c)*(tau - tau(1));
    tauLoop = [(t(i) - timeIntervals(tb))*c + tau0; tauf];
    
    %tau = [-1 + (2/(timeIntervals(tt)-timeIntervals(tb)))*(t(i)-timeIntervals(tb)); 1]; 
    
    [H, ~, ~] = RiccatiBasisCalculator(inputStruct, tauLoop);

    xitrans = reshape(xi(tb,:)', m, nx*nx);
    
    ce(i,:) = (H(1,:) - H(end,:))*xitrans + Sf(tb,:);

end

end