%-------------------------------------------------------------------------%
function [xi] = TFC_MDRE_InitialGuess(CE, setup)
%-------------------------------------------------------------------------%

%{

Authors: Kristofer Drozd
Created: 11/10/2021


Description:

This function obtains the initial guess in the XTFC framework for computing 
the Matrix Differential Riccati Equation.


Inputs:

CE - Vectorized constrained expressions of the upper triangular riccati
     matrix

setup - A structure of setup parameters.


Outputs:

xi - The final computed unkown values. (m x amount of symmetric elements).

Version:

11/10/2021 - Initial completion.
12/02/2021 - Adapted to handle the upper triangular elements of the riccati 
             matrix instead of all of them.

%}

n = setup.tfcVariables.n;
m = setup.tfcVariables.m;
nx = setup.tfcVariables.nx;
nCE = 1/2*nx*(nx+1);

Sfvutr = setup.Sfvutr;

H = setup.tfcVariables.H;

xi = LS_Type((H - ones(n,1)*H(end,:)), CE - Sfvutr, 7);

end