function jac = numeric_jacobian(f, x, setup)
nxx = setup.tfcVariables.nx;
nCE = 1/2*nxx*(nxx+1);
m = setup.tfcVariables.m;
% Calculate Jacobian of function f at given x
epsilon = 1e-14; 
epsilon_inv = 1/epsilon;
f0 = feval(f, x, setup); % caclulate f0, when no perturbation happens
x = reshape(x, nCE*m, 1);
nx = length(x); % Dimension of the input x;
% Do perturbation
for i = 1 : nx
    x_ = x;
    x_(i) =  x(i) + epsilon;
    x__ = reshape(x_, m, nCE);
    jac(:, i) = (feval(f, x__, setup) - f0) .* epsilon_inv;
end