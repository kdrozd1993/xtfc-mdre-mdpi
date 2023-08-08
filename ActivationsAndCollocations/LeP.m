%-------------------------------------------------------------------------%
function [F, FD, FDD, FDDD, FDDDD] = LeP(x, deg)
%-------------------------------------------------------------------------%

% Computes the first deg Legendre orthogonal polynomials (F), with 1st
% (FD), 2nd (FDD), 3rd (FDDD), 4th (FDDDD) derivatives.

x = x(:);
N = length(x);
Zero = zeros(N, 1);
One = ones(N, 1);

% Initialization
if deg == 0
    F = One;
    if nargout > 1, FD   = Zero; end
    if nargout > 2, FDD  = Zero; end
    if nargout > 3, FDDD = Zero; end
    if nargout > 4, FDDDD = Zero; end
    return
elseif deg == 1
    F = [One, x];
    if nargout > 1, FD   = [Zero, One]; end
    if nargout > 2, FDD  = zeros(N, 2); end
    if nargout > 3, FDDD = zeros(N, 2); end
    if nargout > 4, FDDDD = zeros(N, 2); end
else
    F = [One, x, zeros(N, deg-1)];
    if nargout > 1, FD   = [Zero, One, zeros(N, deg-1)]; end
    if nargout > 2, FDD  = zeros(N, deg+1); end
    if nargout > 3, FDDD = zeros(N, deg+1); end
    if nargout > 4, FDDDD = zeros(N, deg+1); end
    
    %% Legendre Polynomials and 1st, 2nd, and 3rd derivatives
    for k = 2:deg
        kk = k-1;
        F(:,k+1) = ((2*kk + 1)*x.*F(:, k) - kk*F(:, k-1))/(kk + 1);
        if nargout > 1, FD(:,k+1) = ((2*kk + 1)*(F(:,k) + x.*FD(:, k)) - kk*FD(:, k-1))/(kk + 1); end
        if nargout > 2, FDD(:,k+1) = ((2*kk + 1)*(2*FD(:,k) + x.*FDD(:, k)) - kk*FDD(:, k-1))/(kk + 1); end
        if nargout > 3, FDDD(:,k+1) = ((2*kk + 1)*(3*FDD(:,k) + x.*FDDD(:, k)) - kk*FDDD(:, k-1))/(kk + 1); end
        if nargout > 4, FDDDD(:,k+1) = ((2*kk + 1)*(4*FDDD(:,k) + x.*FDDDD(:, k)) - kk*FDDDD(:, k-1))/(kk + 1); end
    end
    
end

end