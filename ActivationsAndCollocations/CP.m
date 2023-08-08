%-------------------------------------------------------------------------%
function [F, FD, FDD, FDDD, FDDDD] = CP(x, deg)
%-------------------------------------------------------------------------%

% Computes the first deg Chebyshev Orthogonal Polynomials of the fist kind
% (F), with 1st (FD), 2nd (FDD), 3rd (FDDD), 4th (FDDDD) derivatives.

x = x(:);
N = length(x);
Zero = zeros(N, 1);
One = ones(N, 1);

% Initialization
if deg == 0
    F = ones(N, 1);
    if nargout > 1, FD = Zero; end
    if nargout > 2, FDD = Zero; end
    if nargout > 3, FDDD = Zero; end
    if nargout > 4, FDDDD = Zero; end
    return
elseif deg == 1
    F = [One, x];
    if nargout > 1, FD = [Zero, One]; end
    if nargout > 2, FDD = zeros(N, 2); end
    if nargout > 3, FDDD = zeros(N, 2); end
    if nargout > 4, FDDDD = zeros(N, 2); end
else
    F = [One, x, zeros(N, deg-1)];
    if nargout > 1, FD = [Zero, One, zeros(N, deg-1)]; end
    if nargout > 2, FDD = zeros(N, deg+1); end
    if nargout > 3, FDDD = zeros(N, deg+1); end
    if nargout > 4, FDDDD = zeros(N, deg+1); end
    
    % Chebyshev Polynomials with 1st and 2nd derivatives
    for k = 3:deg+1
        F(:,k) = 2*x.*F(:, k-1) - F(:, k-2);
        if nargout > 1, FD(:,k) = 2*(F(:, k-1) + x.*FD(:, k-1)) - FD(:, k-2); end
        if nargout > 2, FDD(:,k) = 4*FD(:, k-1) + 2*x.*FDD(:, k-1) - FDD(:, k-2); end
        if nargout > 3, FDDD(:,k) = 6*FDD(:, k-1) + 2*x.*FDDD(:, k-1) - FDDD(:, k-2); end
        if nargout > 4, FDDDD(:,k) = 8*FDDD(:, k-1) + 2*x.*FDDDD(:, k-1) - FDDDD(:, k-2); end
    end
    
end

end