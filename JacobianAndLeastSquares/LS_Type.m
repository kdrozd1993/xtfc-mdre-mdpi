function [x] = LS_Type(A, b, type)

%--------------------------------------------------------------------------
% USAGE: [x, CN] = LS_Type(A, b, type)
%--------------------------------------------------------------------------
% SYNOPSIS: Solves the least-squares problem, A*x = b, by various methods
%--------------------------------------------------------------------------
% INPUT:
%    A = Coefficient matrix [nxm, where n>m]
%    b = Known vector [nx1]
% type = 1  -->  x = (A'*A)\A'*b;
% type = 2, -->  x = LS_Scaled(A, b);
% type = 3, -->  x = LS_QR(A, b);
% type = 4, -->  x = LS_SVD(A, b);
% type = 5, -->  x = LS_Cholesky(A, b);
% type = 6, -->  scaled QR
% type = 7, -->  lsqminnorm
%--------------------------------------------------------------------------
% OUTPUT:
%    x = Least-squares solution [mx1]
%   CN = Condition number
%--------------------------------------------------------------------------
%      AUTHOR: Daniele Mortari, Texas A&A University
% LAST UPDATE: December 1, 2016
%--------------------------------------------------------------------------
% WARNING: This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------

%% Least-squares type: 1=Normal, 2=Scaled, 3=QR, 4=SVD, 5=Cholesky, 6=Scaled QR, 7=lsqminnorm
switch type
    case 1, x = (A'*A)\A'*b; 
        CN = cond(A'*A); 
    case 2, [x, CN] = LS_Scaled(A, b);
    case 3, [x, CN] = LS_QR(A, b);
    case 4, [U, S, V] = svd(A);
            s = diag(S);
            [n, m] = size(A); SPI = zeros(m, n);
            tol = max(size(A))*eps(norm(s,inf));
            SingVal = zeros(m, 1);
            for k = 1:m
                if abs(s(k)) > tol
                    SPI(k,k) = 1/s(k); SingVal(k) = s(k); 
                end
            end
            %CN = max(SingVal)/min(SingVal);
            x = V*SPI*U'*b;
    case 5, [x] = LS_Cholesky(A, b);
    case 6, S = 1./sqrt(sum(A.*A, 1));
            [Q, R] = qr(A*diag(S));
            x = S(:).*(R\Q'*b);
            %CN = cond(R);
    case 7, x = lsqminnorm(A, b);
            %CN = cond(A);
    case 8, ss = max(size(A)); 
            co = .0001;
            g = A'*inv(co*eye(ss)+A*A');
            x = g* b;
            %CN = cond(co*eye(ss)+A*A');
    otherwise
        error('From LS_Type: wrong input value "type" = [1,7]');
end

end

%-------------------------------------------------------------------------%
function [x, CN] = LS_QR(A,b)

[Q, R] = qr(A);
x = R\Q'*b;
CN = cond(R);

end

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
function [x, CN] = LS_Scaled(A,b)

S = diag(1./sqrt(sum(A.*A, 1)));
B = A*S;
n = (B'*B)\(B'*b);
x = S*n;
CN = cond(B'*B);

end

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
function [x] = LS_Cholesky(A,b)

% Compute the Cholesky decomposition
L = chol(A'*A); % L is upper triangular and (A'A == L'*L)
d = A'*b;
n = length(d);
% Solve L'z = d (via backwards substitution on a lower triangular matrix)
z = backward_subtitution_assuming_lower_triangular(L',d,n);
% Solve Lx=z (via backwards substitution on an upper triangular matrix)
x = backward_subtitution_assuming_upper_triangular(L,z,n);
CN = cond(L);

end

function x = backward_subtitution_assuming_upper_triangular(T,y,n)
x=zeros(n,1);
x(n)=y(n)/T(n,n);
for i=n-1:-1:1
    x(i)=(y(i)-T(i,i+1:n)*x(i+1:n))/T(i,i);
end

end

function x = backward_subtitution_assuming_lower_triangular(T,y,n)
x=zeros(n,1);
x(1)=y(1)/T(1,1);
for i=2:n
    x(i)=(y(i)-T(i,1:i-1)*x(1:i-1))/T(i,i);
end

end

%-------------------------------------------------------------------------%
