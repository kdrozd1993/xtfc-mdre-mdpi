function D = CGL_DiffMatrixB(tauInterp, tauColloc)

n = length(tauColloc);
m = length(tauInterp);

N = length(tauColloc)-1;

D = zeros(n, m);

for i = 1:n
    
        ci = Cfunc(i, n);
    
    for j = 1:m
        
        cj = Cfunc(j, m);
        
        if i~=j
            D(i,j) = (ci/cj)*(-1)^(j+i)/(tauColloc(i)-tauInterp(j));
        elseif i == j && i == 1
            D(i,j) = -(2*N^2+1)/6;
        elseif i == j && i == n
            D(i,j) = (2*N^2+1)/6;
        else
            D(i,j) = -tauColloc(j)/(2*(1-tauColloc(j)^2));
        end
        
    end
end

end


function c = Cfunc(k, n)

if k == 1 || k == n 
    c = 2;
else
    c = 1;
end

end