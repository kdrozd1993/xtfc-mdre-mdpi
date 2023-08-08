function D = LGL_DiffMatrixB(tauInterp, tauColloc)

n = length(tauColloc);
m = length(tauInterp);

N = length(tauColloc)-1;

D = zeros(n, m);

for i = 1:n
    
    [gColloc, gdColloc, gddColloc] = LeP(tauColloc(i), N);
    
    for j = 1:m
        
        [gInterp, gdInterp, gddInterp] = LeP(tauInterp(j), N);
        
        if i~=j
            D(i,j) = gColloc(end)/(gInterp(end)*(tauColloc(i)-tauInterp(j)));
        elseif i == j && i == 1
            D(i,j) = -N*(N+1)/4;
        elseif i == j && i == n
            D(i,j) = N*(N+1)/4;
        else
            D(i,j) = 0;
        end
        
    end
end

end