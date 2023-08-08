function [xi, iter, lossNorm, CE] = TFC_Riccati_ILS_Debugger(...
    xi, jacFunc, jacCase, maxIter, tol, setup)

err2 = [1e15, 1e15-1];
iter = 0;
nx = setup.tfcVariables.nx;
n = setup.tfcVariables.n;
m = setup.tfcVariables.m;
H = setup.tfcVariables.H;
Sf = setup.Sf;

f1 = figure(1);
CE = zeros(n*nx, nx);
for k = 1:n
    for i = 1:nx
        for j = 1:nx
            indBeg = (((i-1)*nx+j)-1)*m+1;
            indEnd = (((i-1)*nx+j)-1)*m+m;
            CE((k-1)*nx+i, j) = (H(k,:) - H(end,:))*xi(indBeg:indEnd) + Sf(i,j);
        end
    end
end

kk = 1:1:n;
ii = 1;
plot(CE((kk-1)*nx+ii,1), 'r-', 'linewidth', 3); hold on


while abs(err2(2)) > tol && iter < maxIter && abs(err2(1) - err2(2)) > tol

    iter = iter + 1;

    err2(1) = err2(2);
    
    if strcmp(jacCase, 'numeric') == 1
        startJacobianTime = tic;
        J = numeric_jacobian(@Ricatti_Loss, xi, setup);
        computeJacobianTime(iter) = toc(startJacobianTime);
        Lo = Ricatti_Loss(xi, setup);
    elseif strcmp(jacCase, 'adigator') == 1
        startJacobianTime = tic;
        [J, Lo] = jacFunc(xi, setup);
        J = full(J);
        computeJacobianTime(iter) = toc(startJacobianTime);
    end

    startInversionTime = tic;
    [dxi] = LS_Type(J, Lo, 7);
    computeInversionTime(iter) = toc(startInversionTime);

    xi = xi - dxi;

    err2(2) = norm(Lo);
    
    lossNorm(iter) = err2(2);

    CE = zeros(n*nx, nx);
    for k = 1:n
        for i = 1:nx
            for j = 1:nx
                indBeg = (((i-1)*nx+j)-1)*m+1;
                indEnd = (((i-1)*nx+j)-1)*m+m;
                CE((k-1)*nx+i, j) = (H(k,:) - H(end,:))*xi(indBeg:indEnd) + Sf(i,j);
            end
        end
    end
    
    kk = 1:1:n;
    ii = 1;
    plot(CE((kk-1)*nx+ii,1), 'b-'); hold on

end

ylabel('Riccati Matrix Values')
xlabel('Times (sec)')
grid
set(f1, 'Position',  [100, 100, 800, 500])

end