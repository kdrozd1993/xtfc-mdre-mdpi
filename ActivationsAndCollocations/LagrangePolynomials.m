function L = LagrangePolynomials(x, pointx)
%
%LAGRANGE   approx a point-defined function using the Lagrange polynomial interpolation
%
%      LAGRANGE(X,POINTX,POINTY) approx the function definited by the points:
%      P1=(POINTX(1),POINTY(1)), P2=(POINTX(2),POINTY(2)), ..., PN(POINTX(N),POINTY(N))
%      and calculate it in each elements of X
%
%      If POINTX and POINTY have different number of elements the function will return the NaN value
%
%      function wrote by: Calzino
%      7-oct-2001
%

n = length(pointx);
L = ones(n,length(x));

for i = 1:n
    
  for j = 1:n
      
     if (i~=j)
        L(i,:) = L(i,:).*((x-pointx(j))/(pointx(i)-pointx(j)))';
     end
     
  end
  
end

end