clear all
clc

A = [1 3; 3 1];
B1 = [1 2; 4 5];
B2 = [5 1; 3 2];

B = zeros(2,2,2);
B(:,:,1) = B1;
B(:,:,2) = B2;

V = [1; 2];

i = 1;
j = 1;

for k = 1:2
    for q = 1:2    
        M(q,k) = V(q)*B(k,j,q); 
    end
    N(k) = A(k,i)*sum(M(:,k));
end

x1 = sum(N)

x2 = [A(1,1); A(1,2)]'*[V'*reshape(B(1,1,:),2,1); V'*reshape(B(2,1,:),2,1)]