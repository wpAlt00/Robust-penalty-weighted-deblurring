function [ e1,e2,l1,l2 ] = eigendecomposition( M )
%EIGENDECOMPOSITION Summary of this function goes here
%Decomposite the 2x2 matrix

%   Detailed explanation goes here
[m,n,a,b] = size(M);
e1 = zeros(m,n,2);
e2 = zeros(m,n,2);
l1 = zeros(m,n);
l2 = zeros(m,n);

for i = 1:m
    for j = 1:n
    temp(1,1) = M(i,j,1,1);
    temp(1,2) = M(i,j,1,2);
    temp(2,1) = M(i,j,2,1);
    temp(2,2) = M(i,j,2,2);
    [V,D] = eig(temp);
    %arrange it by lamda1 < lambda2
    e1(i,j,:) = V(:,1);
    e2(i,j,:) = V(:,2);
    l1(i,j) = D(1,1);
    l2(i,j) = D(2,2);
    end 
end 

end

