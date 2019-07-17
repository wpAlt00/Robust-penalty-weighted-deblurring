function h = anisotropic_gaussian(f,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

G = grad(f);
[m,n] = size(f);
T0 = zeros(m,n,2,2);

T0(:,:,1,1) = G(:,:,2).^2;
T0(:,:,2,2) = G(:,:,1).^2;
T0(:,:,1,2) = -G(:,:,1).*G(:,:,2);
T0(:,:,2,1) = -G(:,:,1).*G(:,:,2);

%used to normalize the T0
T1 = zeros(m,n);
T1(:,:) = sqrt(G(:,:,1).^2 + G(:,:,2).^2); 
T1(find(T1==0))=1;

TK = zeros(m,n,2,2); 
TK(:,:,1,1) = T0(:,:,1,1)./T1(:,:);
TK(:,:,1,2) = T0(:,:,1,2)./T1(:,:);
TK(:,:,2,1) = T0(:,:,2,1)./T1(:,:);
TK(:,:,2,2) = T0(:,:,2,2)./T1(:,:);

[e1,e2,lambda1,lambda2] = eigendecomposition(TK);
T = tensor_reconstruct(e1,e2,lambda1,lambda2);
gkernel = construct_gaussian(T,t,3);
h = anisconvolve(f,gkernel);

end

