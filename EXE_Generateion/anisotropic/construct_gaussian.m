function gkernel = construct_gaussian(T,t,ksize)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% T is the local structure tensor; t is the time; window size is n x n
temp = 1/4*pi*t;
kn = (ksize-1)/2;
[m,n,a,b]=size(T);
gkernel = zeros(m,n,ksize,ksize);

[dx,dy] = meshgrid(-kn:kn);

for i = 1:m
    for j = 1:n
    b = temp * exp(-(T(i,j,1,1)*dx.^2 + T(i,j,2,1)*dy.*dx...
    + T(i,j,1,2)*dx.*dy + T(i,j,2,2)*dy.^2)/(4*t));
    if sum(b(:))==0
       sumb = 1;
    else
       sumb = sum(b(:));
    end
    b = b./sumb;
    gkernel(i,j,:,:)  = b;
%     gkernel(i,j,:,:) = gkernel(i,j,:,:) ./sum(gkernel(i,j,:)); % normalize
    end
end

end

