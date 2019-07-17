function  T = tensor_reconstruct(e1,e2,l1,l2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
T = zeros( [size(l1),2,2] );
T(:,:,1,1) = l1.*e1(:,:,1).^2 + l2.*e2(:,:,1).^2;
T(:,:,1,2) = l1.*e1(:,:,1).*e1(:,:,2) + l2.*e2(:,:,1).*e2(:,:,2);
T(:,:,2,1) = T(:,:,1,2);
T(:,:,2,2) = l1.*e1(:,:,2).^2 + l2.*e2(:,:,2).^2;

end

