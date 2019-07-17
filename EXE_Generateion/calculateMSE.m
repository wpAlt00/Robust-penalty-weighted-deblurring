function r = calculateMSE(h,hs)
    hsize = size(hs);
    i = size(h)-hsize+1;
	h = h/sum(h(:))*sum(hs(:));
	R = im2col(h,[size(hs,1) size(hs,2)],'sliding');
    
    s = sqrt(sum((R-repmat(reshape(hs,prod(hsize(1:2)),1),1,prod(i(1:2)))).^2,1));
    r = s(ceil(prod(i(1:2))/2));
end