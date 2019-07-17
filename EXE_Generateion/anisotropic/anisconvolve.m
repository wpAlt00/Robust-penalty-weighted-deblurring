function Iopt = anisconvolve(I,kernel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% not sure why there is nan


[m,n,a,b] = size(kernel);
    pb = floor( a/2 );  % padding before
    pf = a-pb-1;            % padding after
    ipad = I;
    
    %flip the edge value to padding;
    %flip the row value
    ipad = [ ipad(pb:-1:1,:); ipad; ipad(end:-1:end-pf+1,:) ];
    %flip the column value
    ipad = [ ipad(:,pb:-1:1), ipad, ipad(:,end:-1:end-pf+1) ];
    Iopt = ipad;
    
    %apply gaussian
for i = 1:m
    for j = 1:n
       tempk = reshape(kernel(i,j,:,:),a,b);
       tempk(isnan(tempk))= 0;
       if sum(tempk(:))==0
           tempk(1+pb,1+pb) = 1;
       end
       needsum = ipad(i:i+2*pb,j:j+2*pb).*tempk;
       Iopt(i+pb,j+pb) = sum(needsum(:));
    end
end
    Iopt = Iopt( (pb + 1):(pb + m), (pb + 1):(pb + n) );
end

