function [I m v] = simpnormimg(G)
% Normalize images in G, m is the minimum value in G and v is the span
lb = min(G(:));
ub = max(G(:));

v = ub-lb;
m = lb;
I = (double(G)-m)/v;
end