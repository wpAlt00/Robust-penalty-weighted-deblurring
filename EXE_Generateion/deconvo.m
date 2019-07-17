function [U H report] = deconvo(G,hsize,hstar,iter)
% load parameters
parameters;
PAR.maxiter = iter;
% get the size (rectangular support) of blurs
if ~exist('hsize','var') || isempty(hsize)
    hsize = [5 5]; % default size of PSFs
end
if iscell(hsize)
  hsize = size(hsize);
else
  H = hsize;
end 

% ground-truth PSF (for comparison)
if ~exist('hstar','var')
    hstar  = [];
end

% normalize images
[G, norm_m, norm_v] = simpnormimg(G);

% PSF ESTIMATION
disp('Estimating PSFs...');

if PAR.doPSFEstimation %I add this to test data

% prepare for multiscale process
L = MSlevels; % number of multiscale levels
if (L<1)
    L = 1;
end
%size of the image selection
sr = [ maxROIsize(1)./(2.^(L-1:-1:0).'), maxROIsize(2)./(2.^(L-1:-1:0).')]; % ROI sizes at each scale level

% precalculate ROI (on which PSF is calculated) for each scale level
ROI = cell(1,L);
hstarP = cell(1,L);
ROI{L} = getROI(G,sr(L,:));
if ~isempty(hstar)
    hstarP{L} = hstar;
end
for i = L-1:-1:1
    ROI{i} = imresize(ROI{i+1},0.5);
    if ~isempty(hstar)
        hstarP{i} = imresize(hstarP{i+1},0.5);
    end
end

% initial PSF size and set them to delta functions
hsize = ceil(hsize/2^(L-1));
cen = floor((hsize+1)/2);
hi = zeros(hsize);
hi(cen(1),cen(2)) = 1; % init PSF as delta impulse

% main estimation loop
report.ms = cell(1,L);
for i = 1:L
	if(PAR.verbose > 0) disp(['hsize: ',num2str(size(hi))]); end
    
	hi = hi/sum(hi(:));
    
    if i == L
        PAR.islaststep = 1;
    end 
	[h u report.ms{i}] = PSFestimaLnoRgrad(ROI{i},hi,PAR,hstarP{i}); % estimate PSF at this scale
    
    hi = imresize(h,2,'lanczos3'); % upsample for new scale
%     figure;
%     imshow(u);
end

% impose constraints on h ???
% h(h<9e-4)=0;
H = h;
H(H<0) = 0;
H = h/sum(H(:));


disp('PSF estimation done.'); % End of PSF estimation
end

% NON-BLIND DECONVOLUTION


% handle circular boundary conditions
[pad_sizem pad_sizen] = size(H);
pad_size = max(pad_sizem,pad_sizen);
pad_size = ceil(pad_size/2);

G = padarray(G, [1 1]*pad_size, 'symmetric', 'both');
for a=1:4
  G = edgetaper(G, H);
end


U = fftCGSRaL(G,H,PAR);

disp('Nonblind deconvolution done.'); % End of PSF estimation code

% denormalize the result
U = U*norm_v + norm_m;

U = U(pad_size+1:end-pad_size, pad_size+1:end-pad_size,:);


report.par = PAR;
end
