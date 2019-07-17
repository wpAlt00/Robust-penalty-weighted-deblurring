function [H, U, Report] = PSFestimaLnoRgrad(G, iH, PAR, Hstar)
% PSFestim
%
% Estimating PSF and image
% solving u and h-step in FT
% in the h-step image gradient is used

Report = [];

gamma = PAR.gamma;
Lp = PAR.Lp;
ccreltol = PAR.ccreltol;
needtobreak = 0;

% size of H
hsize = [size(iH,1) size(iH,2)];
gsize = [size(G,1), size(G,2)];
usize = gsize;

if PAR.verbose > 1
    FIG_HANDLE_H = figure;
    axes('position',[0.25,0.94,0.5,0.01]);
    axis off; 
    title(['\gamma,\lambda,\alpha,\beta = (',num2str([gamma,PAR.alpha_h,PAR.beta_h]),')']);
    FIG_HANDLE_U = figure;
    axes('position',[0.25,0.94,0.5,0.01]);
    axis off; 
    title(['\gamma,\alpha,\beta = (',num2str([gamma,PAR.alpha_u,PAR.beta_u]),')']);
else
    FIG_HANDLE_H = [];
    FIG_HANDLE_U = [];
end

% if true PSF Hstar is provided -> calculate MSE
if exist('Hstar','var') && ~isempty(Hstar)
    doMSE = 1;
    Report.hstep.mse =  zeros(1,PAR.maxiter+1);
else
    doMSE = 0;
end

U = zeros(usize);
H = iH;

% Initialization of variables for min_U step, which do not change
% FU ... FFT of u
FU = fft2(U);
% FDx, FDx ... FFT of x and y derivative operators
FDx = fft2([1 -1],usize(1),usize(2));
FDy = fft2([1; -1],usize(1),usize(2));
DTD = conj(FDx).*FDx + conj(FDy).*FDy;

% FUx, FUy ... FFT of x (y) derivatives of U
FUx = zeros(usize);
FUy = zeros(usize);

% auxiliary variables for image gradient and blurs
% initialize to zeros
Vx = zeros(usize);
Vy = zeros(usize);
Vh = zeros([usize]);
% extra variables for Bregman iterations
Bx = zeros(usize);
By = zeros(usize);
Bh = zeros([usize]);

if doMSE
	Report.hstep.mse(1) = calculateMSE(H,Hstar);
end

eG = edgetaper(G,ones(hsize)/prod(hsize));

FeGu = fft2(eG);
FeGx = FDx.*FeGu;
FeGy = FDy.*FeGu;

% main loop which alternates between u-estimation and h-estimation
% if PAR.islaststep
%     PAR.maxiter = 10;
% end
relconH = zeros(1,2);

for mI = 1:PAR.maxiter
%     counth = 0;
    % u-estimation
    Ustep;
   
	% h-estimation
    Hstep;
    
    
	% reporting
	if doMSE
        Report.hstep.mse(mI+1) = calculateMSE(H,Hstar);
    end
    
    if needtobreak
        break;
    end
     

    % increasing gamma helps
    gamma = gamma*1.5;
end

% PSF centering
if(PAR.centering_threshold > 0)
	H = centerPSF(H, PAR.centering_threshold);
end

% =========================
% U-STEP
% =========================
function Ustep
	FHS = fft2(H,usize(1),usize(2)); % FT of H
	FHTH = conj(FHS).*FHS; % FT of (H^T)H
	FGs = sum(conj(FHS).*FeGu,3); % FT of (H^T)g (Note that we use edgetaper to reduce border effect)

	beta = PAR.beta_u;
	alpha = PAR.alpha_u;

	% main iteration loop, do everything in the FT domain
	for i = 1:PAR.maxiter_u
		FUp = FU;
		b = FGs + beta/gamma*(conj(FDx).*fft2(Vx+Bx) + conj(FDy).*fft2(Vy+By));
		FU = b./( FHTH + beta/gamma*DTD);

		% Prepare my Lp prior
		Pr = asetupLnormPrior(Lp,alpha,beta);
		FUx = FDx.*FU;
		FUy = FDy.*FU;
		xD = real(ifft2(FUx));
		yD = real(ifft2(FUy));
		xDm = xD - Bx;
		yDm = yD - By;
		nDm = sqrt(xDm.^2 + yDm.^2);
		Vy = Pr.fh(yDm,nDm);
		Vx = Pr.fh(xDm,nDm);
		% update Bregman variables
		Bx = Bx + Vx - xD;
		By = By + Vy - yD;

		E = sqrt(Vy.^2+Vx.^2);
		updateFig(FIG_HANDLE_U,{[] By},{i, [], E},{[] Bx});
		% we can beta after every iteration
		% it should help converegence but probably not necessary
		%beta = 2*beta;

		% Calculate relative convergence criterion
		relcon = sqrt(sum(abs(FUp(:)-FU(:)).^2))/sqrt(sum(abs(FU(:)).^2));

		if relcon < ccreltol
			break;
		end
	end

	if PAR.verbose
		disp(['min_U steps: ',num2str(i),' relcon:',num2str([relcon])]);
	end
	U = real(ifft2(FU));
    

	updateFig(FIG_HANDLE_U,{[],U<0},{i, U, E});
end % end of Ustep

% =======================
% min_H STEP
% =======================


function Hstep
    
    FUD = FeGx.*conj(FUx) + FeGy.*conj(FUy); % basically FT of (U^T)g but calculated in image derivatives
	FUTU = conj(FUx).*FUx + conj(FUy).*FUy; % FT of (U^T)U in image derivatives
	FH = fft2(H,usize(1),usize(2));
	
	beta = PAR.beta_h;
	alpha = PAR.alpha_h;
	
%     beta = 1000;
    
	% main loop
	for i = 1:PAR.maxiter_h       
%         beta = beta*10;
		FHp = FH; 
		b = beta/gamma*fft2(Vh+Bh) + FUD;
		FH = b./(FUTU + beta/gamma);
		
		% Calculate relative convergence criterion
		relconH(1) = sqrt(sum(abs(FHp(:)-FH(:)).^2))/sqrt(sum(abs(FH(:)).^2));
		
		Pr = asetupLnormPrior(1,alpha,beta);
		hI = real(ifft2(FH));
		hIm = hI - Bh;
		nIm = abs(hIm);
		Vh = Pr.fh(hIm,nIm); %Lp norm regularization
		Vh(Vh<0) = 0; % Forcing positivity this way is the correct approach!!!
		% force zero values on h ouside its support 
		Vh(hsize(1)+1:end,:,:) = 0; Vh(1:hsize(1),hsize(2)+1:end,:) = 0;
		% update Bregman variables
		Bh = Bh + Vh - hI;

		H = hI(1:hsize(1),1:hsize(2),:); % new H estimate

        if PAR.islaststep&&(i<PAR.maxiter_h)
            H = anisotropic_gaussian(double(H),1);
        end
        
        if PAR.islaststep&&(i==PAR.maxiter_h)

            %try to define the weight;
            num=0;
            [m,n]=size(H);      
            for i = 1: m
                for j = 1:n
                    if H(i,j)~=0
                       num = num +1;
                       temp2(num,1) = i;
                       temp2(num,2) = j;
                    end
                end
            end

            [a,b] = size(temp2);
            for i = 1: m
                for j = 1:n
                    if H(i,j)>0
                        distance = 0;        
                        for k = 1:a
                            distance = distance + abs(i-temp2(k,1))^2+ abs(j-temp2(k,2)^2);
                        end
                        weight2(i,j) = distance;
                    else
                        weight2(i,j)=0;
                    end

                end
            end

            weight2 = 1./weight2;
            weight2(weight2 == inf)=0;
            weight2 = weight2/sum(weight2(:));

            H = H.*(weight2);
            H = H./sum(H(:));
        end

         

		E = abs(Vh);
		updateFig(FIG_HANDLE_H, {[] Bh}, ...
			{i, H, E },...
			[]);

        % convergence test
        if relconH(1) < ccreltol
            break;
        end

    end
    
	if PAR.verbose
		disp(['min_H step ',num2str(i),' relcon:',num2str([relconH(1)])]);
    end
    
    if relconH(2) ~= 0
        if (relconH(1) - relconH(2))/relconH(1) > 0.2
            needtobreak = 1;
        end 
    end 
    
    relconH(2) = relconH(1);

    
end % end of Hstep
end % end of main function

function r = calculateMSE(h,hs)
    hsize = size(hs);
    i = size(h)-hsize+1;
	h = h/sum(h(:))*sum(hs(:));
	R = im2col(h,[size(hs,1) size(hs,2)],'sliding');
    
    s = sqrt(sum((R-repmat(reshape(hs,prod(hsize(1:2)),1),1,prod(i(1:2)))).^2,1));
    r = s(ceil(prod(i(1:2))/2));
end

function [H] = centerPSF(H, thresh)
%CENTERPSF Thresholds and centers nonzero PSF elements within the PSF window
%
% H				a x b x P matrix of PSFs for P images, P=1 for single channel case
% thresh		threshold in the range [0,1] for determining the PSF support (and for actual thresholding, if any). H is constrast-stretched to [0,1] before applying the threshold.
% method		(temporarily removed) thresholding method
%				'none'		(default) PSF is preserved as is, thresholding is used only to determine PSF support
%				'hard'		values < threshold are set to zero
%
% H				a x b x P matrix, each PSF is centered and sums to 1.

hsize = [size(H,1) size(H,2)];

for i=1:size(H,3)
	h = mat2gray(H(:,:,i)); % contrast-stretching to [0,1]	
	
	% nonzero mask
	m = h >= thresh;
	m2 = bwmorph(m, 'clean'); % remove isolated pixels
	if(any(m2(:))) m = m2; end % preserve delta-like PSFs
	
	% determine mask support
	sum1 = sum(m, 1);
	sum2 = sum(m, 2);
	L = [find(sum2, 1, 'first') find(sum1, 1, 'first')];
	R = [find(sum2, 1, 'last') find(sum1, 1, 'last')];
	topleft = fix((L+R+1-hsize)/2); % topleft corner index of the new mask
	
	% indexing (=shifting)
	h = h(max(topleft(1),1):min(topleft(1)+hsize(1)-1,end), max(topleft(2),1):min(topleft(2)+hsize(2)-1,end)); % get the 'existing' data, then pad borders with zeros
	h = padarray(h, max(topleft-[1,1],0), 0, 'post'); % pad with zeros to end up with the same size
	h = padarray(h, max([1,1]-topleft,0), 0, 'pre');
	
	% normalize sum to 1
	H(:,:,i) = h/sum(h(:));
end
end

