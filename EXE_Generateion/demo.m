%% This is the demo code for 
% Ying Qu and Andreas Koschan and Mongi Abidi. 
% "Robust penalty-weighted deblurring via kernel adaption using single image." 
% Journal of Visual Communication and Image Representation.Volume 41,pp. 109-122,2016.
% The code is for research purpose only, all rights reserved.

% Please cite the following paper 
% @article{QU2016109,
%   title = "Robust penalty-weighted deblurring via kernel adaption using single image",
%   author = "Ying Qu and Andreas Koschan and Mongi Abidi",
%   journal = "Journal of Visual Communication and Image Representation",
%   volume={41},
%   pages={109--122},
%   year={2016}
% }


%% load image
close all;
path = 'input/';
imgname = 'paper.png';
% imgname = 'paper.png';
shortname = imgname(1:(find(imgname=='.')-1));

filepose = strcat(path,imgname);
outputfile ='output/'; 

img = im2double(imread(filepose));

%% deconvolution
% (set PAR.verbose = 0 in parameters.m if you do not want to see the progress)
close all;
hsize = [41 41]; % size of the estimated PSF (upper bound of the expected PSF)
iter = 10;
    timeTemp = tic; %time begin

[u h] = deconvo(img, hsize, '', iter); % this is the main deconvolution routine

%% display result
figure; imshow(u);
figure; imshow(h,[]);
outputname = strcat(outputfile,shortname,'_',num2str(hsize(1)),'_out.png');
imwrite(u,outputname);
outputkenelname = strcat(outputfile,shortname,'_',num2str(hsize(1)),'_kernel.png');
imwrite(linscale(h),outputkenelname);

allTimes = toc(timeTemp); %time end
disp(allTimes);


% hsi = cat(3, th, ts, ti); 
% figure;
% result = hsi2rgb(hsi);



