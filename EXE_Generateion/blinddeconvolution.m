function blinddeconvolution(varargin)
varsize =  size(varargin);
fname = varargin{1};
psfx = str2double(varargin{2});
psfy = str2double(varargin{3});
disp(fname);
disp(psfx);
disp(psfy);

% addpath('toolbox_general');
% addpath('toolbox_signal');
% addpath('anisotropic');
% files = dir(strcat('demo/',fname,'/','*.jpg'));
files = dir(strcat(fname,'/','*.jpg'));
 
if isempty(files)
%     files = dir(strcat('demo/',fname,'/','*.png'));
    files = dir(strcat(fname, '/', '*.png'));
 
end
if isempty(files)
%     files = dir(strcat('demo/',fname,'/','*.bmp'));
    files = dir(strcat(fname, '/', '*.bmp'));
 
end
if isempty(files)
%     files = dir(strcat('demo/',fname,'/','*.tif'));
    files = dir(strcat(fname, '/', '*.tif'));
 
end
% path = strcat('demo/',fname,'/');
[n,m]= size(files)
for i = 1:n
    path = strcat(fname,'/');
    disp(path);
    name = files(i).name;
    disp(name);
    filename = name(1:(find(name=='.')-1));
    image = im2double(imread(strcat(path,name)));
    close all;
    hsize = [psfx psfy];
    % hsize = psf;
    disp(hsize);
    iter = 10;
    timeTemp = tic; %time begin
    [u h] = deconvo(image, hsize, '', iter); % this is the main deconvolution routine
    allTimes = toc(timeTemp); %time end
    %% save result
    output = '/output/';
    imwrite(u,strcat(path,output,name,num2str(psfx), '_',num2str(psfy), '_out.jpg'));
    imwrite(linscale(h),strcat(path,output,name,num2str(psfx), '_',num2str(psfy), '_psf_out.jpg'));
    disp(allTimes);
end
