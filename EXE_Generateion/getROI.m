function R = getROI(G,win)
	% selects center part of image for PSF estimation
	
	isize = size(G);
	gsize = isize(1:2);
	if size(G,3) > 1 % RGB image
		cind = 2; %green channel
	else
		cind = 1;
	end
	% if the window size is larger then the image size, set win=isize.
	if any(gsize < win)
		win = gsize;
	end
	margin = floor((gsize-win)/2);
	R = G(margin(1)+1:margin(1)+win(1),margin(2)+1:margin(2)+win(2),cind);
end
