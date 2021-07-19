function ROI = ROIselect(image,varargin)
%calculates ROI mask from image, with optional 2nd argument as the limits of imagesc
width = size(image,1); heigth = size(image,2);

%%%ROI selection
fig = figure;
if nargin == 2
    limits = varargin{1};
    imagesc(image,limits);
else
    imagesc(image);
end
title(['ROI select (double click last point)']);
[xx,yy] = getline('closed');
hold on;
plot(xx,yy,'r','Linewidth',3); 
ROI = poly2mask(xx,yy,width,heigth);
close(fig);
end