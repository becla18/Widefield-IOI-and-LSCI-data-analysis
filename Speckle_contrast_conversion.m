function [contrast_frames] = Speckle_contrast_conversion(speckle_frames,speckleFilter)
%SPECKLE_CONTRAST_CONVERSION converts speckle images acquired with laser illumination to speckle
%contrast maps
n = speckleFilter;
h = fspecial('average', n);
for i = 1:size(speckle_frames,3)
    img = speckle_frames(:,:,i);
    img_std = stdfilt(img,ones(n));
    img_mean = filter2(h,img);
    speckle_frames(:,:,i) = img_std./img_mean;
end
contrast_frames = speckle_frames;


