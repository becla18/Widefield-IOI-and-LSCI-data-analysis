function [Cropped_images, rect] = cropMovie(Images, Rectangle)
%Crop image stack specified as 3D array, nrows x ncols x nframes.
%Rectangle selected on first image is applied to whole stack. If rectangle
%is specified in input, no selection has to be made.

nFrames = size(Images,3);

if nargin == 1
    fig = uifigure;
    uialert(fig,'Select region to be used for motion correction','Motion correction','Icon','info');
    [img, rect] = imcrop(mat2gray(Images(:,:,1)));
    close gcf;
    close(fig);
    Cropped_images = zeros(size(img,1),size(img,2),nFrames);
    for i = 1:nFrames
        Cropped_images(:,:,i) = imcrop(Images(:,:,i),rect);
    end
else
    for i = 1:nFrames
        img = imcrop(Images(:,:,i), Rectangle);
        Cropped_images(:,:,i) = img;
    end
end
end

