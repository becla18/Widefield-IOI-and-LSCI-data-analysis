function [Cropped_images, rect] = cropMovie(Images, Rectangle)
%Crop image stack specified as 3D array, nrows x ncols x nframes.
%Rectangle selected on first image is applied to whole stack. If rectangle
%is specified in input, no selection has to be made.

nFrames = size(Images,3);

if nargin == 1
    imObject = imshow(mat2gray(Images(:,:,1)));
    imObject.Parent.Title.String = 'Select region to be used for motion correction (double click when done)';
    imObject.Parent.Parent.Position = [1855,130,570,430];
    [img, rect] = imcrop(imObject);
    close gcf;
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

