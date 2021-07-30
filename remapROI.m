function I = remapROI(Pixels,ROI)
%remapROI takes a vector of pixels and maps them back to their original locations in 2D image
%format.
%%% INPUTS %%%
% Pixels: vector of the m pixels to be mapped.
% ROI: Binary mask of the size of the original 2D image.
%%% OUTPUT %%%
% M: Matrix of remmapped pixels, size of ROI. Outside ROI elements remain 0.
ImageVector = zeros(size(ROI(:)));
ImageVector(ROI(:)) = Pixels;
I = reshape(ImageVector,size(ROI));
end

