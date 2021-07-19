function [ImMatrix,ROI] = extractImMatrix(path,channel,crop)
%extractImMatrix Extract the image matrix from an IOI/LSCI experiment .dat file
%   ImMatrix is the n x m matrix of the images time series, where n is the number of time points and
%   m the number of pixels.
%   path is the path to the folder that contains the .dat files, first openned from .bin files with Open_IOI_NewSyst.
%   channel indicates which images between 'green','red' or 'speckle' to extract
%   crop (default 0) can be set to 1 to select a sub-region to be analysed with an interractive ROI selector
%   The ROI can be returned for future remapping onto the original image

if nargin == 2
    crop = 0;
end
if strcmp(channel,'green')
    datMatFile = matfile([path filesep 'Data_green.mat']);
elseif strcmp(channel,'red')
    datMatFile = matfile([path filesep 'Data_red.mat']);
elseif strcmp(channel,'speckle')
    datMatFile = matfile([path filesep 'Data_speckle.mat']);
else
    f = errordlg('Channel not found','File Error');
    return
end    
nrows = datMatFile.datSize(1,1);
ncols = datMatFile.datSize(1,2);
nframes = datMatFile.datLength;
mapFile = memmapfile(datMatFile.datFile,...
    'Format',{'single',[double(nrows) double(ncols) double(nframes)],'img'},'Repeat',1,'Writable',true);
threeDimImMatrix = mapFile.Data.img; %Images in 3D format (2D images x time)
%Convert speckle images to speckle contrast
if strcmp(channel,'speckle')
    threeDimImMatrix = Speckle_contrast_conversion(threeDimImMatrix,5);
end
%Format to time x pixels
ImMatrix_perm = permute(threeDimImMatrix,[3 1 2]);
ImMatrix = reshape(ImMatrix_perm,nframes,nrows*ncols);
%Select image sub-region
if crop
    ROI = ROIselect(threeDimImMatrix(:,:,1));
    ImMatrix = ImMatrix(:,ROI(:)'); 
end
end

