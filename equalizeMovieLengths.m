function [movie1,movie2,varargout] = equalizeMovieLengths(movie1,movie2,varargin)
%EQUALIZEMOVIELENGTHS Verifies if movies under 3D array format have the same number of frames.
%If not, the function removes the last frames from the longer movie so that the rest are equal
%length. More than 2 movies can be entered. The use of this function is hopefully temporary and arises due to errors in data
%aquisition.
nextra = length(varargin);
nframes = zeros(1,2+nextra);    
if strcmp('ImagingMovie',class(movie1)) %CASE 1: ImagingMovie objects as arguments
    nframes(1,1) = movie1.nframes;
    nframes(1,2) = movie2.nframes;
    for i = 1:nextra
        nframes(1,i+2) = varargin{i}.nframes;
    end
    if ~isequal(ones(size(nframes))*nframes(1),nframes)
        minLength = min(nframes);
        movie1.data = movie1.data(:,:,1:minLength);
        movie2.data = movie2.data(:,:,1:minLength);
        movie1.nframes = minLength;
        movie2.nframes = minLength;
        for i = 1:nextra
            varargin{i}.data = varargin{i}.data(:,:,1:minLength);
            varargin{i}.nframes = minLength;
        end
    end
else                                   %CASE 2: Data array as arguments
    if ~isequal(ones(size(nframes))*size(movie1,3),nframes)
        minLength = min(nframes);
        movie1 = movie1(:,:,1:minLength);
        movie2 = movie2(:,:,1:minLength);
        for i = 1:nextra
            varargin{i} = varargin{i}(:,:,1:minLength);
        end
    end
end
varargout = varargin;
end

