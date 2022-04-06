function [movie1,movie2] = equalizeMovieLengths(movie1,movie2)
%EQUALIZEMOVIELENGTHS Verifies if two movies under 3D array format have the same number of frames.
%If not, the function removes the last frames from the longer movie so that the two are equal
%length. The use of this function is hopefully temporary and arises due to errors in data
%aquisition.
if strcmp('ImagingMovie',class(movie1)) %CASE 1: ImagingMovie objects as arguments
    if movie1.nframes ~= movie2.nframes
        minLength = min([movie1.nframes,movie2.nframes]);
        movie1.data = movie1.data(:,:,1:minLength);
        movie2.data = movie2.data(:,:,1:minLength);
        movie1.nframes = minLength;
        movie2.nframes = minLength;
    end
else                                   %CASE 2: Data array as arguments
    if size(movie1,3) ~= size(movie2,3)
        minLength = min([size(movie1,3),size(movie2,3)]);
        movie1 = movie1(:,:,1:minLength);
        movie2 = movie2(:,:,1:minLength);
    end
end
end

