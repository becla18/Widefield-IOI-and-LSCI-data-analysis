function movie = mapMatrixToMask(X,mask)
%mapMatrixToMask maps the time series matrix X (time x pixels) to the pixels specified by a binary
%mask and returns a 3D movie array.
%   Detailed explanation goes here
if size(X,2) ~= length(find(mask))
    disp('Error: mask does not correspond to the region in matrix');
    return
end
nframes = size(X,1);
movie = zeros(size(mask,1),size(mask,2),nframes);
for i = 1:nframes
    frame = movie(:,:,i);
    frame(logical(mask(:))) = X(i,:);
    movie(:,:,i) = frame;
end
end

