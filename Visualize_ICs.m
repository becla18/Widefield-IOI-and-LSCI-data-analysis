%Use this script to perform ICA on a single imaging run
% %% Pre-processing
dataFolderPath = 'E:\Data Jeremie\Bilateral Stim and Resting State dataset\M4\M4_rs_10min';
greenMovie = ImagingMovie(dataFolderPath,'green');
greenMovie.correctMotion;
redMovie = ImagingMovie(dataFolderPath,'red');
redMovie.correctMotion;
speckleContrastMovie = ImagingMovie(dataFolderPath,'speckle');
speckleContrastMovie.convertToSpeckleContrast(5);
speckleContrastMovie.binFrames;
speckleContrastMovie.correctMotion;
greenMovie.normalizeByMean;
redMovie.normalizeByMean;
%%
[greenMovie,redMovie] = equalizeMovieLengths(greenMovie,redMovie);
[greenMovie,speckleContrastMovie] = equalizeMovieLengths(greenMovie,speckleContrastMovie);
[HbOMovie, HbRMovie] = convertToHb(greenMovie.data,redMovie.data);
HbOMovie = ImagingMovie(dataFolderPath,'HbO',HbOMovie);
HbRMovie = ImagingMovie(dataFolderPath,'HbR',HbRMovie);
clear greenMovie redMovie
mask = HbOMovie.createMask;
if speckleContrastMovie.nrows ~= HbOMovie.nrows
    maskSC = imresize(mask,[speckleContrastMovie.nrows speckleContrastMovie.ncols]);
else
    maskSC = mask;
end
HbOMovie.applyMask(mask);
HbRMovie.applyMask(mask);
speckleContrastMovie.applyMask(maskSC);
HbOMovie.bandPassFilter(0.01,0.1);
HbRMovie.bandPassFilter(0.01,0.1);
speckleContrastMovie.bandPassFilter(0.01,0.1);
HbOMatrix = HbOMovie.convertTo2DMatrix(mask);
HbRMatrix = HbRMovie.convertTo2DMatrix(mask);
SCMatrix = speckleContrastMovie.convertTo2DMatrix(maskSC);
% %% ICA
% n = 10;
% [componentsHbO, wHbO, tHbO, muHbO] = fastICA(HbOMatrix,n,'negentropy');
% [componentsHbR, wHbR, tHbR, muHbR] = fastICA(HbRMatrix,n,'negentropy');
% [componentsSC, wSC, tSC, muSC] = fastICA(SCMatrix,n,'negentropy');
% A_HbO = tHbO\wHbO'; %A is the matrix that holds the mixing weights
% A_HbO = (A_HbO-mean(A_HbO))./std(A_HbO);
% A_HbR = tHbO\wHbR';
% A_HbR = (A_HbR-mean(A_HbR))./std(A_HbR);
% A_SC = tSC\wSC';
% A_SC = (A_SC-mean(A_SC))./std(A_SC);
% 
% %Display first n ICs
% figure('Color',[1 1 1],'Position',[1531,69,1446,516]);
% map = zeros(size(mask));
% t = tiledlayout(ceil(n/5),5);
% t.TileSpacing = 'compact';
% t.Padding = 'none';
% t.Title.String = 'HbO ICs (0.01 - 0.15 Hz)';
% for i = 1:n
%     nexttile;
%     map(logical(mask(:)))=componentsHbO(i,:);
%     imagesc(map);
%     title(['IC ' num2str(i)]);
%     colorbar;
%     colormap hot;
%     axis image;
%     xticks([]);
%     yticks([]);
% end
% 
% figure('Color',[1 1 1],'Position',[1531,69,1446,516]);
% map = zeros(size(mask));
% t = tiledlayout(ceil(n/5),5);
% t.TileSpacing = 'none';
% t.Padding = 'none';
% t.Title.String = 'HbR ICs (0.01 - 0.15 Hz)';
% for i = 1:n
%     nexttile;    
%     map(find(mask))=componentsHbR(i,:);
%     imagesc(map);
%     title(['IC ' num2str(i)]);
%     colorbar;
%     colormap hot;
%     axis image;
%     xticks([]);
%     yticks([]);
% end
% 
% figure('Color',[1 1 1],'Position',[1531,69,1446,516]);
% map = zeros(size(maskSC));
% t = tiledlayout(ceil(n/5),5);
% t.TileSpacing = 'none';
% t.Padding = 'none';
% t.Title.String = 'Speckle contrast ICs (0.01 - 0.15 Hz)';
% for i = 1:n
%     nexttile;    
%     map(find(maskSC))=componentsSC(i,:);
%     imagesc(map);
%     title(['IC ' num2str(i)]);
%     colorbar;
%     colormap hot;
%     axis image;
%     xticks([]);
%     yticks([]);
% end