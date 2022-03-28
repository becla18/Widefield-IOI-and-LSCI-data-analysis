folderPath = 'E:\Data Jeremie\Bilateral Stim and Resting State dataset\M4\M4_rs_10min';
greenMovie = ImagingMovie(folderPath,'green');
greenMovie.correctMotion;
redMovie = ImagingMovie(folderPath,'red');
redMovie.correctMotion;
speckleContrastMovie = ImagingMovie(folderPath,'speckle');
speckleContrastMovie.correctMotion;
greenMovie.normalizeByMean;
redMovie.normalizeByMean;
[HbOMovie, HbRMovie] = convertToHb(greenMovie.data,redMovie.data);
speckleContrastMovie.convertToSpeckleContrast(5);
HbOMovie = ImagingMovie(folderPath,'HbO',HbOMovie);
HbRMovie = ImagingMovie(folderPath,'HbR',HbRMovie);
clear greenMovie redMovie
%%
HbOMovie.gaussianFilter(1);
HbRMovie.gaussianFilter(1);
speckleContrastMovie.gaussianFilter(1);
%%
HbOMovie.zscore;
HbRMovie.zscore;
speckleContrastMovie.zscore;
%%
maskIOI = HbOMovie.createMask;
maskSpeckle = speckleContrastMovie.createMask;
HbOMovie.applyMask(maskIOI);
HbRMovie.applyMask(maskIOI);
speckleContrastMovie.applyMask(maskSpeckle);
%%
HbOMovie.bandPassFilter(0.01,0.15);
HbRMovie.bandPassFilter(0.01,0.15);
speckleContrastMovie.bandPassFilter(0.01,0.15);
%%
HbOTimeSeries = HbOMovie.convertTo2DMatrix(maskIOI);
HbRTimeSeries = HbRMovie.convertTo2DMatrix(maskIOI);
speckleContrastTimeSeries = speckleContrastMovie.convertTo2DMatrix(maskSpeckle);
%%
%With Matlab's pca function, rows must correspond to observations, and columns to variables
%(dimensions). With the fastICA function, this is inversed.
[scoresHbO] = fastICA(HbOTimeSeries,20,'negentropy');
[scoresHbR] = fastICA(HbRTimeSeries,20,'negentropy');
[scoresSpeckle] = fastICA(speckleContrastTimeSeries,20,'negentropy');

figure;
title('HbO PCs');
map = zeros(size(maskIOI));
tiledlayout(4,5)
for i = 1:20
    nexttile;
    map(find(maskIOI))=scoresHbO(i,:);
    imagesc(map);colorbar;
end

figure;
title('HbR PCs');
map = zeros(size(maskIOI));
tiledlayout(4,5)
for i = 1:20
    nexttile;
    map(find(maskIOI))=scoresHbR(i,:);
    imagesc(map);colorbar;
end

figure;
title('Speckle Contrast PCs');
map = zeros(size(maskSpeckle));
tiledlayout(4,5)
for i = 1:20
    nexttile;
    map(find(maskSpeckle))=scoresSpeckle(i,:);
    imagesc(map);colorbar;
end

