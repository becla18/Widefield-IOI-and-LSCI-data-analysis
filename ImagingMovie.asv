classdef ImagingMovie < handle
    %IMAGINGMOVIE Summary of this class goes here
    
    properties
        channel %Green, red, blue or speckle
        folderPath %Path to folder containing raw data
        metaDataFile
        nrows
        ncols
        nframes
        freq
        data %data in 3D format (nrows by ncols by nframes)
        stim %stim data
        userSpecified
        processedData
        operationList = {} %List of operations that have already been performed
    end
    
    methods
        %Constructor
        function obj = ImagingMovie(folderPath,channel,movie)
            %IMAGINGMOVIE objects must be created by specifying two arguments.
            %folderPath: string identifying path to the folder containing the raw data
            %channel: can refer to either unprocessed ('green','red','blue','speckle') or processed
            %         ('HbO','HbR') data 
            %Processed data must be specified in 3D array format via the movie argument.            
            obj.channel = channel;
            obj.folderPath = folderPath;
            %Find metadata file. For processed data, the metadata file is chosen as the one from the
            %constitutive unprocessed data.
            if any(strcmp(channel,{'HbO','HbR'}))
                obj.metaDataFile = matfile([folderPath filesep 'Data_green.mat']);
                obj.processedData = 1;
            else
                obj.metaDataFile = matfile([folderPath filesep 'Data_' channel '.mat']);
                obj.processedData = 0;
            end
            obj.nrows = obj.metaDataFile.datSize(1,1);
            obj.ncols = obj.metaDataFile.datSize(1,2);
            obj.nframes = obj.metaDataFile.datLength;
            obj.freq = obj.metaDataFile.Freq;
            obj.stim = obj.metaDataFile.Stim;
            if nargin == 2
                obj.data = getUnprocessedMovie(obj);
                obj.userSpecified = 0;
            elseif nargin == 3
                obj.data = movie;
                obj.userSpecified = 1;
            end               
        end
        
        function movie = getUnprocessedMovie(obj)
            %returns 3D movie matrix of unprocessed data (nrows by ncols by nframes)
            mapFile = memmapfile(obj.metaDataFile.datFile,...
                'Format',{'single',[double(obj.nrows) double(obj.ncols) double(obj.nframes)],'img'},'Repeat',1,'Writable',true);
            movie = mapFile.Data.img;
        end
        function matrix = convertTo2DMatrix(obj,mask)
            %returns 2D data matrix (nframes by npixels) from movie. Optional 2nd argument specifies
            %a binary mask so that only points inside the mask are put in the matrix
            if nargin == 1
                mask = ones(size(obj.data(:,:,1)));
            end
            maskIndexes = find(mask);
            matrix_perm = permute(obj.data,[3 1 2]);
            matrix = reshape(matrix_perm,obj.nframes,obj.nrows*obj.ncols);
            matrix = matrix(:,maskIndexes); %Mask indexes are mapped as the column number in the (time by pixel) time series matrix
        end
        function mask = createMask(obj)
            %creates a user specified 2D mask that can be applied to obj frames
            mask = ROIselect(obj.data(:,:,1));
        end
        function applyMask(obj,mask)
            %applies binary mask to images. If there are NaN points outside the mask, the function
            %converts these to 0.
            obj.data = obj.data.*mask;
            nanInd = find(isnan(obj.data));
            obj.data(nanInd) = 0;
            obj.operationList = [obj.operationList 'applyMask'];
        end
        function correctMotion(obj,translate)
            %corrects motion using dftregistration algorithm. The second argument
            %translate can be specified as a nframes by 2 array, with the first column representing
            %x translations and the second y translations.
            %Frames from the speckle channel are corrected by first applying a 3x3 std kernel to
            %reveal vessel structure (Miao et al., 2010, IEEE Trans. Biomed eng.). 
            if any(strcmp(obj.operationList,'correctMotion'))
                disp('Operation has already been performed');
                return
            end
            if ~exist('transform','var')
                translate = [];
            end
            if strcmp(obj.channel,'speckle')
                movingImages = zeros(size(obj.data));
                for i = 1:obj.nframes
                    movingImages(:,:,i) = stdfilt(obj.data(:,:,i),ones(3));
                end
            else
                movingImages = obj.data;
            end                     
            croppedFrames = cropMovie(movingImages); %The correction is made on user selected ROIs to improve performance by using sharp contrast-rich regions.
            referenceFrame = croppedFrames(:,:,1);
            usfac = 1; %precision level for dftregistration function
            if isempty(translate)
                for j = 1:size(croppedFrames,3)
                    output = dftregistration(fft2(referenceFrame),fft2(croppedFrames(:,:,j)),usfac);
                    tx = output(4); 
                    ty = output(3);
                    obj.data(:,:,j) = imtranslate(obj.data(:,:,j),[tx,ty],'OutputView','same','method','bicubic');
                end
            else
                for j = 1:size(croppedFrames,3)
                    obj.data(:,:,j) = imtranslate(obj.data(:,:,j),translate(j,:),'OutputView','same','method','bicubic');
                end                                        
            end
            obj.operationList = [obj.operationList 'correctMotion'];
        end
        function zscore(obj)
            %Converts each pixel into its z-score
            if any(strcmp(obj.operationList,'zscore'))
                disp('Operation has already been performed');
                return
            end
            obj.data = (obj.data-mean(obj.data,3))./std(obj.data,0,3);
            obj.operationList = [obj.operationList 'zscore'];
        end
        function normalizeByMean(obj)
            if any(strcmp(obj.operationList,'normalizeByMean'))
                disp('Operation has already been performed');
                return
            end
            %Normalizes each pixel by its mean value
            obj.data = obj.data./mean(obj.data,3);
            obj.operationList = [obj.operationList 'normalizeByMean'];
        end
        function gaussianFilter(obj,sigma)
            %Applies gaussian filter with standard deviation sigma to movie frames
            if any(strcmp(obj.operationList,'gaussianFilter'))
                disp('Operation has already been performed');
                return
            end
            for i = 1:obj.nframes
                obj.data(:,:,i) = imgaussfilt(obj.data(:,:,i),sigma);
            end
            obj.operationList = [obj.operationList 'gaussianFilter'];
        end
        function bandPassFilter(obj,cutOnFreq,cutOffFreq)
            if any(strcmp(obj.operationList,'bandPassFilter'))
                disp('Operation has already been performed');
                return
            end
            if cutOffFreq > obj.freq/2
                cutOffFreq = obj.freq/2;
                disp(['Cut off frequency did not respect Nyquist criteria and was reduced to ' num2str(obj.freq) 'Hz']);
            end
            signalMatrix = convertTo2DMatrix(obj);
            [b,a] = cheby1(1,3,[cutOnFreq cutOffFreq]/obj.freq);
            filteredMatrix = filtfilt(b,a,double(signalMatrix)); %same method as Cramer 2019, NeuroImage
            %Convert back to 3D
            reshapedMatrix = reshape(filteredMatrix,obj.nframes,obj.nrows,obj.ncols);
            obj.data = permute(reshapedMatrix,[2,3,1]);
            obj.operationList = [obj.operationList 'bandPassFilter'];
        end
        function convertToSpeckleContrast(obj,kernelSize)
            if any(strcmp(obj.operationList,'convertToSpeckleContrast'))
                disp('Operation has already been performed');
                return
            end
            h = fspecial('average', kernelSize);
            for i = 1:obj.nframes
                img_std = stdfilt(obj.data(:,:,i),ones(kernelSize));
                img_mean = filter2(h,obj.data(:,:,i));
                obj.data(:,:,i) = img_std./img_mean;
            end
            obj.operationList = [obj.operationList 'convertToSpeckleContrast'];
        end
        function averageResponse = averageResponses(obj,preStim,postStim,indexes)
            %UNTITLED Summary of this function goes here
            %   Detailed explanation goes here
            stimOnsetsInd = find(diff(obj.stim)==1)+1;
            stimOffsetsInd = find(diff(obj.stim)==-1);
            nbStim = length(stimOnsetsInd);
            if nargin == 3
                indexes = 1:nbStim;
            end
            stimOnsetsInd = stimOnsetsInd(indexes);
            stimOffsetsInd = stimOffsetsInd(indexes);
            windowLength = length(obj.stim(stimOnsetsInd(1)-round(preStim*obj.freq):stimOffsetsInd(1)+round(postStim*obj.freq)));
            averageResponse = zeros(obj.nrows,obj.ncols,windowLength);
            for i = 1:length(indexes)
                baseline = mean(obj.data(:,:,stimOnsetsInd(i)-round(preStim*obj.freq):stimOnsetsInd(i)-1),3);
                averageResponse = averageResponse+obj.data(:,:,stimOnsetsInd(i)-round(preStim*obj.freq):stimOffsetsInd(i)+round(postStim*obj.freq))./baseline;
            end
            averageResponse = averageResponse/length(indexes);
        end
    end
    
end

