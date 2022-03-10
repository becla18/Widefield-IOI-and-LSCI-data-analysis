classdef FunctionalImagingMovie < handle
    %FUNCTIONALIMAGINGMATRIX Summary of this class goes here
    %   The object must be created by specifying 2 arguments:
    %   folderPath: string identifying path to the folder containing the raw data
    %   channel: can be either 'green','red','blue' or 'speckle'
    
    properties
        channel %Green, red, blue or speckle
        folderPath %Path to folder containing raw data
        metaDataFile
        nrows
        ncols
        nframes
        freq
        data %data in 3D format (nrows by ncols by nframes) 
    end
    
    methods
        function obj = FunctionalImagingMovie(folderPath,channel)
            %FUNCTIONALIMAGDATA Construct an instance of this class
            %   Detailed explanation goes here
            obj.channel = channel;
            obj.folderPath = folderPath;
            obj.metaDataFile = matfile([folderPath filesep 'Data_' channel '.mat']);
            obj.nrows = obj.metaDataFile.datSize(1,1);
            obj.ncols = obj.metaDataFile.datSize(1,2);
            obj.nframes = obj.metaDataFile.datLength;
            obj.freq = obj.metaDataFile.Freq;
            obj.data = getMovie(obj);
        end
        
        function movie = getMovie(obj)
            %getMovie returns 3D movie matrix (nrows by ncols by nframes)
            mapFile = memmapfile(obj.metaDataFile.datFile,...
                'Format',{'single',[double(obj.nrows) double(obj.ncols) double(obj.nframes)],'img'},'Repeat',1,'Writable',true);
            movie = mapFile.Data.img;
        end
        function matrix = get2DMatrix(obj)
            %get2DMatrix returns 2D data matrix (nframes by npixels)
            mapFile = memmapfile(obj.metaDataFile.datFile,...
                'Format',{'single',[double(obj.nrows) double(obj.ncols) double(obj.nframes)],'img'},'Repeat',1,'Writable',true);
            movie = mapFile.Data.img;
            matrix_perm = permute(threeDimImMatrix,[3 1 2]);
            matrix = reshape(matrix_perm,obj.nframes,obj.nrows*obj.ncols); 
        end
        function mask = createMask(obj)
            %createMask creates a user specified 2D mask that can be applied to obj frames
            mask = ROIselect(obj.data(:,:,1));
        end
        function applyMask(obj,mask)
            %applyMask multi
            obj.data = obj.data.*mask;
        end
        function correctMotion(obj,translate)
            %correctMotion corrects motion using dftregistration algorithm. The second argument
            %transform can be specified as a nframes by 2 array, with the first column representing
            %x translations and the second y translations.
             if ~exist('transform','var')
                 translate = [];
             end
             croppedFrames = cropMovie(obj.data); %The correction is made on user selected ROIs to improve performance by using sharp contrast-rich regions.
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
        end
        function zscore(obj)
            %Converts movie data into zscores
            obj.data = (obj.data-mean(obj.data,3))./std(obj.data,0,3);
        end        
            
    end
end

