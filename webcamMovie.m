classdef webcamMovie < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        path
        videoObject
        syncSignal
        frames
        operationList
    end
    
    methods
        %Constructor
        function obj = webcamMovie(moviePath,syncFilePath)
            %webcamMovie object is built by specifying the path to a video file. The object can
            %access the properties of a VideoReader object through the videoObject property. A
            %second path syncFilePath that leads to a file containing synchronization signals can be
            %specified to assign webcam frames to the frames from an ImagingMovie object.
            obj.path = moviePath;
            obj.videoObject = VideoReader(obj.path);
            if nargin == 2
                load(syncFilePath,'GreenLED'); %We are here still dependent on the precise naming of the LED trigger signal 
                obj.syncSignal = GreenLED; %LED trigger signal
            end
        end
        
        function convertToGrayscaleFrames(obj)
            %Reads the video file and converts it to a 3D array of gray scale frames, stored under
            %the frames property.
            grayFrames = zeros(obj.videoObject.height,obj.videoObject.width,obj.videoObject.NumFrames);
            for i = 1:obj.videoObject.NumFrames
                frame = read(obj.videoObject,i);
                grayFrames(:,:,i) = double(rgb2gray(frame)); %Convert to proper gray scale format
            end
            obj.frames = grayFrames;
        end
        function signal = extractROIsignal(obj,mask)
            %Extract the average signal from the ROI specified a binary mask.
            signal = zeros(obj.videoObject.NumFrames,1);
            npixelROI = sum(sum(mask));
            for i = 1:obj.videoObject.NumFrames
                frame = read(obj.videoObject,i);
                frame = double(rgb2gray(frame)); %Convert to proper gray scale format
                signal(i) = sum(sum(frame.*mask))/npixelROI;
            end
        end         
        function mask = createMask(obj)
            %Create binary mask by hand-selecting a ROI
            frame = readFrame(obj.videoObject);
            mask = ROIselect(frame); %Binary mask
        end
        function signal = synchronizeSignal(obj,initialSignal)
            if isempty(obj.syncSignal)
                disp('Specify path to synchronization signal in object definition');
                return
            end
            threshold = 4.5; %voltage threshold to consider the LED ON
            diffSyncSignal = [0 diff(obj.syncSignal)]; %a zero is added so that the trigger signal is the same length as the original signal, 
            LEDtrigs = diffSyncSignal>threshold;
            signal = initialSignal(LEDtrigs);
        end
        function zscore(obj)
            %Converts each pixel into its z-score
            if any(strcmp(obj.operationList,'zscore'))
                disp('Operation has already been performed');
                return
            end
            obj.frames = (obj.frames-mean(obj.frames,3))./std(obj.frames,0,3);
            obj.operationList = [obj.operationList 'zscore'];
        end
    end
end

