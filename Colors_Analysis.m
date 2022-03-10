function Colors_Analysis(FolderName,preStim,postStim,Binning,speckleFilter,NewAnalysis,Registration,CorrectSCMovement,Aberrant_Events,spatialFilt,spatialFiltStd,timeMedFilt,timeFiltRange,SGFilt,SGFiltOrder,LowpassFilter,PassbandFreq,StopbandFreq,Detrending,singleTrials,substractBG,Fluo)
%Computes mean response images for Green, Red, HbO, HbR, fluorescence and Speckle
%contrast

f = uifigure;
d = uiprogressdlg(f,'Title','Processing Data','Indeterminate','on');

load([FolderName filesep 'StimParameters.mat']);

%list of good events
GoodEvents = [];
for i = 1:NbStim
    if isempty(find(i==Aberrant_Events,1))
        GoodEvents = [GoodEvents i];
    end
end
NbStim = NbStim-length(Aberrant_Events);

%% Openning Frames
if NewAnalysis && exist([FolderName filesep 'Mean_Frames.mat']) == 2
    delete([FolderName filesep 'Mean_Frames.mat']);
end

if exist([FolderName filesep 'Mean_Frames.mat']) == 2
    load([FolderName filesep 'Mean_Frames.mat']);
    Mean_Frames_loaded = 1;
else       
    [greenMap, redMap, speckleMap, fluoMap, Ref_Frame] = openFrames(FolderName,Binning);
    Mean_Frames_loaded = 0;
    %Dark image for background substraction. No binning for speckle images.
if substractBG    
    if Binning
        [Background, folderName] = selectBackground(Binning);
        Background_speckle = selectBackground(0,folderName); 
    else
        Background = selectBackground(Binning);
        Background_speckle = Background;
    end
    Dat_Gptr = matfile([FolderName filesep 'Data_green.mat']);
    Dat_Sptr = matfile([FolderName filesep 'Data_speckle.mat']);
    if size(Background,1)~=Dat_Gptr.datSize(1,1) || size(Background,2)~=Dat_Gptr.datSize(1,2) %Correction if background image size does not correspond to other images
        Background = imresize(mat2gray(Background,[0 4095]),Dat_Gptr.datSize);
        Background = Background*4095;
        Background_speckle = imresize(mat2gray(Background,[0 4095]),Dat_Sptr.datSize);
        Background_speckle = Background_speckle*4095;
    end
end
    responseList = {}; %List containg the names of variables to save at the end
    if ~isempty(Ref_Frame)
        responseList = [responseList 'Ref_Frame'];
    end
end

if ~Mean_Frames_loaded
    IsThereGreen = false; 
    IsThereRed = false;
    IsThereSpeckle = false;
    IsThereFluo = false;

    if ~isempty(greenMap)
    IsThereGreen = true;
    end
    if ~isempty(redMap)
        IsThereRed = true;
    end
    if ~isempty(speckleMap)
        IsThereSpeckle = true;
    end
    if ~isempty(fluoMap)
        IsThereFluo = true;
    end
    ncolors = IsThereGreen+IsThereRed+IsThereSpeckle+IsThereFluo;
    
    %% Speckle channel detected %%
    if IsThereSpeckle
        Dat_Sptr = matfile([FolderName filesep 'Data_speckle.mat']);
        freq = Dat_Sptr.Freq;
        StimG = Dat_Sptr.Stim; 
        SC_Frames = speckleMap.Data.img;
        if substractBG
            SC_Frames = SC_Frames - Background_speckle;
        end
        nrows = size(SC_Frames,1);
        ncols = size(SC_Frames,2);
        nframes = size(SC_Frames,3);
        %Indexes of stim onsets
        onsets = find(diff(StimG)==1)+1;
        onsets = onsets(GoodEvents);
        if Detrending                       
            SC_Frames = detrendImageSeries(SC_Frames,2,1);      
        end 
        if Registration
            refImg = stdfilt(SC_Frames(:,:,1),ones(3));
            for i = 1:nframes
                movingImg = stdfilt(SC_Frames(:,:,i),ones(3));
                [output,~] = dftregistration(fft2(refImg),fft2(movingImg),1);
                tx = output(4);
                ty = output(3);
                SC_Frames(:,:,i) = imtranslate(SC_Frames(:,:,i),[tx,ty],'OutputView','same','method','bicubic');
            end
        end
        SC_Frames = Speckle_contrast_conversion(SC_Frames,speckleFilter); %Convert speckle to speckle contrast
        if CorrectSCMovement           
            %Choose ROI for movement regression
            ROI_movement  = ROI_movementRegression(SC_Frames(:,:,1));
            npixelROI = sum(ROI_movement(:));
            %%%Regress movement signal            
            g = zeros(nframes,1); %g is the movement regressor           
            %SC_Frames = SC_Frames - Background_speckle; %Substract background         
            g = squeeze(sum(sum(SC_Frames.*ROI_movement))/npixelROI);
            g = g-mean(g);
            SC_Frames = distributed(SC_Frames); %Parallel computing to accelerate permutation operations
            SC_Frames_p = permute(SC_Frames,[3 1 2]);%Format for regression
            SC_Frames_r = reshape(SC_Frames_p,nframes,nrows*ncols); 
            SC_Frames_reg = regressMovement(SC_Frames_r,g);
            SC_Frames = permute(reshape(SC_Frames_reg,nframes,nrows,ncols),[2 3 1]);%Reformat back to image format
            SC_Frames = gather(SC_Frames);
            clear SC_Frames_p SC_Frames_r SC_Frames_reg speckleMap           
        end        
    
        %Preallocation for frames of response window (pre-stim+stim+post-stim) & baseline frames
        windowSize = floor(freq*(preStim+StimLength+postStim))+1;
        SC_indResponse = zeros(nrows,ncols,windowSize);
        rSC_Response = zeros(nrows,ncols,windowSize); %relative speckle contrast response;
        absSC_Response = zeros(nrows,ncols,windowSize); %absolute speckle contrast response
        preStim_Frames = zeros(nrows,ncols,length(1:floor(preStim*freq)));
        baseline = zeros(nrows,ncols);
        if singleTrials
            SC_Trials = cell(1,NbStim);
        end
        for i = 1:NbStim
            SC_indResponse = SC_Frames(:,:,(onsets(i)-floor(preStim*freq)):(onsets(i)-floor(preStim*freq))+windowSize-1);        
            if spatialFilt
              for j = 1:size(SC_indResponse,3)
                   SC_indResponse(:,:,j) = imgaussfilt(SC_indResponse(:,:,j),spatialFiltStd);
              end 
            end
            if LowpassFilter
                SC_indResponse = reshape(permute(SC_indResponse,[3 1 2]),windowSize,[]);
                SC_indResponse = filtfilt(f,SC_indResponse);
                SC_indResponse = reshape(permute(SC_indResponse,[2 1]),nrows,ncols,windowSize);
            end
            if timeMedFilt
                SC_indResponse = timeMedianFilter(SC_indResponse,timeFiltRange,freq);
            end
            if SGFilt
                SC_indResponse = timeSGFilter(SC_indResponse,SGFiltOrder,timeFiltRange,freq);
            end      
            %Calculate relative response
            preStim_Frames = SC_indResponse(:,:,1:floor(preStim*freq));
            baseline = mean(preStim_Frames,3);
            rSC_Response = rSC_Response+SC_indResponse./baseline;
            absSC_Response = absSC_Response+SC_indResponse;
            if singleTrials
                SC_Trials{i} = SC_indResponse./baseline;
            end
        end
        %Full relative response
        rSC_Response = rSC_Response/NbStim;
        %Full absolute SC response
        absSC_Response = absSC_Response/NbStim;   
        responseList = [responseList 'rSC_Response' 'absSC_Response'];
        clear Dat_Sptr SC_Frames SC_indResponse
    end
        %% Green channel detected %%
   if IsThereGreen
        Dat_Gptr = matfile([FolderName filesep 'Data_green.mat']);
        freq = Dat_Gptr.Freq;
        StimG = Dat_Gptr.Stim;
        nrows = Dat_Gptr.datSize(1,1);
        ncols = Dat_Gptr.datSize(1,2);
        Green_Frames = greenMap.Data.img;
        if substractBG
            Green_Frames = Green_Frames - Background;
        end
        %Indexes of stim onsets
        onsets = find(diff(StimG)==1)+1;
        onsets = onsets(GoodEvents);
        if Detrending
            Green_Frames = detrendImageSeries(Green_Frames,2,1);
        end      
        %Preallocation for frames of response window (pre-stim+stim+post-stim) & baseline frames
        windowSize = floor(freq*(preStim+StimLength+postStim))+1;
        Green_indResponse = zeros(nrows,ncols,windowSize);
        Rel_Green_Response = zeros(size(Green_indResponse));
        preStim_Frames = zeros(nrows,ncols,length(1:floor(preStim*freq)));
        baseline = zeros(nrows,ncols);
        if singleTrials
            Green_Trials = cell(1,NbStim);
        end
        for i = 1:NbStim
            Green_indResponse = Green_Frames(:,:,(onsets(i)-floor(preStim*freq)):(onsets(i)-floor(preStim*freq))+windowSize-1);
            if Registration
                if i == 1
                    [mvmtcorrFrames,rect] = Crop(Green_indResponse); %identify region to use for motion correction computations
                else
                    mvmtcorrFrames = Crop(Green_indResponse,rect);
                end
                ref = mvmtcorrFrames(:,:,1); %reference frame for motion correction
                usfac = 1; %precision level for dftregistration function
                for j = 1:size(mvmtcorrFrames,3)
                    output = dftregistration(fft2(ref),fft2(mvmtcorrFrames(:,:,j)),usfac);
                    tx = output(4); %Keep movement variables in array for possible future use in GLM
                    ty = output(3);
                    Green_indResponse(:,:,j) = imtranslate(Green_indResponse(:,:,j),[tx,ty],'OutputView','same','method','bicubic');
                end
            end
            if spatialFilt
              for j = 1:size(Green_indResponse,3)
                      Green_indResponse(:,:,j) = imgaussfilt(Green_indResponse(:,:,j),spatialFiltStd);
              end 
            end
            if LowpassFilter
                f = designfilt('lowpassfir','PassbandFrequency', PassbandFreq, 'StopbandFrequency', StopbandFreq, ...
                    'StopbandAttenuation', 60,'SampleRate', freq);
                Green_indResponse = reshape(permute(Green_indResponse,[3 1 2]),windowSize,[]);
                Green_indResponse = filtfilt(f,Green_indResponse);
                Green_indResponse = reshape(permute(Green_indResponse,[2 1]),nrows,ncols,windowSize);
            end
            if timeMedFilt
                Green_indResponse = timeMedianFilter(Green_indResponse,timeFiltRange,freq);
            end
            if SGFilt
                Green_indResponse = timeSGFilter(Green_indResponse,SGFiltOrder,timeFiltRange,freq);
            end
            %Calculate relative response
            preStim_Frames = Green_indResponse(:,:,1:floor(preStim*freq));
            baseline = mean(preStim_Frames,3);
            Rel_Green_Response = Rel_Green_Response+Green_indResponse./baseline;
            if singleTrials
                Green_Trials{i} = Green_indResponse./baseline;
            end
        end
        %Full average relative response
        Rel_Green_Response=Rel_Green_Response/NbStim;
        responseList = [responseList 'Rel_Green_Response'];
        %Time series time axis
        timeInc = 1/freq;
        Timeline = -floor(preStim*freq)*timeInc:timeInc:(size(Rel_Green_Response,3)-floor(preStim*freq)-1)*timeInc;
        time_corr = 6*(freq*ncolors)^(-1)/ncolors; %Correction for time offset due to sequential color acquisition
        Timeline = Timeline+time_corr;
        clear Dat_Gptr Green_Frames Green_indResponse
   end   
           %% Fluo channel detected %%
   if IsThereFluo && Fluo
        Dat_Gptr = matfile([FolderName filesep 'Data_green.mat']);
        Dat_Fptr = matfile([FolderName filesep 'Data_yellow.mat']);
        freq = Dat_Gptr.Freq;
        StimG = Dat_Gptr.Stim;
        StimF = Dat_Fptr.Stim;
        nrows = Dat_Gptr.datSize(1,1);
        ncols = Dat_Gptr.datSize(1,2);
        Green_Frames = greenMap.Data.img;
        Fluo_Frames = fluoMap.Data.img;
        if substractBG
            Green_Frames = Green_Frames - Background;
            Fluo_Frames = Fluo_Frames - Background;
        end
        %Indexes of stim onsets
        onsetsG = find(diff(StimG)==1)+1;
        onsetsG = onsets(GoodEvents);
        %Indexes of stim onsets
        onsetsF = find(diff(StimF)==1)+1;
        onsetsF = onsets(GoodEvents);
        if Detrending
            Green_Frames = detrendImageSeries(Green_Frames,2,1);
            Fluo_Frames = detrendImageSeries(Fluo_Frames,2,1);
        end      
        %Preallocation for frames of response window (pre-stim+stim+post-stim) & baseline frames
        windowSize = floor(freq*(preStim+StimLength+postStim))+1;
        Rel_Fluo_Response = zeros(nrows,ncols,windowSize);
        if singleTrials
            Green_Trials = cell(1,NbStim);
            Fluo_Trials = cell(1,NbStim);
        end
        for i = 1:NbStim
            Green_indResponse = Green_Frames(:,:,(onsetsG(i)-floor(preStim*freq)):(onsetsG(i)-floor(preStim*freq))+windowSize-1);
            Fluo_indResponse = Fluo_Frames(:,:,(onsetsF(i)-floor(preStim*freq)):(onsetsF(i)-floor(preStim*freq))+windowSize-1);           
            if Registration
                if i == 1
                    [mvmtcorrFramesG,rect] = Crop(Green_indResponse);
                    mvmtcorrFramesF = Crop(Fluo_indResponse); %identify region to use for motion correction computations
                else
                    mvmtcorrFramesG = Crop(Green_indResponse,rect);
                    mvmtcorrFramesF = Crop(Fluo_indResponse,rect);
                end
                refG = mvmtcorrFramesG(:,:,1); %reference frame for motion correction
                refF = mvmtcorrFramesF(:,:,1); %reference frame for motion correction
                usfac = 1; %precision level for dftregistration function
                for j = 1:size(mvmtcorrFramesG,3)
                    output = dftregistration(fft2(refG),fft2(mvmtcorrFramesG(:,:,j)),usfac);
                    tx = output(4); %Keep movement variables in array for possible future use in GLM
                    ty = output(3);
                    Green_indResponse(:,:,j) = imtranslate(Green_indResponse(:,:,j),[tx,ty],'OutputView','same','method','bicubic');
                end
                for j = 1:size(mvmtcorrFramesF,3)
                    output = dftregistration(fft2(refF),fft2(mvmtcorrFramesF(:,:,j)),usfac);
                    tx = output(4); %Keep movement variables in array for possible future use in GLM
                    ty = output(3);
                    Fluo_indResponse(:,:,j) = imtranslate(Fluo_indResponse(:,:,j),[tx,ty],'OutputView','same','method','bicubic');
                end
            end
            if spatialFilt
              for j = 1:size(Green_indResponse,3)
                      Green_indResponse(:,:,j) = imgaussfilt(Green_indResponse(:,:,j),spatialFiltStd);
              end
              for j = 1:size(Fluo_indResponse,3)
                      Fluo_indResponse(:,:,j) = imgaussfilt(Fluo_indResponse(:,:,j),spatialFiltStd);
              end
            end
            if LowpassFilter
                f = designfilt('lowpassfir','PassbandFrequency', PassbandFreq, 'StopbandFrequency', StopbandFreq, ...
                    'StopbandAttenuation', 60,'SampleRate', freq);
                Green_indResponse = reshape(permute(Green_indResponse,[3 1 2]),windowSize,[]);
                Green_indResponse = filtfilt(f,Green_indResponse);
                Green_indResponse = reshape(permute(Green_indResponse,[2 1]),nrows,ncols,windowSize);
                Fluo_indResponse = reshape(permute(Fluo_indResponse,[3 1 2]),windowSize,[]);
                Fluo_indResponse = filtfilt(f,Fluo_indResponse);
                Fluo_indResponse = reshape(permute(Fluo_indResponse,[2 1]),nrows,ncols,windowSize);
            end
            if timeMedFilt
                Green_indResponse = timeMedianFilter(Green_indResponse,timeFiltRange,freq);
                Fluo_indResponse = timeMedianFilter(Fluo_indResponse,timeFiltRange,freq);
            end
            if SGFilt
                Green_indResponse = timeSGFilter(Green_indResponse,SGFiltOrder,timeFiltRange,freq);
                Fluo_indResponse = timeSGFilter(Fluo_indResponse,SGFiltOrder,timeFiltRange,freq);
            end
            %Calculate relative response
            preStim_Frames = Green_indResponse(:,:,1:floor(preStim*freq));
            baseline = mean(preStim_Frames,3);
            Rel_Green_Response_corr = Green_indResponse./baseline;
            
            preStim_Frames = Fluo_indResponse(:,:,1:floor(preStim*freq));
            baseline = mean(preStim_Frames,3);
            Rel_Fluo_Response = Rel_Fluo_Response+(Fluo_indResponse./baseline)./Rel_Green_Response_corr;            
            if singleTrials
                Fluo_Trials{i} = (Fluo_indResponse./baseline)./Rel_Green_Response_corr;
            end
        end
        %Full average relative response
        Rel_Fluo_Response = Rel_Fluo_Response/NbStim;
        responseList = [responseList 'Rel_Fluo_Response'];
        clear Dat_Gptr Green_Frames Green_indResponse greenMap Dat_Fptr Fluo_Frames Fluo_indResponse fluoMap... 
                Rel_Green_Response_corr
    end
    
%        %% Fluorescence channel detected %%
%     %Fluorescence images corrected for Hb absorption with green reflectance
%     %images
%     if IsThereFluo
%         Dat_Fptr = matfile([FolderName filesep 'Data_yellow.mat']); %Blue light replaces amber light
%         Dat_Gptr = matfile([FolderName filesep 'Data_green.mat']);
%         nrows = Dat_Fptr.datSize(1,1);
%         ncols = Dat_Fptr.datSize(1,2);
%         nframes = Dat_Fptr.datLength;
%         freq = Dat_Fptr.Freq;
%         StimFluo = Dat_Fptr.Stim;
%         StimGreen = Dat_Gptr.Stim;
% 
% %         %Background substraction
% %         fluoMap.Data.img = fluoMap.Data.img - Background;
% 
%         %Stimulation on and off times
%         onsetsFluo = find(diff(StimFluo)==1)+1;
%         onsetsFluo = onsetsFluo(GoodEvents);
%         offsetsFluo = find(diff(StimFluo)==-1);
%         offsetsFluo = offsetsFluo(GoodEvents);
% 
%         onsetsGreen = find(diff(StimGreen)==1)+1;
%         onsetsGreen = onsetsFluo(GoodEvents);
%         offsetsGreen = find(diff(StimGreen)==-1);
%         offsetsGreen = offsetsFluo(GoodEvents);
% 
%         %Preallocation for frames of response window (including pre and post stim) & baseline frames 
%         Fluo_Frames = zeros(nrows,ncols,windowSize);
%         Rel_Fluo_Response = zeros(size(Fluo_Frames));
%         preStim_Frames = zeros(nrows,ncols,length(1:floor(preStim*freq)));
%         Fbaseline = zeros(nrows,ncols);
%         Gbaseline = zeros(nrows,ncols);
% 
%         for i = 1:NbStim  
%             Fluo_Frames = fluoMap.Data.img(:,:,(onsetsFluo(i)-floor(preStim*freq)):(onsetsFluo(i)-floor(preStim*freq))+windowSize-1);
%             Green_Frames = greenMap.Data.img(:,:,(onsetsGreen(i)-floor(preStim*freq)):(onsetsGreen(i)-floor(preStim*freq))+windowSize-1);
% 
%             if Registration
%                 mvmtcorrFrames = Crop(Fluo_Frames,rect); %identify region to use for motion correction computations
%                 ref = mvmtcorrFrames(:,:,1); %reference frame for motion correction
%                 usfac = 100; %precision level for dftregistration function
%                 for j = 1:size(mvmtcorrFrames,3)
%                     output = dftregistration(fft2(ref),fft2(mvmtcorrFrames(:,:,j)),usfac);
%                     tx = output(4);
%                     ty = output(3);
%                     Fluo_Frames(:,:,j) = imtranslate(Fluo_Frames(:,:,j),[tx,ty]);
%                 end
%                 mvmtcorrFrames = Crop(Green_Frames,rect); %identify region to use for motion correction computations
%                 ref = mvmtcorrFrames(:,:,1); %reference frame for motion correction
%                 usfac = 100; %precision level for dftregistration function
%                 for j = 1:size(mvmtcorrFrames,3)
%                     output = dftregistration(fft2(ref),fft2(mvmtcorrFrames(:,:,j)),usfac);
%                     tx = output(4);
%                     ty = output(3);
%                     Green_Frames(:,:,j) = imtranslate(Green_Frames(:,:,j),[tx,ty]);
%                 end
%             end
% 
%             if spatialFilt
%               for j = 1:size(Fluo_Frames,3)
%                       Fluo_Frames(:,:,j) = imgaussfilt(Fluo_Frames(:,:,j),spatialFiltStd);
%               end 
%               for j = 1:size(Green_Frames,3)
%                       Green_Frames(:,:,j) = imgaussfilt(Green_Frames(:,:,j),spatialFiltStd);
%               end
%             end
% 
%             if LowpassFilter
%                 Fluo_Frames = reshape(permute(Fluo_Frames,[3 1 2]),windowSize,[]);
%                 Fluo_Frames = filtfilt(f,Fluo_Frames);
%                 Fluo_Frames = reshape(permute(Fluo_Frames,[2 1]),nrows,ncols,windowSize);
%                 Green_Frames = reshape(permute(Green_Frames,[3 1 2]),windowSize,[]);
%                 Green_Frames = filtfilt(f,Green_Frames);
%                 Green_Frames = reshape(permute(Green_Frames,[2 1]),nrows,ncols,windowSize);
%             end
% 
%             if timeMedFilt
%                 Fluo_Frames = timeMedianFilter(Fluo_Frames,timeFiltRange,freq);
%                 Green_Frames = timeMedianFilter(Green_Frames,timeFiltRange,freq);
%             end
% 
%             if SGFilt
%                 Fluo_Frames = timeSGFilter(Fluo_Frames,SGFiltOrder,timeFiltRange,freq);
%                 Green_Frames = timeSGFilter(Green_Frames,SGFiltOrder,timeFiltRange,freq);
%             end
% 
%             %Calculate relative response
%             preStim_Frames = Fluo_Frames(:,:,1:floor(preStim*freq));
%             Fbaseline = mean(preStim_Frames,3); 
%             preStim_Frames = Green_Frames(:,:,1:floor(preStim*freq));
%             Gbaseline = mean(preStim_Frames,3);
%             Rel_Fluo_Response = Rel_Fluo_Response+(Fluo_Frames./Fbaseline)./(Green_Frames./Gbaseline);
%         end   
%         %Full average relative response
%         Rel_Fluo_Response = Rel_Fluo_Response/NbStim;    
%         responseList = [responseList 'Rel_Fluo_Response'];  
%         clear Dat_Fptr Dat_Gptr Fbaseline Gbaseline onsetsFluo onsetsGreen offsetsFluo offsetsGreen Fluo_Frames Green_Frames      
%     end
    %% Red channel detected %%
    if IsThereRed
        Dat_Rptr = matfile([FolderName filesep 'Data_red.mat']);
        freq = Dat_Rptr.Freq;
        StimG = Dat_Rptr.Stim;
        nrows = Dat_Rptr.datSize(1,1);
        ncols = Dat_Rptr.datSize(1,2);
        Red_Frames = redMap.Data.img;
        if substractBG
            Red_Frames = Red_Frames - Background;
        end
        %Indexes of stim onsets
        onsets = find(diff(StimG)==1)+1;
        onsets = onsets(GoodEvents);
        if Detrending
            Red_Frames = detrendImageSeries(Red_Frames,2,1);
        end
        
        %Preallocation for frames of response window (pre-stim+stim+post-stim) & baseline frames
        windowSize = floor(freq*(preStim+StimLength+postStim))+1;
        Red_indResponse = zeros(nrows,ncols,windowSize);
        Rel_Red_Response = zeros(size(Red_indResponse));
        preStim_Frames = zeros(nrows,ncols,length(1:floor(preStim*freq)));
        baseline = zeros(nrows,ncols);
        if singleTrials
            Red_Trials = cell(1,NbStim);
        end
        for i = 1:NbStim
            Red_indResponse = Red_Frames(:,:,(onsets(i)-floor(preStim*freq)):(onsets(i)-floor(preStim*freq))+windowSize-1);
            if Registration
                if i == 1
                    [mvmtcorrFrames,rect] = Crop(Red_indResponse); %identify region to use for motion correction computations
                else
                    mvmtcorrFrames = Crop(Red_indResponse,rect);
                end
                ref = mvmtcorrFrames(:,:,1); %reference frame for motion correction
                usfac = 100; %precision level for dftregistration function
                for j = 1:size(mvmtcorrFrames,3)
                    output = dftregistration(fft2(ref),fft2(mvmtcorrFrames(:,:,j)),usfac);
                    tx = output(4); %Keep movement variables in array for possible future use in GLM
                    ty = output(3);
                    Red_indResponse(:,:,j) = imtranslate(Red_indResponse(:,:,j),[tx,ty],'OutputView','same','method','bicubic');
                end
            end
            if spatialFilt
              for j = 1:size(Red_indResponse,3)
                      Red_indResponse(:,:,j) = imgaussfilt(Red_indResponse(:,:,j),spatialFiltStd);
              end 
            end
            if LowpassFilter
                f = designfilt('lowpassfir','PassbandFrequency', PassbandFreq, 'StopbandFrequency', StopbandFreq, ...
                    'StopbandAttenuation', 60,'SampleRate', freq);
                Red_indResponse = reshape(permute(Red_indResponse,[3 1 2]),windowSize,[]);
                Red_indResponse = filtfilt(f,Red_indResponse);
                Red_indResponse = reshape(permute(Red_indResponse,[2 1]),nrows,ncols,windowSize);
            end
            if timeMedFilt
                Red_indResponse = timeMedianFilter(Red_indResponse,timeFiltRange,freq);
            end
            if SGFilt
                Red_indResponse = timeSGFilter(Red_indResponse,SGFiltOrder,timeFiltRange,freq);
            end
            %Calculate relative response
            preStim_Frames = Red_indResponse(:,:,1:floor(preStim*freq));
            baseline = mean(preStim_Frames,3);
            Rel_Red_Response = Rel_Red_Response + Red_indResponse./baseline;
            if singleTrials
                Red_Trials{i} = Red_indResponse./baseline;
            end
        end
        %Full average relative response
        Rel_Red_Response=Rel_Red_Response/NbStim;
        responseList = [responseList 'Rel_Red_Response'];    
        clear Dat_Gptr Red_Frames Red_indResponse redMap
    end
    %% 
    if IsThereGreen && IsThereRed 
        if singleTrials
            HbO_Trials = cell(1,NbStim);
            HbR_Trials = cell(1,NbStim);
            for i = 1:NbStim
                [HbO_Trials{i}, HbR_Trials{i}] = Hb_Maps(Green_Trials{i},Red_Trials{i},Fluo);
            end
        else                
            [HbO_Response, HbR_Response] = Hb_Maps(Rel_Green_Response,Rel_Red_Response,Fluo);    
            responseList = [responseList 'HbO_Response' 'HbR_Response'];
        end
    end
    if singleTrials
        %Save single trial images time series
        responseList = {'HbO_Trials' 'HbR_Trials' 'SC_Trials' 'Fluo_Trials' 'Timeline' 'freq' 'preStim' 'postStim' 'StimLength' 'Ref_Frame'};
        savingFile = [FolderName filesep 'Single_Trial_Frames'];
        save(savingFile,responseList{1});
        for responseIndex = 2:length(responseList)
            save(savingFile,responseList{responseIndex},'-append');
        end
    else
        %Save mean images time series
        responseList = [responseList 'Timeline' 'freq' 'preStim' 'postStim' 'StimLength' 'Ref_Frame'];
        savingFile = [FolderName filesep 'Mean_Frames'];
        save(savingFile,responseList{1});
        for responseIndex = 2:length(responseList)
            save(savingFile,responseList{responseIndex},'-append');
        end
    end
close(d);
close(f);
end    