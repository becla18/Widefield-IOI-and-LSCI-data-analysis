function out = OpenIOI_NewSyst(FolderName, Binning, Version, DEF_STIMSLAVE)

%%%%DEFINES -> THESE CONSTANTS ARE HARDCODED!!!! DO NOT CHANGE THEM.
NOFPF = 256; %Max number of frame per .bin file
DEF_VISUEL = 0;

m_fFeedBack = figure('Position',[750 100 600 600],...             
    'Name', 'Output Display', 'NumberTitle', 'off');
OStream = uicontrol('Parent', m_fFeedBack, 'Style', 'edit',...
           'String', '', 'Max', 100, 'Min', 1,...
           'HorizontalAlignment', 'left',...
           'Enable', 'off', 'Position', [50 50 500 500]);
       
       
%OStream = m_OutStream;

if( exist('Config.m','file') )
    Config;
else
    DEF_VISUEL = 0;
    DEF_FLUO = 0;
    DEF_SPECKLE = 1;
end
%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% Acq. Info file:
%%%%%%%%%%%%%%%%%%%%%
AcqInfoStream = readtable([FolderName filesep 'info.txt'],...
    'Delimiter',':','ReadVariableNames',false, 'ReadRowNames',true);

if( isempty(OStream) )
    disp('Recovering stimulation parameters')
    disp('**************************');
else
    OStream.String = sprintf('%s\r%s',...
        'Recovering stimulation parameters',...
        OStream.String);
    drawnow;
end
if( DEF_VISUEL )
    tStim = AcqInfoStream{'Stimulation',1};
    tAIChan = AcqInfoStream{'AINChannels',1};
else
    if( Version > 1)
        if( sum(ismember(AcqInfoStream.Properties.RowNames,'Stimulation1')) )
            tStim = AcqInfoStream{'Stimulation1',1};
        elseif ( sum(ismember(AcqInfoStream.Properties.RowNames,'Simulation') ) )
            tStim = AcqInfoStream{'Simulation',1};
            Version = -1;
        end
    else
        if( sum(ismember(AcqInfoStream.Properties.RowNames,'Stimulation')) )
            
            tStim = AcqInfoStream{'Stimulation',1};
        elseif ( sum(ismember(AcqInfoStream.Properties.RowNames,'Simulation') ) )
            tStim = AcqInfoStream{'Simulation',1};
            Version = -1;
        end
    end
    tAIChan = AcqInfoStream{'AINChannels',1};
end

if( DEF_STIMSLAVE )
    tStim = 1;
end

if( iscell(tAIChan) )   
    tAIChan =  str2double(cell2mat(tAIChan));
end

if ( iscell(tStim) )
    tStim = str2double(cell2mat(tStim));
end

%%%%%%%%%%%%%%%%%%%%%
% Stimulation detected
%%%%%%%%%%%%%%%%%%%%%
if( tStim )
    IOIReadStimFile_NS(FolderName, Version, tAIChan, DEF_STIMSLAVE);
else
    if( isempty(OStream) )
        fprintf('No stimulation detected. \n');
    else
        OStream.String = sprintf('%s\r%s',...
            'No stimulation detected.',...
            OStream.String);
        drawnow;
    end
    
    %TODO: resting state.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Images sequence validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgFilesList = dir([FolderName filesep 'img_0*.bin']);
aiFilesList = dir([FolderName filesep 'ai_*.bin']);

%Version Management:
if( Version == 3)
    hWima = 5;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint64', 3, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2 + 3*8;
elseif( Version == 2)
    hWima = 5;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint64', 1, 'framej';'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2 + 8;
elseif( Version == 1 )
    %TODO: To be validated
    hWima = 4;
    hWai = 5;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', hWima, 'header'; 'uint64', 1, 'frame'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2 + 8;
elseif( Version == -1 )
    %TODO: To be validated
    hWima = 0;
    hWai = 4;
    header = memmapfile([FolderName filesep imgFilesList(1).name], ...
        'Offset', 0, 'Format', {'int32', 4, 'header'}, 'repeat', 1);
    
    nx=header.Data.header(2);
    ny=header.Data.header(3);
    frameFormat = {'uint16', [double(nx), double(ny)], 'imgj'};
    ImRes_XY = [nx, ny];
    SizeImage = nx*ny*2;
else
    disp(['Error! System version ' int2str(Version) ' is not suported by this software']);
end

%Find number of images
NombreImage = 0;
for ind = 1:size(imgFilesList,1)
    data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
    NombreImage = NombreImage+size(data.Data,1);
end

%Create image id vector idImg
if( Version == 3 )
   idImg = zeros(NombreImage, 3);
else
   idImg = zeros(NombreImage, 1);
end

%Fill idImg with image indexes
if(Version>0)
    for ind = 1:size(imgFilesList,1)
        data = memmapfile([FolderName filesep imgFilesList(ind).name],'Offset',hWima*4,'Format',frameFormat,'repeat',inf);
        idImg((NOFPF*(ind-1)+1):(NOFPF*(ind-1)+size(data.Data,1)),:) = cell2mat(arrayfun(@(x) data.Data(x).framej, 1:size(data.Data,1),'UniformOutput',false))';
    end
else
    idImg(1:end)=linspace(1,size(idImg,1),size(idImg,1));
end

clear nx ny data header ind;

% Verbose
if( isempty(OStream) )
    disp(['Opening of: ' FolderName]);
    disp(['Number of Frames acquired: ' int2str(NombreImage)]);
    disp(['Frames'' resolution: ' int2str(ImRes_XY(1)) ' pix X ' int2str(ImRes_XY(2)) ' pix']);
else
    OStream.String = sprintf('%s%s\r%s%s\r%s%s%s%s%s\r%s',...
        'Opening of: ', FolderName,...
        'Number of Frames acquired: ', int2str(NombreImage),...
        'Frames'' resolution: ', int2str(ImRes_XY(1)),...
        ' pix X ', int2str(ImRes_XY(2)), ' pix',...
        OStream.String);
    drawnow;
end
% end of Verbose

%%%%%%%%%%%%%%%%%%%%%
% Binning
%%%%%%%%%%%%%%%%%%%%%
if( Binning )
    %Verbose
    if( isempty(OStream) )
        disp('Binning option is ON');
    else
        OStream.String = sprintf('%s\r%s',...
            'Binning option is ON',...
            OStream.String);
        drawnow;
    end
    %end of Verbose
    
    Rx = round(ImRes_XY(1)/2);
    Ry = round(ImRes_XY(2)/2);
   
else
    Rx = ImRes_XY(1);
    Ry = ImRes_XY(2);
end


%%%%%%%%%%%%%%%%%%%%%
% Analog inputs
%%%%%%%%%%%%%%%%%%%%%
AnalogIN = [];
if(Version > 0)
    for ind = 1:size(aiFilesList,1)
        data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', hWai*4, 'Format', 'double', 'repeat', inf);
        tmp = data.Data;
        tmp = reshape(tmp, 1e4, tAIChan, []);
        tmp = permute(tmp,[1 3 2]);
        tmp = reshape(tmp,[],tAIChan);
        AnalogIN = [AnalogIN; tmp];
    end
else
    % Prototype version of system
    for ind = 1:size(aiFilesList,1)
        data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', hWai*4, 'Format', 'double', 'repeat', inf);
        tmp = data.Data;
        tmp = reshape(tmp, 1e4, 8, []);
        tmp = permute(tmp,[1 3 2]);
        tmp = reshape(tmp,[],8);
        flip_tmp = zeros(size(tmp));
        flip_tmp(:,1)=tmp(:,4);
        flip_tmp(:,2)=tmp(:,1);
        flip_tmp(:,3)=tmp(:,2);
        flip_tmp(:,4)=tmp(:,3);
        flip_tmp(:,5:8)=tmp(:,5:8);
        AnalogIN = [AnalogIN; flip_tmp];
    end
end
clear tmp ind data;

%%%%
%Stimulation Params
%%%%
Str = [];
if( exist([FolderName filesep 'StimParameters.mat'], 'file') )
    load([FolderName filesep 'StimParameters.mat']);
    if( NbStim > 0 )
        Str = sprintf('%s\r%s%s\r%s%s%s',...
            'Stim detected: yes',...
            'Number of events: ', int2str(NbStim),...
            'Length of each event: ', int2str(StimLength), ' sec');
        bStim = 1;
    else
        Str = sprintf('%s',...
            'Stim detected: no');
        bStim = 0;
    end
else
    bStim = 0;
end
if( isempty(OStream) )
    fprintf(Str);
    fprintf('\n');
else
    OStream.String = sprintf('%s\r%s',...
        Str,...
        OStream.String);
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%
% Camera Trigs
%%%%%%%%%%%%%%%%%%%%%%%
CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5))+1;
StartDelay = round(CamTrig(1)/10);
EndDelay = round((length(AnalogIN(:,1)) - CamTrig(end))/10);

% Verbose
Str = sprintf('%s\r%s\r%s',...
    ['Camera Trigs detected: ' int2str(length(CamTrig))],...
    ['Recording of analog inputs starts ' int2str(StartDelay) ' ms before the first trigger.'],... 
    ['Recording of analog inputs ends ' int2str(EndDelay) ' ms after the last trigger.']); 
if( isempty(OStream) )
    fprintf(Str);
    fprintf('\n');
else
    OStream.String = sprintf('%s\r%s',...
        Str,...
        OStream.String);
    drawnow;
end
% end of Verbose
clear StartDelay EndDelay
%Missing frames?
%Not taken care of for now...

%Less trig than images... something's wrong!
if( length(CamTrig) < NombreImage  ) 
    disp('IOI Error: Analog recordings and Image files don''t match. Impossible to continue further.');
    out = 'Error';
    return
end

%%%%%%%%%%%%%%%%%%%%%%%
% Illumination Sequence
%%%%%%%%%%%%%%%%%%%%%%%
% Red           = 0001;
% Yellow/Amber  = 0010;
% Green         = 0100;
% Other         = 1000;

tColor = AcqInfoStream{'Illumination',1};
if(Version > 0)
    if( iscell(tColor) )
        tColor = str2double(cell2mat(tColor));
    end
    if( DEF_VISUEL )
        bFluo = (tColor > 3); tColor = mod(tColor,4);
        bYellow = (tColor > 1); tColor = mod(tColor,2);
        bRed = (tColor > 0);
        bGreen = 0;
        clear fInfo tColor;
    else
        bFluo = (tColor > 7); tColor = mod(tColor,8);
        bGreen = (tColor > 3); tColor = mod(tColor,4);
        bYellow = (tColor > 1); tColor = mod(tColor,2);
        bRed = (tColor > 0);
        clear fInfo tColor;
    end
else
    % Prototype system, only red and green
    bGreen = 1;
    bRed = 1;
    bFluo = 0;
    bYellow = 0;
end


Freq = AcqInfoStream{'FrameRateHz',1};
if( iscell(Freq) )
    Freq = str2double(cell2mat(Freq));
end

nbColors = (bFluo + bGreen + bYellow + bRed);
if( bFluo )
    if( DEF_FLUO )
       
        if( exist([FolderName filesep 'Data_Fluo.mat'],'file') )
            delete([FolderName filesep 'Data_Fluo.mat']);
        end
        fSpeckle = matfile([FolderName filesep 'Data_Fluo.mat'],'Writable',true);
        fSpeckle.datFile = [FolderName filesep 'fChan.dat'];
        fSpeckle.datSize = [2*Rx, 2*Ry];
        fSpeckle.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
        fSpeckle.Freq = Freq/nbColors;
        cSpeckle = 1;
        fidS = fopen([FolderName filesep 'fChan.dat'],'w');
    elseif( DEF_SPECKLE )
        
        if( exist([FolderName filesep 'Data_speckle.mat'],'file') )
            delete([FolderName filesep 'Data_speckle.mat']);
        end
        fSpeckle = matfile([FolderName filesep 'Data_speckle.mat'],'Writable',true);%Create Data_speckle.mat file, which contains sChan.dat
        fSpeckle.datFile = [FolderName filesep 'sChan.dat'];
        if Binning
            fSpeckle.datSize = [2*Rx, 2*Ry];
        else
            fSpeckle.datSize = [Rx, Ry];
        end
        fSpeckle.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
        fSpeckle.Freq = Freq/nbColors;
        cSpeckle = 1;
        fidS = fopen([FolderName filesep 'sChan.dat'],'w');
    end
end
if( bRed ) %Create Data_red.mat file, which contains rChan.dat
   
    if( exist([FolderName filesep 'Data_red.mat'],'file') )
        delete([FolderName filesep 'Data_red.mat']);
    end
    fRed = matfile([FolderName filesep 'Data_red.mat'],'Writable',true);
    fRed.datFile = [FolderName filesep 'rChan.dat'];
    fRed.datSize = [Rx, Ry];
    fRed.Stim = zeros(floor(NombreImage/nbColors), 1, 'single'); %Stim vector during red frames
    fRed.Freq = Freq/nbColors;
    cRed = 1;
    fidR = fopen([FolderName filesep 'rChan.dat'],'w');
end
if( bYellow )
    
    if( exist([FolderName filesep 'Data_yellow.mat'],'file') )
        delete([FolderName filesep 'Data_yellow.mat']);
    end
    fYellow = matfile([FolderName filesep 'Data_yellow.mat'],'Writable',true);
    fYellow.datFile = [FolderName filesep 'yChan.dat'];
    fYellow.datSize = [Rx, Ry];
    fYellow.Stim = zeros(floor(NombreImage/nbColors),1, 'single');
    fYellow.Freq = Freq/nbColors;
    cYellow = 1;
    fidY = fopen([FolderName filesep 'yChan.dat'],'w');
end
if( bGreen )
       
    if( exist([FolderName filesep 'Data_green.mat'],'file') )
        delete([FolderName filesep 'Data_green.mat']);
    end
    fGreen = matfile([FolderName filesep 'Data_green.mat'],'Writable',true);
    fGreen.datFile = [FolderName filesep 'gChan.dat'];
    fGreen.refImgFile = [FolderName filesep 'gChan_Ref.dat'];
    fGreen.datSize = [Rx, Ry];
    fGreen.Stim = zeros(floor(NombreImage/nbColors), 1, 'single');
    fGreen.Freq = Freq/nbColors;
    cGreen = 1;
    fidG = fopen([FolderName filesep 'gChan.dat'],'w');
    fidG_Ref = fopen([FolderName filesep 'gChan_Ref.dat'],'w');
end

%Interpolation for bad or missing frames
[~, idxOri] = unique(idImg(:,1)); %returns vector of indexes of idImg without repetitions and in ascending order
%badFrames = find(~ismember(1:NombreImage, uniqueFramesID));
Conseq = conv(idImg(:,1),[1 1 -2],'same') == 3; %If 3 frames are consecutive,conv(idImg(:,1),[1 1 -2],'same') gives 3 at the index of the 1rst frame 
idxE = idImg(find(diff(Conseq,1,1)==1) + 1) -1;
Conseq = conv(idImg(:,1),[0 0 1 1 -2],'same') == 3;
idxS = idImg(find(diff(Conseq,1,1)==-1)) + 1;

badFrames = [];
for ind = 1:length(idxS)
    badFrames = [badFrames, idxS(ind):idxE(ind)];
end
idxOri(ismember(idImg(idxOri,1),badFrames)) = [];

%%% Lookup Table For missing frames
InterpLUT = zeros(8,1);
if( ~isempty(badFrames) )
    InterpLUT = zeros(8,size(badFrames,2));
    % 1: Frame before
    % 2: File Number where to find this frame
    % 3: Image ID in this file
    % 4: Frame After
    % 5: File Number where to find this frame
    % 6: Image ID in this file
    % 7: Ratio between these frame and the one missing
    % 8: Frame tag id
    for ind = 1:size(badFrames,2)
        tmpID = badFrames(ind);
        tmpBefore = tmpID - (nbColors:nbColors:(tmpID-1));
        idx = find(ismember(tmpBefore,idImg(:,1))&~ismember(tmpBefore,badFrames),1,'first');
        tmpBefore = tmpBefore(idx);
        InterpLUT(1,ind) = tmpBefore;
        idx = find(tmpBefore == idImg);
        InterpLUT(2,ind) = floor((idx-1)/256) + 1;
        InterpLUT(3,ind) = rem((idx-1),256) + 1;
        
        tmpAfter = tmpID + (nbColors:nbColors:(NombreImage));
        idx = find(ismember(tmpAfter,idImg(:,1))&~ismember(tmpAfter,badFrames),1,'first');
        tmpAfter = tmpAfter(idx);
        InterpLUT(4,ind) = tmpAfter;
        idx = find(tmpAfter == idImg);
        InterpLUT(5,ind) = floor((idx-1)/256) + 1;
        InterpLUT(6,ind) = rem((idx-1),256) + 1;
        
        tmpRatio =  (tmpID - InterpLUT(1,ind))./...
            (InterpLUT(4,ind) - InterpLUT(1,ind));
        InterpLUT(7,ind) = tmpRatio;
        InterpLUT(8,ind) = badFrames(ind);
    end
    clear tmpRatio tmpAfter tmpBefore;
    %%% Interpolation of missing frames
    TmpFrames = struct('framej',[], 'imgj',[]);
    for ind = 1:size(InterpLUT,2)
        dBefore = memmapfile([FolderName filesep...
            imgFilesList(InterpLUT(2,ind)).name],...
            'Offset', hWima*4 + (InterpLUT(3,ind) - 1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        dAfter = memmapfile([FolderName filesep...
            imgFilesList(InterpLUT(5,ind)).name],...
            'Offset', hWima*4 + (InterpLUT(6,ind) - 1)*SizeImage,...
            'Format', frameFormat, 'repeat', 1);
        TmpFrames(ind).imgj = uint16(round(InterpLUT(7,ind)*double(dAfter.Data.imgj - dBefore.Data.imgj))) + dBefore.Data.imgj;
        TmpFrames(ind).framej = uint64([InterpLUT(8,ind), 1, 1]);
    end
    
    fid = fopen([FolderName filesep 'img_interp.bin'],'w');
    for ind = 1:size(InterpLUT,2)
        fwrite(fid, TmpFrames(ind).framej, 'uint64');
        fwrite(fid, TmpFrames(ind).imgj, 'uint16');
    end
    fclose(fid);
end

%Rebuilding addresses for each frames...
ImAddressBook = zeros(NombreImage,2);
for ind = 1:NombreImage
    if( ismember(ind, InterpLUT(8,:)) )
        fidx = find( ind == InterpLUT(8,:), 1, 'first');
        ImAddressBook(ind,1) = size(imgFilesList,1) + 1;
        ImAddressBook(ind,2) = fidx;
    elseif( ismember(ind, idImg(idxOri)) )
        fidx = find( ind == idImg, 1, 'first');
        ImAddressBook(ind,1) = floor((fidx-1)/256) + 1;
        ImAddressBook(ind,2) = rem(fidx-1, 256) + 1;
    end
end %ImAddressBook contains in column 1 no. of .bin file and in column 2 the number of the associated good frame

%Saving infos...
if( ~strcmp(FolderName(end), filesep) )
    FolderName = [FolderName filesep];
end
save([FolderName 'ImagesLUT.mat'], 'ImAddressBook');

%%%%
% Images Classification and filtering
%%%%
ind = 1; 
if( bFluo ) %Data_speckle filling
    %%%
    ind = nbColors; %Hypothesis: Laser is last in illumination sequence
    %%%
    OStream.String = sprintf('%s\r%s',...
        'Auxiliary channel classification:',...
        OStream.String);
    drawnow;
    
    tags = ind:nbColors:NombreImage;
    Images = zeros(2*Rx, 2*Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    StaticStr = OStream.String;
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if ImAddressBook(indF,1)>0
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        end
        
        %No Binning for speckle?%%
%         if( Binning )
%             img = interp2(double(dat.Data.imgj), (1:2:size(dat.Data.imgj,2))',...
%                 (1:2:size(dat.Data.imgj,1)));
%         else
            img = dat.Data.imgj;
%         end

        Images = img;
        fwrite(fidS, Images, 'single');
        
        if( bStim )
            fSpeckle.Stim(cSpeckle,1) = Stim(indF);
        else
            fSpeckle.Stim(cSpeckle,1) = 0;
        end
        cSpeckle = cSpeckle + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
            if( isempty(OStream) )
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            else
                OStream.String = sprintf('%s\r%s',...
                    ['Completion: ' int2str(P) '%'],...
                    StaticStr);
                drawnow;
            end
            indT = indT + 1;
        end
    end
    ind = ind + 1;
    fSpeckle.datLength = cSpeckle - 1;
    fclose(fidS);
    OStream.String = sprintf('%s\r%s',...
        'Done.',...
        StaticStr);
    drawnow;
end
if( bRed )
    %%%
    ind = 1; %Hypothesis: Red is first in illumination sequence
    %%%
   OStream.String = sprintf('%s\r%s',...
        'Red channel classification:',...
        OStream.String);
    drawnow;
    
    if( Version < 0 )
        tags = 2:nbColors:NombreImage;
    else
        tags = ind:nbColors:NombreImage;
    end
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    StaticStr = OStream.String;
    indT = 1;  
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if ImAddressBook(indF,1)>0
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        end
        
        if( Binning )
           img = interp2(double(dat.Data.imgj), (1:2:size(dat.Data.imgj,2))',...
               (1:2:size(dat.Data.imgj,1)));
        else
           img = dat.Data.imgj;
        end
        Images = img;
        fwrite(fidR, Images, 'single');
        
        if( bStim )
            fRed.Stim(cRed,1) = Stim(indF);
        else
            fRed.Stim(cRed,1) = 0;
        end
        cRed = cRed + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
            if( isempty(OStream) )
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            else
                OStream.String = sprintf('%s\r%s',...
                    ['Completion: ' int2str(P) '%'],...
                    StaticStr);
                drawnow;
            end
            indT = indT + 1;
        end
    end
    
    ind = ind + 1;
    fRed.datLength = cRed - 1;
    fclose(fidR);
     OStream.String = sprintf('%s\r%s',...
        'Done.',...
        StaticStr);
    drawnow;
end
if( bYellow )
    OStream.String = sprintf('%s\r%s',...
        'Yellow channel classification:',...
        OStream.String);
    drawnow;
    
    tags = ind:nbColors:NombreImage;
    Images = zeros(Rx, Ry, 'single');
    
    PrcTag = round(linspace(0, length(tags), 20));
    StaticStr = OStream.String;
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
        
        if ImAddressBook(indF,1)>0
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        end
               
        if( Binning )
           img = interp2(double(dat.Data.imgj), (1:2:size(dat.Data.imgj,2))',...
               (1:2:size(dat.Data.imgj,1)));
        else
           img = dat.Data.imgj;
        end
        Images = img;
        fwrite(fidY, Images, 'single');
        
        if( bStim )
            fYellow.Stim(cYellow,1) = Stim(indF);
        else
            fYellow.Stim(cYellow,1) = 0;
        end
        cYellow = cYellow + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
            if( isempty(OStream) )
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            else
                OStream.String = sprintf('%s\r%s',...
                    ['Completion: ' int2str(P) '%'],...
                    StaticStr);
                drawnow;
            end
            indT = indT + 1;
        end
    end
   
    ind = ind + 1;
    fYellow.datLength = cYellow - 1;
    fclose(fidY);
     OStream.String = sprintf('%s\r%s',...
        'Done.',...
        StaticStr);
    drawnow;
end
if( bGreen )

    OStream.String = sprintf('%s\r%s',...
        'Green channel classification:',...
        OStream.String);
    drawnow;
    
    if( Version < 0 )
       tags = 1:nbColors:NombreImage;
    else
       tags = ind:nbColors:NombreImage;
    end
    Images = zeros(Rx, Ry, 'single');
    if Binning                                 %Create reference image for speckle ROI selection (19/06/19)
        Ref_imgage = zeros(Rx,Ry,'single');
    else
        Ref_image = Images;
    end
    
    PrcTag = round(linspace(0, length(tags), 20));
    StaticStr = OStream.String;
    indT = 1;
    for indI = 1:length(tags)
        indF = tags(indI);
   
        if ImAddressBook(indF,1)>0
        if( ImAddressBook(indF,1) <= size(imgFilesList,1) )
            dat =   memmapfile([FolderName filesep...
                imgFilesList(ImAddressBook(indF,1)).name],...
                'Offset', hWima*4 + (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        else
            dat =   memmapfile([FolderName filesep 'img_interp.bin'],...
                'Offset', (ImAddressBook(indF,2)-1)*SizeImage,...
                'Format', frameFormat, 'repeat', 1);
        end
        end
       
        if( Binning )
           img = interp2(double(dat.Data.imgj), (1:2:size(dat.Data.imgj,2))',...
               (1:2:size(dat.Data.imgj,1)));
        else
           img = dat.Data.imgj;
        end
        
        if indI == 1
            Ref_image = dat.Data.imgj;
            fwrite(fidG_Ref,Ref_image,'single');
            fclose(fidG_Ref);
        end
        
        Images = img;
        fwrite(fidG, Images, 'single');     
        
        if( bStim )
            fGreen.Stim(cGreen,1) = Stim(indF);
        else
            fGreen.Stim(cGreen,1) = 0;
        end
        cGreen = cGreen + 1;
        
        if( indI >= PrcTag(indT) )
            P = round((100*PrcTag(indT))/length(tags));
            if( isempty(OStream) )
                fprintf('%d%% .. ', P);
                if( indT == 10 )
                    fprintf('\n');
                end
                
            else
                OStream.String = sprintf('%s\r%s',...
                    ['Completion: ' int2str(P) '%'],...
                    StaticStr);
                drawnow;
            end
            indT = indT + 1;
        end
    end

    ind = ind + 1;
    fGreen.datLength = cGreen - 1;
    fclose(fidG);
     OStream.String = sprintf('%s\r%s',...
        'Done.',...
        StaticStr);
    drawnow;
end

fprintf('\n');
%Verbose
if( isempty(OStream) )
    fprintf(['Done with file ' FolderName]);
    fprintf('\n');
else
    OStream.String = sprintf('%s\r%s\r%s',...
        ['Done with file ' FolderName],...
        '************* ',...
        OStream.String);
    drawnow;
end
%end of Verbose
close(m_fFeedBack);
end
