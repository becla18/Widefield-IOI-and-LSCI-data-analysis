function IOIReadStimFile_NS(FolderName, Version, NbChan, bSlave)

aiFilesList = dir([FolderName filesep 'ai_*.bin']);

AnalogIN = [];
if Version > 0
    for ind = 1:size(aiFilesList,1)
        data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', 5*4, 'Format', 'double', 'repeat', inf);
        tmp = data.Data;
        tmp = reshape(tmp, 1e4, NbChan, []);
        tmp = permute(tmp,[1 3 2]);
        tmp = reshape(tmp,[],NbChan);
        AnalogIN = [AnalogIN; tmp];
    end
else
    % Prototype version of system
    for ind = 1:size(aiFilesList,1)
        data = memmapfile([FolderName filesep aiFilesList(ind).name], 'Offset', 4*4, 'Format', 'double', 'repeat', inf);
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
clear tmp flip_tmp ind data;
CamTrig = find((AnalogIN(1:(end-1),1) < 2.5) & (AnalogIN(2:end,1) >= 2.5))+1;

if( ~bSlave )
    sChan = 2;
else
    sChan = 3;
end
maxStimVoltage = max(AnalogIN(:,sChan));
StimON = find((AnalogIN(1:(end-1), sChan) < maxStimVoltage/2) & (AnalogIN(2:end, sChan) >= maxStimVoltage/2))+1;
%StimON = find(diff(AnalogIN(1:end,sChan))>4.5)+1; %mod 17/06/20 to decrease sensitivity to noise
%StimLim = find(diff(AnalogIN(1:end,sChan))<-4.5);
StimLim = find((AnalogIN(1:(end-1), sChan) > maxStimVoltage/2) & (AnalogIN(2:end, sChan) <= maxStimVoltage/2));
%Fix to avoid detection of single point voltage bursts
VburstNoise = find(StimON==StimLim);
StimON(VburstNoise)=[];
StimLim(VburstNoise)=[];
if( ~isempty(StimON) )
    Period = median(StimON(2:end)-StimON(1:(end-1)))/10000;
    Freq = 1/Period;
    %Width = sum(AnalogIN(StimON(1):StimON(2),sChan) >2.5)/(Period*10000);
    Width = (length(AnalogIN(StimON(1):StimLim(1),sChan)))/(10000); %mod 17/06/19 to account for smaller voltage changes
    StimStarts = find(diff(StimON)>20000)+1;
    StimStarts = [1;StimStarts]; %Indexes of StimON vector that correspond to event onsets
    NbStim = length(StimStarts);
    StimStops = find(diff(StimLim)>20000); %Indexes of StimLim vector that correspond to event offsets
    StimStops = [StimStops;length(StimLim)];
    if( NbStim == length(StimON) ) %Single Pulse trigged Stims
        %StimLim = find((AnalogIN(1:(end-1), sChan) > 0.8) & (AnalogIN(2:end, sChan) <= 0.8))+1;
        StimLength = mean(StimLim - StimON)./1e4;
        InterStim_min = min((StimON(2:end) - StimLim(1:(end-1)))./10000);
        InterStim_max = max((StimON(2:end) - StimLim(1:(end-1)))./10000);
        Stim = zeros(length(AnalogIN(:,2)),1);
        for indS = 1:NbStim
           Stim(StimON(indS):StimLim(indS)+0.9*1e4) = 1; 
        end
    else %Pulse train Stim
        StimLength = round(length(StimON)/(NbStim*Freq));
        InterStim_min = min((StimON(StimStarts(2:end)) - StimON(StimStarts(1:end-1)))./10000);
        InterStim_max = max((StimON(StimStarts(2:end)) - StimON(StimStarts(1:end-1)))./10000);
        InterStim_min = InterStim_min - StimLength;
        InterStim_max = InterStim_max - StimLength;
    
        Stim = zeros(length(AnalogIN(:,2)),1);
        for indS = 1:length(StimStarts)
            Stim(StimON(StimStarts(indS)):StimLim(StimStops(indS))) = 1;
        end
    end
    
    Stim = Stim(CamTrig);
    save([FolderName filesep 'StimParameters.mat'], 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
else
    disp('No Stimulations detected. Resting State experiment?');
    Stim = 0;
    StimLength = 0;
    NbStim = 0;
    InterStim_min = 0;
    InterStim_max = 0;
    save([FolderName filesep 'StimParameters.mat'], 'Stim', 'StimLength', 'NbStim', 'InterStim_min', 'InterStim_max');
end
end