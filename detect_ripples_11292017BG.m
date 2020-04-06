function detect_ripples_11292017BG
%--------------------------------------------------------------------------
% Detects ripples based on the RMS of the bandpassed LFP when it goes above
% a threshold.
% Ripple edges are detected based on a second threshold. Only segments
% corresponding to NREM and IS sleep states are considered.
% This method follows Zugaro's methodology as defined in the FMA Toolbox
% and Girardeau et al, Nature Neuro., 2009.
%--------------------------------------------------------------------------
%% Select Stage Scored File:
working_dir=pwd;
current_dir='C:\SleepData';
cd(current_dir);
scoredCheck = 0;
while isequal(scoredCheck, 0)
    [scoredFile, scoredPath] = uigetfile({'*.xls','Excel 1997-2003 File (*.xls)'},...
        'Select the Sleep Scored File');
    if isequal(scoredFile,0) || isequal(scoredPath,0)
        uiwait(errordlg('You need to select a file. Please try again',...
            'ERROR','modal'));
    else
        cd(working_dir);
        stageScoredFile= fullfile(scoredPath, scoredFile);
        %Load sleep scored file:
        try
            [numData, stringData] = xlsread(stageScoredFile);
            scoredCheck = 1;
        catch %#ok<*CTCH>
            % If file fails to load, it will notify user and prompt to
            % choose another file.
            uiwait(errordlg('Check if the scored file is saved in Microsoft Excel format.',...
             'ERROR','modal'));
         scoredCheck = 0;
        end

    end
end

%% Detect if states are in number or 2-letter format:
if isequal(size(numData,2),3)
    scoredStates = numData(:,2:3);
    clear numData stringData
else
    scoredStates = numData(:,2);
    clear numData
    stringData = stringData(3:end,3);
    [stateNumber] = stateLetter2NumberConverter(stringData);
    scoredStates = [scoredStates stateNumber];
    clear stateNumber stringData
end
epochInSeconds = scoredStates(2,1) - scoredStates(1,1);
startTime = scoredStates(1,1) * 10^6;
endTime = (scoredStates(end,1) + epochInSeconds) * 10^6;

%% Load Neuralynx continuous file
[CSCFilename, CSCFilePath] = uigetfile({'*.ncs',...
        'Pick CSC files.'},'Select Continuously Sampled Channel File');
cscFile = fullfile(CSCFilePath, CSCFilename);
[Timestamps, Samples] = Nlx2MatCSC(cscFile, [1 0 0 0 1], 0, 4, [startTime endTime]);

% Reshape LFP array:
[m1,n1]=size(Samples);  
newM = m1*n1;
nsamp = m1;
LFP = reshape(Samples, newM, 1);
clear Samples

% Fill in the time stamps:
eelen=length(Timestamps);
precise_timestamps = zeros(eelen*nsamp, 1);
idx = 1;
for i = 1:eelen
  if i < eelen
    t1 = Timestamps(i);
    t2 = Timestamps(i+1);
    interval = (t2-t1)/nsamp;
    trange =([t1 : interval : t2-interval]);
    precise_timestamps(idx:idx+nsamp-1,1) = trange;
  else
    t1 = Timestamps(i);
    t2 = t1+interval*nsamp;
    trange =([t1 :interval : t2-interval]);
    precise_timestamps(idx:idx+nsamp-1,1) = trange;
  end
  idx = idx + nsamp;
end
clear Timestamps eelen t1 t2 interval trange
% Convert from usec to seconds:
timeStamps = precise_timestamps/1000000;
clear precise_timestamps
Fs = 1/(timeStamps(2)-timeStamps(1));

%% Filter Signal for EEG/LFP:
if Fs > 1050
    [z, p, k] = ellip(7,1,60, 300/(Fs/2),'low');
    [sos, g] = zp2sos(z,p,k);
    LFP = filtfilt(sos, g, LFP);
end

%% Downsample Signal:
targetFs = 1000;
DsFactor = floor(Fs/targetFs); % Determine downsampling factor to get to target sampling frequency
  
timeStamps = timeStamps(1:DsFactor:end);
LFP = LFP(1:DsFactor:end);
Fs = 1/(timeStamps(2)-timeStamps(1));  % New downsampled sampling frequency

%% Filter Signal for RIPPLE Frequency:
rippleBand = [100 200]; % Ripple band as defined by Zugaro
[z, p, k] = ellip(7,1,60, rippleBand/(Fs/2),'bandpass');
[sos, g] = zp2sos(z,p,k);
filtRipple = filtfilt(sos, g, LFP);
% clear LFP

%% Assign states to each data point as a new vector:
lengthSignal = length(filtRipple);
sleepsamp = zeros(lengthSignal,1);
for i = 1:size(scoredStates, 1)
    if isequal(i, size(scoredStates, 1))
        subIntervalIndx = find(timeStamps >= scoredStates(i,1) & timeStamps < (scoredStates(i,1) + 10));
    else
        subIntervalIndx = find(timeStamps >= scoredStates(i,1) & timeStamps < scoredStates(i+1,1));
    end
    if ~isempty(subIntervalIndx)
        sleepsamp(subIntervalIndx) = scoredStates(i,2);
        clear subIntervalIndx
    end
end

%% Calculate the root mean square at each point with a 6 msec window size:
squaredSignal = filtRipple.^2;
% clear filtRipple
windowSeconds = 0.006; % 6 millisecond window
windowSize = windowSeconds * Fs;
halfWindow = floor(windowSize/2);
rmsStart = halfWindow+1;
rmsEnd = lengthSignal - halfWindow;
rmsSignal = zeros(lengthSignal,1);
for i = rmsStart:rmsEnd
    rmsSignal(i) = sqrt(sum(squaredSignal(i-halfWindow:i+halfWindow))/(2*halfWindow+1));
end
clear squaredSignal

%% Determine upper and lower thresholds
% meanNremRMS = mean(rmsSignal(sleepsamp==2 | sleepsamp==6));
stdDevRMS = std(rmsSignal(sleepsamp==2 | sleepsamp==6));
upperCrit = 5;     
lowerCrit = 2;
upperThreshold = upperCrit * stdDevRMS;
lowerThreshold = lowerCrit * stdDevRMS;

%% Find where RMS signal crosses threshold
aboveThreshold = rmsSignal > upperThreshold & (sleepsamp==2 | sleepsamp==6);
belowThreshold = rmsSignal < lowerThreshold & (sleepsamp==2 | sleepsamp==6);
% clear rmsSignal

%% Find the beginning and end of each upper threshold crossing:
startUpper = find(diff(aboveThreshold)==1)+1;
endUpper = find(diff(aboveThreshold)==-1);
if endUpper(1) < startUpper(1)                                          % Corrects for ripples that start at the beginning of scored period
    startUpper = [0; startUpper];
end
if startUpper(end) > endUpper(end)                                      % Corrects for ripples that end at the end of scored period
    endUpper = [endUpper; lengthSignal];
end

%% Remove segments at beginning and/or end if at the start and/or end of the recording:
if isequal(startUpper,0)
    startUpper(1)=[];
    endUpper(1)=[];
end

if isequal(endUpper,lengthSignal)
    startUpper(end)=[];
    endUpper(end)=[];
end

%% Find the beginning and end of each lower threshold crossing:
startLower = find(diff(belowThreshold)==1)+1;
endLower = find(diff(belowThreshold)==-1);
if endLower(1) < startLower(1)                                          % Corrects for ripples that start at the beginning of scored period
    startLower = [0; startLower];
end
if startLower(end) > endLower(end)                                      % Corrects for ripples that end at the end of scored period
    endLower = [endLower; lengthSignal];
end 

%% Find any upper crossings that are greater than the final lower crossing:
removeUpperIdx = find(endUpper > startLower(end));
if ~isempty(removeUpperIdx)
    startUpper(removeUpperIdx) = [];
    endUpper(removeUpperIdx) = [];
end

%% Use the nearestpoint.m function to find beginning and end of each ripple:
rippleStartIdx = endLower(nearestpoint(startUpper,endLower,'previous'));
rippleStopIdx = startLower(nearestpoint(endUpper,startLower,'next'));

%% Remove all replicates:
rippleIdx = unique([rippleStartIdx rippleStopIdx],'rows');

%% Combine close ripples:
minTimeBetween = 0.03;                                                      % Maximum distance between ripples to combine into same ripple
mindist = minTimeBetween * Fs;

distBetweenRipples = rippleIdx(2:end,1) - rippleIdx(1:(end-1),2);           % Find the distances between all ripples
k = (distBetweenRipples <= mindist);                                        % Find when the time between ripples is less than the maximum
clear distBetweenRipples maxdist
k1 = [k; 0];
k2 = [0; k];                                                                % Add a zero for the first ripple that was not included
clear k
rippleIdx(k1==1,2) = rippleIdx(k2==1,2);                                    % Unite them (use the ending of the second ripple as end of the first one...
rippleIdx(k2==1,:) = [];
clear k1 k2

%% Remove ripples that are too short or long in duration:
mindur = 0.02;                                                              % Minimum accepted ripple duration (20 msec)
minIdx = mindur * Fs;                                                       % Convert duration to # of indices
maxdur = 0.1;                                                               % Maximum accepted ripple duration (100 msec; Buzsaki and Zugaro)
maxIdx = maxdur * Fs;                                                       % Convert duration to # of indices
rippleDurIdx = rippleIdx(:,2) - rippleIdx(:,1) + 1;                             % Find the length of ripples in # of indices
rippleBounds = (rippleDurIdx < minIdx) | (rippleDurIdx > maxIdx);           % Find ripples lasting less than 20 msec or greater than 300 msec
rippleIdx(rippleBounds,:) = [];                                                   % Discard ripples that are too short or long 
clear minIdx maxIdx rippleDurIdx rippleBounds

%% Find the maximum RMS within each ripple:
ripple.idxMaxRMS = zeros(size(rippleIdx,1),1);
ripple.maxRMS = zeros(size(rippleIdx,1),1);
for i = 1:size(rippleIdx,1)
    [ripple.maxRMS(i),idxM] = max(rmsSignal(rippleIdx(i,1):rippleIdx(i,2)));
    ripple.idxMaxRMS(i) = rippleIdx(i,1) + idxM - 1;
end
clear idxM

%% Assign state to each ripple:
ripple.state = sleepsamp(ripple.idxMaxRMS);
clear sleepsamp

%% Find the ripple start and end time stamps:
ripple.TS = zeros(size(rippleIdx));
ripple.TS = [timeStamps(rippleIdx(:,1)) timeStamps(rippleIdx(:,2))];

ripple.idx = rippleIdx;
clear rippleIdx

%% SAVE DATA AND FIGURE
%Request user to name output file:
prompt = {'Enter the filename you want to save it as: (just the name)'};
def = {'Rat#_Day'};
dlgTitle = 'Save .MAT file';
lineNo = 1;
answer = inputdlg(prompt,dlgTitle,lineNo,def);
filename = char(answer(1,:));
save(fullfile(CSCFilePath,['Ripples', filename, '.mat']), 'scoredFile', 'CSCFilename', 'rippleBand',...
    'windowSeconds', 'upperCrit', 'lowerCrit', 'minTimeBetween', 'mindur',...
    'maxdur', 'ripple');

