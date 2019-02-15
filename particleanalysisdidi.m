% Hopefully this can be used to analyze the calcium imaging data. 
% February 2019 by Didi

% let's set all the variables obtained until now from now.

% First, after you obtained all the ROI particle and pasted in excel: 
[~,~,excel] = xlsread('C:\Users\Dieds\Google Drive\PhD PCDH19\data\calcium imaging mosaic\test2.xlsx'); % fill out the name 
greennumber = 1; % fill out the number of green ROIs, make sure they are second in the excel file
rednumber = 7; % fill out the number of red ROIs, make sure they are third in the excel file
imagingperiod = 0.12415810148978500000;

%then you need the LFP data from Zebra
LFPstartgalvo = 3.6865; % indicate the start time you obtained in the eventdetection tab 'save start imaging'
LFPstopgalvo = 189.9205; % indicate the stop time you obtained in the eventdetection tab 'save start imaging'
load('L:\Data to be analyzed\Calcium imaging\26-11-18 calcium imaging with LFP\UStable_26-11-18 calcium imaging with LFP_lfp_0008.abf.mat');
% here add the path to the saved table from eventdetection

% check if import table up states is correct
if isempty(tmp)
    error('US table not loaded properly')
end

% Now let's set the time of each frame 
[frame, column_roi] = size(excel);
frames = 1:frame;
timeframes = (frames-1)*imagingperiod; % a vector where each row contains the time when this frame was taken

% Then set the number of ROIs
numberROI = column_roi/5;
wholefieldROIs = 1;
greenROIs = 2:greennumber+1;
redROIs = greenROIs(end)+1:numberROI;

% test if okay
if length(greenROIs) ~= greennumber || length(redROIs) ~= rednumber
    error('you calculated the green and red ROIs wrong');
end

% find the up state start and end times as determined in eventdetection,
% from the start of the ephys recording
NUS = length(tmp);
USrawstart = zeros(NUS,1);
USrawend = zeros(NUS,1);
for i = 1:NUS
    st = tmp{i,2};
    en = tmp{i,3};
    USrawstart(i) = st;
    USrawend(i) = en;
end

% then we have to remove those up states that were outside of the imaging
% interval. I do this by creating a vector called validUS that contains
% info on whether or not the up state is valid

validUS = zeros(NUS,1);
for i = 1:NUS
    if USrawend(i) <= LFPstartgalvo || USrawstart >= LFPstopgalvo
        validUS(i) = 0;
    else
        validUS(i) = 1;
    end
end
      
% find the first and last up states, if they are on the edge of the imaging
% time, change the start/end time to start/end time of imaging. 
validUS2 = find(validUS > 0);
NUSvalid = length(validUS2);
firstUS = validUS2(1);
lastUS = validUS2(end);
if tmp{firstUS, 2} < LFPstartgalvo
    USrawstart(firstUS) = LFPstartgalvo;
end
if tmp{lastUS, 3} > LFPstopgalvo
    USrawend(lastUS) = LFPstopgalvo;
end

% Now we can set the correct start time in terms of time in seconds
USstart = zeros(NUS, 1);
USend = zeros(NUS, 1);
for i = firstUS:lastUS
    USstart(i) = USrawstart - LFPstartgalvo;
    USend(i) = USrawend - LFPstartgalvo;
end
    
% Done with all the setting up.

% Start with findings 