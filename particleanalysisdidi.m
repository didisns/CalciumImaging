% Hopefully this can be used to analyze the caclium imaging data. 
% February 2019 by Didi

% let's set all the variables obtained until now from now.

% First, after you obtained all the ROI particle and pasted in excel:
% import it using matlab 'import data'. 
name = test2; % fill out the name 
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
[frame, column_roi] = size(test2);
frames = 1:frame;
timeframes = (frames-1)*imagingperiod;

% Then set the number of ROIs
numberROI = column_roi/5;
wholefieldROIs = 1;
greenROIs = 2:greennumber+1;
redROIs = greennumber+2:numberROI;

% find the real up state start and end times
NUS = length(tmp);
USrawstart = zeros(NUS,1);
USrawend = zeros(NUS,1);
for i = 1:NUS
    st = tmp{i,2};
    en = tmp{i,3};
    USrawstart(i) = st;
    USrawend(i) = en;
end


