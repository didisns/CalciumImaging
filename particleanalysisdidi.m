
% February 2019 by Didi. In order to analyze the calcium imaging data

% let's set all the variables obtained until now from now.

% First, after you obtained all the ROI particle and pasted in excel: 
[~,~,excel] = xlsread('C:\Users\Dieds\Google Drive\PhD PCDH19\data\calcium imaging mosaic\test2.xlsx'); % fill out the name 
greennumber = 1; % fill out the number of green ROIs, make sure they are second in the excel file
rednumber = 7; % fill out the number of red ROIs, make sure they are third in the excel file
imagingperiod = 0.12415810148978500000;

%then you need the LFP data from Zebra
LFPstartgalvo = 3.6865; % indicate the start time you obtained in the eventdetection tab 'save start imaging'
LFPstopgalvo = 189.9205; % indicate the stop time you obtained in the eventdetection tab 'save start imaging'
load('H:\Data to be analyzed\Calcium imaging\26-11-18 calcium imaging with LFP\UStable_26-11-18 calcium imaging with LFP_lfp_0008.abf.mat');
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
    USrawstart(i) = tmp{i,2};
    USrawend(i) = tmp{i,3};    
end

% then we have to remove those up states that were outside of the imaging
% interval. I do this by creating a vector called validUS that contains
% info on whether or not the up state is valid

validUS = zeros(NUS,1);
for i = 1:NUS
    if USrawend(i) <= LFPstartgalvo || USrawstart(i) >= LFPstopgalvo
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

% Now we can set the correct start time in terms of time in seconds, all
% nonvalid US have times of 0
USstart = zeros(NUS, 1);
USend = zeros(NUS, 1);
for i = firstUS:lastUS
    USstart(i) = USrawstart(i) - LFPstartgalvo;
    USend(i) = USrawend(i) - LFPstartgalvo;
end

% Final part of setting up: the output struct
    
output.wFOV = struct('total area',1,'time active', 1, 'number transients',1,'total area transients',1,...
    'duration transients',1, 'peak transients',1, 'percentage in up states',1);

output.red = struct('total area',1,'time active', 1, 'number transients',1,'total area transients',1,...
    'duration transients',1, 'peak transients',1, 'percentage in up states',1);

output.green = struct('total area',1,'time active', 1, 'number transients',1,'total area transients',1,...
    'duration transients',1, 'peak transients',1, 'percentage in up states',1);

% Done with all the setting up.

% Start with findings 

%% First part: total area involved in calcium transients. 
%# pixels per minute of imaging in whole field of view or per cell.Each cell is kept
% separately to look at later (otherwhise different amount of cells red and
% green cells in the same ROI mess up the statistics)

% whole field of view
minutesimaging = timeframes(end)/60;
areawFOV = sum([excel{:,3}])/minutesimaging;

%  For green cells
areagreen = zeros(greennumber,1);
for i = greenROIs    
    cnumber = ((i-1)*5)+3;
    areagreen(i) = sum([excel{:,cnumber}])/minutesimaging;
end

%  For red cells
areared = zeros(rednumber,1);
for i = redROIs    
    cnumber = ((i-1)*5)+3;
    areared(i) = sum([excel{:,cnumber}])/minutesimaging;
end

%% Second part: total time involved in calcium imaging

% For whole field of view
timewFOV = find([excel{:,3}]>0);
perc_timewFOV = (length(timewFOV)/frame)*100;

% for green cells
timegreen = zeros(greennumber,1);
perc_timegreen = zeros(greennumber,1);
for i = greenROIs
    cnumber = ((i-1)*5)+3;
    timegreen(i) = find([excel{:,cnumber}]>0);  
    perc_timegreen(i) = (length(timegreen(i))/frame)*100;
end

% for red cells
timered = zeros(rednumber,1);
perc_timered = zeros(rednumber,1);
for i = redROIs
    cnumber = ((i-1)*5)+3;
    timered(i) = find([excel{:,cnumber}]>0);
    perc_timered(i) = (length(timered(i))/frame)*100;
end

%% Third part: find the number of transients, their area, duration, and peak

% for whole field of view

for i=2:length(timewFOV)
    if notranswFOV
    

