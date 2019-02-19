
% February 2019 by Didi. In order to analyze the calcium imaging data

% let's set all the variables obtained until now from now.

% First, after you obtained all the ROI particle and pasted in excel: 
[~,~,excel] = xlsread('C:\Users\SNS\Google Drive\PhD PCDH19\data\calcium imaging mosaic\test2.xlsx'); % fill out the name 
greennumber = 1; % fill out the number of green ROIs, make sure they are second in the excel file
greenROIarea = [166]; % fill out the area's of the green ROIs measured in imagej ROI manager
redROIarea = [151, 140, 133, 112, 112, 140, 151]; % fill out the area's of the red ROIs measured in imagej ROI manager
rednumber = 7; % fill out the number of red ROIs, make sure they are third in the excel file
imagingperiod = 0.12415810148978500000; % this I always use, but always check whether correct
umpixel = 1.18803166994996; % this I always use, but always check whether correct
pixelnumber = 256; % this I always use, but always check whether correct

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
    
output.wFOV = struct('area',1,'time', 1, 'number_transients',1,'area_transients',1,...
    'duration_transients',1, 'peak_transients',1, 'percentage_up_states',1);

output.green = struct('area',1,'time', 1, 'number_transients',1,'area_transients',1,...
    'duration_transients',1, 'peak_transients',1, 'percentage_up_states',1);

output.red = struct('area',1,'time', 1, 'number_transients',1,'area_transients',1,...
    'duration_transients',1, 'peak_transients',1, 'percentage_up_states',1);

% Done with all the setting up.

% Start with findings 

%% First part: permil area involved in calcium transients. 
% permil pixels involved in calcium activity of imaging in whole field of view or per cell.Each cell is kept
% separately to look at later (otherwhise different amount of cells red and
% green cells in the same ROI mess up the statistics)

% whole field of view
totalareawFOV = sum([excel{:,3}]);
areawFOV = (totalareawFOV/(pixelnumber*pixelnumber*frame))*1000;

%  For green cells
totalareagreen = zeros(greennumber, 1);
areagreen = zeros(greennumber, 1);
for i = 1:greennumber
    ROInumber = greenROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    totalareagreen(i) = sum([excel{:,cnumber}]);  
    areagreen(i) = (totalareagreen(i)/(greenROIarea(i)*frame))*1000;
end

%  For red cells
totalareared = zeros(rednumber, 1);
areared = zeros(rednumber, 1);
for i = 1:rednumber
    ROInumber = redROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    totalareared(i) = sum([excel{:,cnumber}]);
    areared(i) = (totalareared(i)/(redROIarea(i)*frame))*1000;
end

% Now add it all to the output struct

output.wFOV.area = areawFOV;
output.green.area = areagreen;
output.red.area = areared;

%% Second part: percentage time involved in calcium imaging

% For whole field of view
timewFOV = find([excel{:,3}]>0);
perc_timewFOV = (length(timewFOV)/frame)*100;

% for green cells
perc_timegreen = zeros(greennumber,1);
for i = 1:greennumber
    ROInumber = greenROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    timegreen = find([excel{:,cnumber}]>0);  
    perc_timegreen(i) = (length(timegreen)/frame)*100;
end

% for red cells
perc_timered = zeros(rednumber,1);
for i = 1:rednumber
    ROInumber = redROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    timered = find([excel{:,cnumber}]>0);
    perc_timered(i) = (length(timered)/frame)*100;
end

% Put it in output struct
output.wFOV.time = perc_timewFOV;
output.green.time = perc_timegreen;
output.red.time = perc_timered;

%% Third part: find the number of transients, their area, duration, and peak

% for whole field of view

% count the number of transients and register their start time

countwFOV = 0;
if timewFOV(1) == 1
    countwFOV = countwFOV+1;
    savei(countwFOV) = 1; % needed to determine the end time in the next loop
    starttransients_wFOV(countwFOV) = 0;
end

for i=2:length(timewFOV)
    if timewFOV(i)-1 ~= timewFOV(i-1);
        countwFOV = countwFOV+1;
        savei(countwFOV) = i;
        currentframe = timewFOV(i);
        starttransients_wFOV(countwFOV) = timeframes(currentframe);
    end
end

% now find their end times 
endtransients_wFOV = zeros(countwFOV, 1);

for i = 1:(countwFOV-1)
    indexend = savei(i+1)-1;
    currentframe = timewFOV(indexend);
    endtransients_wFOV(i) = timeframes(currentframe);
end

% then the end of the last transient is the last value of timewFOV
currentframe = timewFOV(end);
endtransients_wFOV(countwFOV) = timeframes(currentframe);

% Then calculate metrics of transients
duration_wFOV = zeros(countwFOV);
areatr_wFOV = zeros(countwFOV);
peak_wFOV = zeros(countwFOV);
for i = 1:countwFOV
    duration_wFOV(i) = endtransients_wFOV(i)-starttransients_wFOV(i);
end
    
    
    

