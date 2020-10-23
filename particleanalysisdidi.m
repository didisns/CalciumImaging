
% February 2019 by Didi. Program to analyze calcium transients obtained
% from ImageJ, and correspond them to up states detected in the LFP trace. 

% The program takes as an input an excel file containing the particles calculated using ImageJ's
% "Analyze Particles". It contains on each column information for one ROI, each row represents one
% time point, a 0 denotes no particle detected, a 1 that a particle has been detected.
% The input excel file requires the particles counted in the entire field of view
% in the first column, then the green ROI(s) in the following column(s), the red
% ROI(s) in the last column(s). No header info should be provided.

% THe second input file is the UP state table, which is the output of ZebraExplore's 
% eventdetection module.  

% The program outputs some quantifications in terms of the calcium 'transients'.
% Results are copied to the clipboard, which I pasted in excel files. 

% Output:
% Area: 	in cubic micron, average per frame
% Areapermille:	permille of all pixels involved in calcium acitivity, for ROIs based on area ROI
% Time:	percentage of frames where there is some form of calcium activity
% Number_transients:	Transients (framenumber with calcium acivity) in Hz (SO NOT INIDIVDUAL SPOTS!!)
% Duration_transients:	Duration in seconds of above-mentioned transients
% Area_transients:	Sum all pixels in transients in um
% Areapermille_transients:	permille of all pixels in transients involved in calcium acitivity, for ROIs based on area ROI
% Peak_transients:	Highest value of transients, in um
%Peakpermille_transients:	Highest value of transients, in permille pixels, for ROIs based on area ROI
% Percentage_up_states:	Percentage pixels in binary file that fall within up states


%% Import information

% Reading of the excel input file and other information
[~,~,excel] = xlsread('L:\Data to be analyzed\Calcium imaging\20190803\lfp018.xlsx'); % fill out the name 
greennumber = 5; % the number of green ROIs, make sure they are second in the excel file
greenROIarea = [137, 177, 128, 140, 104]; % the area's of the green ROIs measured in imagej ROI manager
redROIarea = [219, 137, 136, 102, 102]; % the area's of the red ROIs measured in imagej ROI manager
rednumber = 5; % the number of red ROIs, make sure they are third in the excel file
imagingperiod = 0.12415810148978500000; % period of one frame (in seconds)
umpixel = 1.18803166994996; % um per pixel of the imaging data
pixelnumber = 256; % always a square: i.e. 356*256 pixels.

%then the LFP data from ZebraExplore
LFPstartgalvo = 3.6035; % indicate the start time obtained in the eventdetection tab 'save start imaging'
LFPstopgalvo = 189.8365; % indicate the stop time obtained in the eventdetection tab 'save start imaging'
load('L:\Data to be analyzed\Calcium imaging\20190803\UStable_20190803_lfp_0018.abf.mat');
% here add the path to the saved table from eventdetection

%%  

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
[NUS value] = size(tmp);
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
    
output.wFOV = struct('area',1, 'areapermille', 1, 'time', 1, 'number_transients',1,'area_transients',1,...
    'areapermille_transients', 1, 'duration_transients',1, 'peak_transients',1, ...
    'peakpermille_transients', 1, 'percentage_up_states',1);

output.green = struct('area',1, 'areapermille', 1, 'time', 1, 'number_transients',1,'area_transients',1,...
    'areapermille_transients', 1, 'duration_transients',1, 'peak_transients',1, ...
    'peakpermille_transients', 1, 'percentage_up_states',1);

output.red = struct('area',1, 'areapermille', 1, 'time', 1, 'number_transients',1,'area_transients',1,...
    'areapermille_transients', 1, 'duration_transients',1, 'peak_transients',1, ...
    'peakpermille_transients', 1, 'percentage_up_states',1);

% Done with all the setting up.

% Start with findings 

%% First part: permil area involved in calcium transients. 
% permil pixels involved in calcium activity of imaging in whole field of view or per cell.Each cell is kept
% separately to look at later (otherwhise different amount of cells red and
% green cells in the same ROI mess up the statistics)

% whole field of view
totalareawFOV = sum([excel{:,3}]);
areawFOV_permille = (totalareawFOV/(pixelnumber*pixelnumber*frame))*1000;
areawFOV = (totalareawFOV/frame)*(umpixel*umpixel);

%  For green cells
totalareagreen = zeros(greennumber, 1);
areagreen = zeros(greennumber, 1);
areagreen_permille = zeros(greennumber, 1);
for i = 1:greennumber
    ROInumber = greenROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    totalareagreen(i) = sum([excel{:,cnumber}]);    
    areagreen_permille(i) = (totalareagreen(i)/(greenROIarea(i)*frame))*1000;
    areagreen(i) = (totalareagreen(i)/frame)*(umpixel*umpixel);
end

%  For red cells
totalareared = zeros(rednumber, 1);
areared = zeros(rednumber, 1);
areared_permille = zeros(rednumber, 1);
for i = 1:rednumber
    ROInumber = redROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    totalareared(i) = sum([excel{:,cnumber}]);
    areared_permille(i) = (totalareared(i)/(redROIarea(i)*frame))*1000;
    areared(i) = (totalareared(i)/frame)*(umpixel*umpixel);
end

% Now add it all to the output struct

output.wFOV.area = areawFOV;
output.wFOV.areapermille = areawFOV_permille;
output.green.area = areagreen;
output.green.areapermille = areagreen_permille;
output.red.area = areared;
output.red.areapermille = areared_permille;

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

%% Third part: find the number of transients (in Hz), their mean area (as permille), mean duration (in s), and mean peak (in permille): wFOV.

% for whole field of view

% count the number of transients and register their start time
countwFOV = 1;        
savei(countwFOV) = 1;
currentframe = timewFOV(1);
starttransients_wFOV(countwFOV) = timeframes(currentframe);
        
if length(timewFOV)>1 
    for i=2:length(timewFOV)
        if timewFOV(i)-1 ~= timewFOV(i-1);
            countwFOV = countwFOV+1;
            savei(countwFOV) = i;
            currentframe = timewFOV(i);
            starttransients_wFOV(countwFOV) = timeframes(currentframe);
        end
    end
end

% now find their end times 
endtransients_wFOV = zeros(countwFOV, 1);       
if countwFOV < 2
    currentframe = timewFOV(end);
    endtransients_wFOV(1) = timeframes(currentframe);
else
    for i = 1:(countwFOV-1)
        indexend = savei(i+1)-1;
        currentframe = timewFOV(indexend);
        endtransients_wFOV(i) = timeframes(currentframe);
    end
    
    % then the end of the last transient is the last value of timegreen
    currentframe = timewFOV(end);
    endtransients_wFOV(countwFOV) = timeframes(currentframe);
end

% Then calculate metrics of transients
duration_wFOV = zeros(countwFOV, 1);
areatr_wFOV = zeros(countwFOV, 1);
peak_wFOV = zeros(countwFOV, 1);
for i = 1:countwFOV
    duration_wFOV(i) = (endtransients_wFOV(i)-starttransients_wFOV(i))+imagingperiod;
    indexstart = find(timeframes==starttransients_wFOV(i));
    indexend = find(timeframes==endtransients_wFOV(i));
    areatr_wFOV(i) = sum([excel{indexstart:indexend,3}]);
    peak_wFOV(i) = max([excel{indexstart:indexend,3}]);
end
    
meanduration_wFOV = mean(duration_wFOV);
meandurationI_wFOV = meanduration_wFOV/imagingperiod;
meanareatr_wFOV = mean(areatr_wFOV);
meanareatr_wFOV_um = meanareatr_wFOV*(umpixel*umpixel);
meanareatr_permille_wFOV = (meanareatr_wFOV/(pixelnumber*pixelnumber*meandurationI_wFOV))*1000;
meanpeak_wFOV = mean(peak_wFOV);
meanpeak_wFOV_um = meanpeak_wFOV*(umpixel*umpixel);
meanpeak_permille_wFOV = (meanpeak_wFOV/(pixelnumber*pixelnumber))*1000;

% save in output struct
frequencywFOV = countwFOV/(imagingperiod*frame);

output.wFOV.number_transients = frequencywFOV;
output.wFOV.duration_transients = meanduration_wFOV;
output.wFOV.area_transients = meanareatr_wFOV_um;
output.wFOV.areapermille_transients = meanareatr_permille_wFOV;
output.wFOV.peak_transients = meanpeak_wFOV_um;
output.wFOV.peakpermille_transients = meanpeak_permille_wFOV;

%%  Fourth part: find the number of transients, their area, duration, and peak. Green Cells.

countgreen = zeros(greennumber, 1);
frequencygreen = zeros(greennumber, 1);
meanduration_green = zeros(greennumber, 1);
meandurationI_green = zeros(greennumber, 1);
meanareatr_green = zeros(greennumber, 1);
meanareatr_green_um = zeros(greennumber, 1);
meanareatr_permille_green = zeros(greennumber, 1);
meanpeak_green = zeros(greennumber, 1);
meanpeak_green_um = zeros(greennumber, 1);
meanpeak_permille_green = zeros(greennumber, 1);

for gr = 1:greennumber
    % first reset some vectors to make sure they don't end up wrong
    savei = [];
    starttransients_green = [];
    endtransients_green = [];
    
    ROInumber = greenROIs(gr);
    cnumber = ((ROInumber-1)*5)+3;
    timegreen = find([excel{:,cnumber}]>0);
    
    % if no transients at all
    if isempty (timegreen)
        countgreen(gr) = 0;
        frequencygreen (gr) = 0;
        meanduration_green(gr) = 0;
        meandurationI_green (gr) = 0;
        meanareatr_green(gr) = 0;
        meanareatr_green_um(gr) = 0;
        meanareatr_permille_green (gr) = 0;
        meanpeak_green(gr) = 0;
        meanpeak_green_um(gr) = 0;
        meanpeak_permille_green (gr) = 0;
    else    
        countgreen(gr) = 1;        
        savei(countgreen(gr)) = 1;
        currentframe = timegreen(1);
        starttransients_green(countgreen(gr)) = timeframes(currentframe);
        
        if length(timegreen)>1 
            for i=2:length(timegreen)
                if timegreen(i)-1 ~= timegreen(i-1);
                    countgreen(gr) = countgreen(gr)+1;
                    savei(countgreen(gr)) = i;
                    currentframe = timegreen(i);
                    starttransients_green(countgreen(gr)) = timeframes(currentframe);
                end
            end
        end
    
        % now find their end times 
        endtransients_green = zeros(countgreen(gr), 1);       
        if countgreen(gr)<2
            currentframe = timegreen(end);
            endtransients_green(1) = timeframes(currentframe);
        else
            for i = 1:(countgreen(gr)-1)
                indexend = savei(i+1)-1;
                currentframe = timegreen(indexend);
                endtransients_green(i) = timeframes(currentframe);
            end
    
            % then the end of the last transient is the last value of timegreen
            currentframe = timegreen(end);
            endtransients_green(countgreen(gr)) = timeframes(currentframe);
        end

        % Then calculate metrics of transients
        duration_green = zeros(countgreen(gr),1);
        areatr_green = zeros(countgreen(gr),1);
        peak_green = zeros(countgreen(gr),1);
        for i = 1:countgreen(gr)
            duration_green(i) = (endtransients_green(i)-starttransients_green(i))+imagingperiod;
            indexstart = find(timeframes==starttransients_green(i));
            indexend = find(timeframes==endtransients_green(i));
            areatr_green(i) = sum([excel{indexstart:indexend,cnumber}]);
            peak_green(i) = max([excel{indexstart:indexend,cnumber}]);
        end
        
        frequencygreen (gr) = countgreen (gr)/(imagingperiod*frame);
        meanduration_green(gr) = mean(duration_green);        
        meandurationI_green (gr) = meanduration_green(gr)/imagingperiod;       
        meanareatr_green(gr) = mean(areatr_green);
        meanareatr_green_um(gr) = meanareatr_green(gr)*(umpixel*umpixel);
        meanareatr_permille_green (gr) = (meanareatr_green(gr)/(greenROIarea(gr)*meandurationI_green (gr)))*1000;
        meanpeak_green(gr) = mean(peak_green);
        meanpeak_green_um(gr) = meanpeak_green(gr)*(umpixel*umpixel);
        meanpeak_permille_green (gr) = (meanpeak_green(gr)/greenROIarea(gr))*1000;
    end     
end

% save in output struct


output.green.number_transients = frequencygreen;
output.green.duration_transients = meanduration_green;
output.green.area_transients = meanareatr_green_um;
output.green.areapermille_transients = meanareatr_permille_green;
output.green.peak_transients = meanpeak_green_um;
output.green.peakpermille_transients = meanpeak_permille_green;

%%  Fifth part: find the number of transients, their area, duration, and peak. Red Cells.

countred = zeros(rednumber, 1);
frequencyred = zeros(rednumber, 1);
meanduration_red = zeros(rednumber, 1);
meandurationI_red = zeros(rednumber, 1);
meanareatr_red = zeros(rednumber, 1);
meanareatr_red_um = zeros(rednumber, 1);
meanareatr_permille_red = zeros(rednumber, 1);
meanpeak_red = zeros(rednumber, 1);
meanpeak_red_um = zeros(rednumber, 1);
meanpeak_permille_red = zeros(rednumber, 1);

for rr = 1:rednumber
    % first reset some vectors to make sure they don't end up wrong
    savei = [];
    starttransients_red = [];
    endtransients_red = [];
    
    ROInumber = redROIs(rr);
    cnumber = ((ROInumber-1)*5)+3;
    timered = find([excel{:,cnumber}]>0);
    
    % if no transients at all
    if isempty (timered)
        countred(rr) = 0;
        frequencyred (rr) = 0;
        meanduration_red(rr) = 0;
        meandurationI_red (rr) = 0;
        meanareatr_red(rr) = 0;
        meanareatr_red_um(rr) = 0;
        meanareatr_permille_red (rr) = 0;
        meanpeak_red(rr) = 0;
        meanpeak_red_um(rr) = 0;
        meanpeak_permille_red (rr) = 0;
        
    else
        countred(rr) = 1;        
        savei(countred(rr)) = 1;
        currentframe = timered(1);
        starttransients_red(countred(rr)) = timeframes(currentframe);
        
        if length(timered)>1 
            for i=2:length(timered)
                if timered(i)-1 ~= timered(i-1);
                    countred(rr) = countred(rr)+1;
                    savei(countred(rr)) = i;
                    currentframe = timered(i);
                    starttransients_red(countred(rr)) = timeframes(currentframe);
                end
            end
        end    
           
        % now find their end times
        endtransients_red = zeros(countred(rr), 1);       
        if countred(rr)<2
            currentframe = timered(end);
            endtransients_red(1) = timeframes(currentframe);
        else
            for i = 1:(countred(rr)-1)
                indexend = savei(i+1)-1;
                currentframe = timered(indexend);
                endtransients_red(i) = timeframes(currentframe);
            end
    
            % then the end of the last transient is the last value of
            % timered
            currentframe = timered(end);
            endtransients_red(countred(rr)) = timeframes(currentframe);
        end
        
        % Then calculate metrics of transients
        duration_red = zeros(countred(rr),1);
        areatr_red = zeros(countred(rr),1);
        peak_red = zeros(countred(rr),1);
        for i = 1:countred(rr)
            duration_red(i) = (endtransients_red(i)-starttransients_red(i))+imagingperiod;
            indexstart = find(timeframes==starttransients_red(i));
            indexend = find(timeframes==endtransients_red(i));
            areatr_red(i) = sum([excel{indexstart:indexend,cnumber}]);
            peak_red(i) = max([excel{indexstart:indexend,cnumber}]);
        end
        
        frequencyred (rr) = countred (rr)/(imagingperiod*frame);
        meanduration_red(rr) = mean(duration_red);        
        meandurationI_red(rr) = meanduration_red(rr)/imagingperiod;        
        meanareatr_red(rr) = mean(areatr_red);
        meanareatr_red_um(rr) = meanareatr_red(rr)*(umpixel*umpixel);
        meanareatr_permille_red(rr) = (meanareatr_red(rr)/(redROIarea(rr)*meandurationI_red (rr)))*1000;
        meanpeak_red(rr) = mean(peak_red);
        meanpeak_red_um(rr) = meanpeak_red(rr)*(umpixel*umpixel);
        meanpeak_permille_red (rr) = (meanpeak_red(rr)/redROIarea(rr))*1000;
    end    
end

% save in output struct

output.red.number_transients = frequencyred;
output.red.duration_transients = meanduration_red;
output.red.area_transients = meanareatr_red_um;
output.red.areapermille_transients = meanareatr_permille_red;
output.red.peak_transients = meanpeak_red_um;
output.red.peakpermille_transients = meanpeak_permille_red;

%% Finally: determine percentage of total pixels occuring in up states

% for wFOV
totalareawFOV = sum([excel{:,3}]); % already calculated but better to be sure

us_area_wFOV = 0;
for us = validUS2(1):validUS2(end)
    [value, idxst_us] = min(abs(timeframes-USstart(us)));
    [value2, idxen_us] = min(abs(timeframes-USend(us)));
    currentarea = sum([excel{idxst_us:idxen_us,3}]);
    us_area_wFOV = us_area_wFOV + currentarea;    
end

perc_us_wFOV = (us_area_wFOV/totalareawFOV)*100;

% For green cells

totalareagreen = zeros(greennumber, 1);
us_area_green = zeros(greennumber, 1);
perc_us_green = zeros(greennumber, 1);

for i = 1:greennumber
    ROInumber = greenROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    totalareagreen(i) = sum([excel{:,cnumber}]); 
    if totalareagreen(i) == 0
        perc_us_green(i) = NaN;
    else        
        for us = validUS2(1):validUS2(end)
            [value, idxst_us] = min(abs(timeframes-USstart(us)));
            [value2, idxen_us] = min(abs(timeframes-USend(us)));
            currentarea = sum([excel{idxst_us:idxen_us,cnumber}]);
            us_area_green(i) = us_area_green(i) + currentarea; 
        end 
        perc_us_green(i) = (us_area_green(i)/totalareagreen(i))*100;
    end
end

% For red cells

totalareared = zeros(rednumber, 1);
us_area_red = zeros(rednumber, 1);
perc_us_red = zeros(rednumber, 1);

for i = 1:rednumber
    ROInumber = redROIs(i);
    cnumber = ((ROInumber-1)*5)+3;
    totalareared(i) = sum([excel{:,cnumber}]);    
    if totalareared(i) == 0
        perc_us_red(i) = NaN;
    else
        for us = validUS2(1):validUS2(end)
            [value, idxst_us] = min(abs(timeframes-USstart(us)));
            [value2, idxen_us] = min(abs(timeframes-USend(us)));
            currentarea = sum([excel{idxst_us:idxen_us,cnumber}]);
            us_area_red(i) = us_area_red(i) + currentarea; 
        end 
        perc_us_red(i) = (us_area_red(i)/totalareared(i))*100;
    end
end

% Add everything to the output struct

output.wFOV.percentage_up_states = perc_us_wFOV;
output.green.percentage_up_states = perc_us_green;
output.red.percentage_up_states = perc_us_red;

%% Finally: let's copy stuff into the clipboard for easy workflow

% start with a first row of text
fields = fieldnames(output.wFOV);
column = length(fields)+1;
str = ' ';

% the second row will contain the wFOV 
for f = 1:column
    if f == 1
        str = sprintf('%s\t', str, 'wFOV');
    elseif f == column        
        str = sprintf('%s%f', str, output.wFOV.(fields{f-1}));
    else str = sprintf('%s%f\t', str, output.wFOV.(fields{f-1}));
    end
end
str = sprintf('%s\n', str);

% then the green cells

for row = 1:length(output.green.peak_transients)
    for f = 1:column        
        if f == 1
            str = sprintf('%s\t', str, 'green');
        elseif f == column
            vector = output.green.(fields{f-1});
            str = sprintf('%s%f', str, vector(row));
        else vector = output.green.(fields{f-1});
            str = sprintf('%s%f\t', str, vector(row));
        end
    end
    str = sprintf('%s\n', str);
end

% Finally the red cells
for row = 1:length(output.red.peak_transients)
    for f = 1:column        
        if f == 1
            str = sprintf('%s\t', str, 'red');
        elseif f == column
            vector = output.red.(fields{f-1});
            str = sprintf('%s%f', str, vector(row));
        else vector = output.red.(fields{f-1});
            str = sprintf('%s%f\t', str, vector(row));
        end
    end
    str = sprintf('%s\n', str);
end
    
        
clipboard ('copy',str);
'finished!'