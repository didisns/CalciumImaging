function varargout = eventDetectionMelies(varargin)
%EVENTDETECTIONMELIES M-file for eventDetectionMelies.fig
%      EVENTDETECTIONMELIES, by itself, creates a new EVENTDETECTIONMELIES or raises the existing
%      singleton*.
%
%      H = EVENTDETECTIONMELIES returns the handle to a new EVENTDETECTIONMELIES or the handle to
%      the existing singleton*.
%
%      EVENTDETECTIONMELIES('Property','Value',...) creates a new EVENTDETECTIONMELIES using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to eventDetectionMelies_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      EVENTDETECTIONMELIES('CALLBACK') and EVENTDETECTIONMELIES('CALLBACK',hObject,...) call the
%      local function named CALLBACK in EVENTDETECTIONMELIES.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eventDetectionMelies

% Last Modified by GUIDE v2.5 11-May-2018 10:00:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eventDetectionMelies_OpeningFcn, ...
                   'gui_OutputFcn',  @eventDetectionMelies_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end

% --- Executes just before eventDetectionMelies is made visible.
function eventDetectionMelies_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for eventDetectionMelies
handles.output = hObject;

% establish the link with the calling process
handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

% local structures and variables
%handles.USlist = struct('N',1,'tfrom',1,'tto',1,'delta',1);
handles.USlist = {1 1 1 1 1 1};

handles.upstates = struct('fromI',1,'fromT',1,'toI',1,'toT',1,'deltaT',1,'medianV',1,'SD',1,'pw',1,'MeanCa', 1, 'IntCa',1,'enable',1);
handles.downstates = struct('fromI',1,'fromT',1,'toI',1,'toT',1,'deltaT',1,'medianV',1,'SD',1,'pw',1,'enable',1);

% Update handles structure
guidata(hObject, handles);
end

% UIWAIT makes eventDetectionMelies wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = eventDetectionMelies_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function baselineThr_Callback(hObject, eventdata, handles)
end

function baselineThr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function minDuration_Callback(hObject, eventdata, handles)
end

function minDuration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function maxInterruption_Callback(hObject, eventdata, handles)
end

function maxInterruption_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function eventThr_Callback(hObject, eventdata, handles)
end

function eventThr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function extractEvents_Callback(hObject, eventdata, handles)

handles.parentHandles = guidata(handles.parentObject);    % this to refresh the local values of parentHandles

handles.SPG = [];
handles.LFP = [];
handles.timeSPG = [];
handles.upstates.fromI = [];
handles.upstates.toI = [];
handles.downstates.fromI = [];
handles.downstates.toI = [];
handles.USlist = [];
handles.USlist = {1 1 1 1 1 1 1};

% set the threshold for event detection based on the used method
if get(handles.autoThr,'Value')
    lowThr = handles.autoblThr;
    highThr = handles.autoevThr;
else
    lowThr = str2double(get(handles.baselineThr,'String'));
    highThr = str2double(get(handles.eventThr,'String'));
end

cCh = handles.parentHandles.wCh;    % Work channel identified by loadEphys
len = handles.parentHandles.dtaLen; % Lenght of the Ephys data
spp = handles.parentHandles.sp;

handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));

handles.LFP(1:len) = handles.parentHandles.LFP (cCh,1:len);

sl = length(handles.parentHandles.BPpwrTime);
handles.timeSPG = handles.parentHandles.BPpwrTime;
handles.SPG(1:sl) = handles.parentHandles.logBP(cCh,1:sl);

downStatesBolean = handles.SPG < lowThr;
downStatesI = find(handles.SPG< lowThr);
upStatesBolean = handles.SPG > highThr;
upStatesI = find(handles.SPG>highThr);
limboBolean = ~(downStatesBolean | upStatesBolean);

% upStatesI contains the indexes of the events. 
% This has to be converted to the indexes of the workLFP file

lus = length (upStatesI);
lds = length (downStatesI);

% process high energy states first
if lus>0
    firstI = upStatesI (1);
    firstTime = handles.timeSPG (firstI); 
    handles.upstates.fromT (1) = firstTime;
    handles.upstates.fromI (1) = 1+firstTime/spp;
    nus = 1;
    for i=2:lus
        nextI = upStatesI(i);
        if nextI==firstI+1
            % continue the old US
            firstI=nextI;
            if i == lus
                % we have reached the end of the track. close the US
                endTime = handles.timeSPG (firstI);
                handles.upstates.toT (nus) = endTime;
                handles.upstates.toI (nus) = int32(1+endTime/spp);
                handles.upstates.deltaT (nus) = endTime-handles.upstates.fromT (nus);
                handles.upstates.enable (nus) = 1;
            end
        else
            firstTime = handles.timeSPG (firstI);
            % beginning of a new US. Close the old US
            handles.upstates.toT (nus) = firstTime;
            handles.upstates.toI (nus) = int32(1+firstTime/spp);
            handles.upstates.deltaT (nus) = firstTime-handles.upstates.fromT (nus);
            handles.upstates.enable (nus) = 1;
            nus = nus + 1;
            % open new US
            handles.upstates.fromT (nus) = handles.timeSPG (nextI);
            handles.upstates.fromI (nus) =  int32(1+handles.timeSPG (nextI)/spp);
            firstI = nextI;
        end    
    end
end

nus = length (handles.upstates.toI);      % just to make sure...

% process low energy states
if nus>0
    firstI = downStatesI (1);
    firstTime = handles.timeSPG (firstI); 
    handles.downstates.fromT (1) = firstTime;
    handles.downstates.fromI (1) = 1+firstTime/spp;
    nds = 1;

    for i=2:lds
        nextI = downStatesI(i);
        if nextI==firstI+1
            % continue the old DS
            firstI=nextI;
        else
            firstTime = handles.timeSPG (firstI);
            % beginning of a new DS. Close the old DS
            handles.downstates.toT (nds) = firstTime;
            handles.downstates.toI (nds) = int32(1+firstTime/spp);
            handles.downstates.deltaT (nds) = firstTime-handles.downstates.fromT (nds);
            handles.downstates.enable (nds) = 1;
            nds = nds + 1;
            % open new DS
            handles.downstates.fromT (nds) = handles.timeSPG (nextI);
            handles.downstates.fromI (nds) =  int32(1+handles.timeSPG (nextI)/spp);
            firstI = nextI;
        end
    end
end    
nds = length(handles.downstates.toI);

guidata(hObject, handles);

% defragment DS and US across short state gap. The max size of the filled
% gap is defined in the GUI
% These are the rules: a brief interruption shorter than maxInterruption is
% filled in by assigning to the neighboring US that are fused in a longer
% one. Isolated US briefer than minDuration are attributed to the limbo.

maxGap = str2double(get(handles.maxInterruption,'String'));
minUS = str2double(get(handles.minDuration,'String'));

if get(handles.removeCk,'Value');
    % process US first
    usi = 2
    while usi <= nus
        % perform the fusions first since brief US might get fused together
        % to form a longer and legal US.
        interval = [handles.upstates.fromT(usi) handles.upstates.toT(usi-1)];
        distance = interval(1) - interval(2);
        if (distance<=maxGap), 
            % first fuse the two upstates
            handles.upstates.toT(usi-1) = handles.upstates.toT(usi);
            handles.upstates.toI(usi-1) = handles.upstates.toI(usi);
            handles.upstates.deltaT(usi-1) = handles.upstates.toT(usi-1) - handles.upstates.fromT(usi-1);
            % second remove the US pointed to by usi
            for k = usi:nus-1
                % shift all USs
                handles.upstates.toT(k) = handles.upstates.toT(k+1);
                handles.upstates.toI(k) = handles.upstates.toI(k+1); 
                handles.upstates.fromT(k) = handles.upstates.fromT(k+1);
                handles.upstates.fromI(k) = handles.upstates.fromI(k+1); 
                handles.upstates.deltaT(k) = handles.upstates.deltaT (k+1);
            end
            nus = nus - 1;
        else
            usi = usi+1;
        end    
    end
    % remove brief states
    for usi=1:nus
        if handles.upstates.deltaT(usi) <= minUS, handles.upstates.enable(usi) = 0; 
        end 
    end
    % remove the disabled states from the list
    cnt = nus;
    for usi=nus:-1:1
        if handles.upstates.enable(usi);
            % do nothing
        else
            cnt = cnt - 1;
            % shift down if the US is not the last one of the track
            if usi < nus
                handles.upstates.fromI (cnt:-1:usi) = handles.upstates.fromI (cnt+1:-1:usi+1);
                handles.upstates.toI (cnt:-1:usi) = handles.upstates.toI (cnt+1:-1:usi+1);
                handles.upstates.fromT (cnt:-1:usi) = handles.upstates.fromT (cnt+1:-1:usi+1);
                handles.upstates.toT (cnt:-1:usi) = handles.upstates.toT (cnt+1:-1:usi+1);
                handles.upstates.deltaT (cnt:-1:usi) = handles.upstates.deltaT (cnt+1:-1:usi+1);                
                handles.upstates.enable (cnt:-1:usi) = handles.upstates.enable (cnt+1:-1:usi+1);
            end
        end 
    end
    nus = cnt;

    % process DS
    for usi=2:nds
        interval = [handles.downstates.fromT(usi) handles.downstates.toT(usi-1)];
        distance = interval(1) - interval(2);
        if (distance<=maxGap), 
            % fuse the contiguous downstates
            handles.downstates.toT(usi-1) = handles.downstates.toT(usi);
            handles.downstates.toI(usi-1) = handles.downstates.toI(usi);
            handles.downstates.deltaT(usi-1) = handles.downstates.toT(usi-1) - handles.downstates.fromT(usi-1);
            handles.downstates.enable(usi) = 0;    % mark the US for future removal
        end    
    end

    for usi=1:nds
        if handles.downstates.deltaT(usi) <= minUS, handles.downstates.enable(usi) = 0; 
        end
    end

    cnt = nds;
    for usi=nds:-1:1
        if handles.downstates.enable(usi);
            % do nothing
        else
            cnt = cnt - 1;
            % shift down
            if usi < nds
                handles.downstates.fromI (cnt:-1:usi) = handles.downstates.fromI (cnt+1:-1:usi+1);
                handles.downstates.toI (cnt:-1:usi) = handles.downstates.toI (cnt+1:-1:usi+1);
                handles.downstates.fromT (cnt:-1:usi) = handles.downstates.fromT (cnt+1:-1:usi+1);
                handles.downstates.toT (cnt:-1:usi) = handles.downstates.toT (cnt+1:-1:usi+1);
                handles.downstates.deltaT (cnt:-1:usi) = handles.downstates.deltaT (cnt+1:-1:usi+1);                
                handles.downstates.enable (cnt:-1:usi) = handles.downstates.enable (cnt+1:-1:usi+1);                
            end
        end
    end
    nds = cnt;
end

handles.NUS = nus;
handles.NDS = nds;

% At this stage all low power (DS) ad high power events (US) have been extracted and we can compute the
% relative metrics.
% Compute median, PW and SD of each US and DS.
for usi=1:nds   % compute DS first since you need this to compute US size
   i1 = handles.downstates.fromI(usi);
   i2 = handles.downstates.toI(usi);
   handles.downstates.medianV(usi) = median(handles.LFP(i1:i2)); 
   handles.downstates.SD(usi) = std(handles.LFP(i1:i2)); 
   handles.downstates.pw(usi) = rms(handles.LFP(i1:i2)); 
end

DSsearch = 1;
dsi1 = handles.downstates.fromI(DSsearch);
dsi2 = handles.downstates.toI(DSsearch);
for usi=1:nus
   i1 = handles.upstates.fromI(usi);
   i2 = handles.upstates.toI(usi);
   handles.upstates.medianV(usi) = median(handles.LFP(i1:i2));
   handles.upstates.SD(usi) = std(handles.LFP(i1:i2));
   handles.upstates.pw(usi) = rms(handles.LFP(i1:i2));   
   % now search for the closest DS
   while dsi2 < i1 && DSsearch <= nds
       DSsearch = DSsearch + 1;
       if DSsearch <= nds
           dsi2 = handles.downstates.toI(DSsearch);
       end
   end
   % now DSsearch points to the DS immediatey after unless the data ends
   % with a US. A second exception is when the data starts with an US
   if DSsearch == 1     % data begins with US
       baseline = handles.downstates.medianV(1);
   else
       if DSsearch > nds    % data ends up with a US
           baseline = handles.downstates.medianV(nds);
       else
           % OK, this is a middle of the road US!
           baseline = (handles.downstates.medianV(DSsearch)+handles.downstates.medianV(DSsearch-1))/2;
       end
   end
   handles.USlist (usi,5) = num2cell(handles.upstates.pw(usi));   % WARNING: this is the power of the event, not the energy!!!
end

% Creation of the list of high energy events
%handles.USlist (1,1) = num2cell(1);
%handles.USlist (1,2) = num2cell(handles.upstates.fromT(1));
handles.USlist (1:nus,1) = num2cell(1:nus);
handles.USlist (1:nus,2) = num2cell(handles.upstates.fromT(1:nus));
handles.USlist (1:nus,3) = num2cell(handles.upstates.toT(1:nus));
delta(1:nus) = handles.upstates.toT(1:nus) - handles.upstates.fromT(1:nus);
handles.USlist (1:nus,4) = num2cell(delta(1:nus));  % Duration of the event
handles.USlist (1:nus,6) = num2cell(handles.upstates.pw(1:nus).*delta(1:nus));   % energy
set(handles.UStable,'Data',handles.USlist);

% Copy the high power state descriptors to the parent process variable
% space

handles.parentHandles.events = zeros(nus,10);   % Reset of the array and preallocation
handles.parentHandles.nEvents = nus;
handles.parentHandles.events(1:nus,1:3) = cell2mat(handles.USlist(1:nus,1:3));
handles.parentHandles.events(1:nus,4) = cell2mat(handles.USlist (1:nus,4));             % Duration of the event
handles.parentHandles.events(1:nus,5) = cell2mat(handles.USlist (1:nus,6));       % Energy of the event

guidata(handles.parentObject,handles.parentHandles);

% plot the states IDs
plotStates (hObject, handles)

guidata(hObject, handles);
end

%%
function plotStates (hObject, handles)

cCh = handles.parentHandles.wCh;    % Work channel identified by loadEphys
len = handles.parentHandles.dtaLen; % Lenght of the Ephys data
spp = handles.parentHandles.sp;
handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));

axes(handles.plotStates);
cla
hold on
time = 0:spp:spp*(len-1);
plot (time,handles.LFP,'black');

for usi=1:handles.NUS
   i1 = handles.upstates.fromI(usi);
   i2 = handles.upstates.toI(usi);
   plot(time(i1:i2),handles.LFP(i1:i2),'red'); 
end

for usi=1:handles.NDS
   i1 = handles.downstates.fromI(usi);
   i2 = handles.downstates.toI(usi);
   plot(time(i1:i2),handles.LFP(i1:i2),'green'); 
end
axis ([handles.tmin handles.tmax -inf inf])

axes(handles.statesTrack);
set(gca, 'XTick', []);
cla
hold on
plot (handles.timeSPG,handles.SPG);
axis ([handles.tmin handles.tmax -inf inf])
hold off
end
%%

function analyzeData_Callback(hObject, eventdata, handles)
% This function computes the distribution of the spectral power in the
% 'envelope' band width. The distribution is fitted by a double gaussian
% and the limits for the event detection are computed

handles.parentHandles = guidata(handles.parentObject);    % this to refresh the local values of parentHandles

handles.Hcenters = [];
handles.cdist = [];
handles.cdistN = [];
handles.Hnelements = [];
handles.normNelements = [];

% first, extract the data to be analyzed
sl = length(handles.parentHandles.BPpwrTime);
cCh = handles.parentHandles.wCh;
handles.timeSPG = handles.parentHandles.BPpwrTime;
handles.SPG(1:sl) = handles.parentHandles.logBP(cCh,1:sl);

leftH = str2double(get(handles.fromH, 'String'));
rightH = str2double(get(handles.toH, 'String'));
bins = str2double(get(handles.binNumber, 'String'));


% check about auto scaling or not
if leftH==0, [nelements,centers] = hist(temp,bins);
else
   xValues = leftH:(rightH-leftH)/bins:rightH;
   [nelements,centers] = hist(handles.SPG,xValues);
end

% compute the 1) normalized histo, 2) cumulative distribution, 3)
% normalized cumulative distribution
totalN = sum (nelements);
handles.normNelements = nelements/totalN;
for j=1:bins+1
    handles.cdist (j) = sum(nelements(1:j));
end
handles.cdistN = handles.cdist / totalN;

handles.Hcenters = centers;
handles.Hnelements = nelements;

% plot the histogram
axes (handles.histoPlot);
cla
bar(centers,nelements);
axis ([leftH rightH 0 inf], 'auto y');
handles.meanPW(cCh) = mean(handles.SPG);
handles.stdPW(cCh) = std(handles.SPG);
nsamples = length(centers);
nmeasures = sum(nelements);

tll = handles.meanPW(cCh)-handles.stdPW(cCh)/3;
set(handles.baselineThr,'String',num2str(tll));
tll = handles.meanPW(cCh)+handles.stdPW(cCh)/3;
set(handles.eventThr,'String',num2str(tll));

%f = fit(centers.',nelements.','gauss2');

% fit the histogram with a double gaussian. 10.0266 = 2*2*sqr(2*pi)
param(2)= handles.meanPW(cCh)-handles.stdPW(cCh); % mean first component
param(3)= handles.stdPW(cCh)/2;                      % SD first component
amplitude(1)= nmeasures/(10.0266*param(3));                 % amplitude first component
param(5)= handles.meanPW(cCh)+handles.stdPW(cCh); % mean second component
param(6)= handles.stdPW(cCh)/2;                      % SD second component
amplitude(2)= nmeasures/(10.0266*param(6));                 % amplitude first component
gaussFnc = [];

% July 10, 2017. First, fit with fixed means and STDs to obtain a
% preliminary estimates of the gussian amplitudes.
ampliOut = fminsearch(@gaussFit1,amplitude);
param(1) = ampliOut(1);
param(4) = ampliOut(2);
% next, recompute everything.
paramOut = fminsearch(@gaussFit,param);

% paramOut contains the 6 parameters that define the best fit.
st = ['Gauss 1; mean:' num2str(paramOut(2),'%2.2f') ' ±SD: ' num2str(paramOut(3),'%2.2f')];
set (handles.gauss1param, 'String', st);
st = ['Gauss 2; mean:' num2str(paramOut(5),'%2.2f') ' ±SD: ' num2str(paramOut(6),'%2.2f')];
set (handles.gauss2param, 'String', st);

hold on
plot (centers,gaussFnc,'red','LineWidth',2);
yl = ylim;
line([paramOut(2)+paramOut(3) paramOut(2)+paramOut(3)],[yl(1) yl(2)]);
line([paramOut(5)-paramOut(6) paramOut(5)-paramOut(6)],[yl(1) yl(2)]);

% now compute the state thresholds
xinterc = fzero(@gaussDiff,(paramOut(2)+paramOut(5))/2);

handles.autoblThr = xinterc;
handles.autoevThr = xinterc;

guidata(gcbo, handles);

    function res = gaussFit(a)
        % gmodel = a(1)*exp(-(centers-a(2))^2/2*a(3).^2) + a(4)*exp(-(centers-a(5))^2/2*a(6).^2)
        % nelements is the array that contains the histogram bin counts: this
        % is the function to fit
        gaussFnc = a(1)*exp(-(centers-a(2)).^2/(2*(a(3)^2))) + a(4)*exp(-(centers-a(5)).^2/(2*(a(6)^2)));
        res = sum((nelements-gaussFnc).^2);
    end
    function res = gaussFit1(a)
        % gmodel = a(1)*exp(-(centers-a(2))^2/2*a(3).^2) + a(4)*exp(-(centers-a(5))^2/2*a(6).^2)
        % nelements is the array that contains the histogram bin counts: this
        % is the function to fit
        gaussFnc = a(1)*exp(-(centers-param(2)).^2/(2*(param(3)^2))) + a(2)*exp(-(centers-param(5)).^2/(2*(param(6)^2)));
        res = sum((nelements-gaussFnc).^2);
    end

    function y = gaussDiff(x)
        y = paramOut(1)*exp(-(x-paramOut(2)).^2/(2*(paramOut(3)^2)));
        y = y - paramOut(4)*exp(-(x-paramOut(5)).^2/(2*(paramOut(6)^2)));
    end
end

function binNumber_Callback(hObject, eventdata, handles)
end

function binNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function leftTime_Callback(hObject, eventdata, handles)
end

function leftTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function rightTime_Callback(hObject, eventdata, handles)
end

function rightTime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function replot_Callback(hObject, eventdata, handles)
plotStates (hObject, handles)
end

function erodeCk_Callback(hObject, eventdata, handles)
end

function removeCk_Callback(hObject, eventdata, handles)
end

function autoThr_Callback(hObject, eventdata, handles)
end

function toH_Callback(hObject, eventdata, handles)
end

function toH_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function fromH_Callback(hObject, eventdata, handles)
end

function fromH_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function traslateLeft_Callback(hObject, eventdata, handles)
handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));
delta = handles.tmax - handles.tmin;

set(handles.leftTime,'String',handles.tmin-delta);
set(handles.rightTime,'String',handles.tmin);
guidata(hObject, handles);

plotStates (hObject, handles)
end

function traslateRight_Callback(hObject, eventdata, handles)
handles.tmin = str2double(get(handles.leftTime,'String'));
handles.tmax = str2double(get(handles.rightTime,'String'));
delta = handles.tmax - handles.tmin;

set(handles.leftTime,'String',handles.tmax);
set(handles.rightTime,'String',handles.tmax+delta);
guidata(hObject, handles);

plotStates (hObject, handles)
end

function saveHCD_Callback(hObject, eventdata, handles)
% Save the distributions of the integrated spectral power
% March 22, 2017.

tmpOut1 = [];

handles.parentHandles = guidata(handles.parentObject);    % this to refresh the local values of parentHandles
cTr = handles.parentHandles.currentTrial;
cCh = handles.parentHandles.currentCh;

fileOut = [handles.parentHandles.dir_in 'PWSdistribution Ch' num2str(cCh) ' Trial' num2str(cTr) '.dat'];
tmpOut1 (1,:) = handles.Hcenters;        
tmpOut1 (2,:) = handles.Hnelements;    
tmpOut1 (3,:) = handles.normNelements;
tmpOut1 (4,:) = handles.cdist;    
tmpOut1 (5,:) = handles.cdistN;

% create the format descriptor
fmtSt=['%12.10f %12.10f %12.10f %12.10f %12.10f\n'];  % 5 columns and start new line
fid = fopen(fileOut,'w');
fprintf (fid,fmtSt,tmpOut1);
fclose(fid);

end
