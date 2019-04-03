function varargout = loadEPhys(varargin)
% LOADEPHYS MATLAB code for loadEPhys.fig
%      LOADEPHYS, by itself, creates a new LOADEPHYS or raises the existing
%      singleton*.
%
%      H = LOADEPHYS returns the handle to a new LOADEPHYS or the handle to
%      the existing singleton*.
%
%      LOADEPHYS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOADEPHYS.M with the given input arguments.
%
%      LOADEPHYS('Property','Value',...) creates a new LOADEPHYS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before loadEPhys_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to loadEPhys_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help loadEPhys

% Last Modified by GUIDE v2.5 14-May-2018 11:01:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @loadEPhys_OpeningFcn, ...
                   'gui_OutputFcn',  @loadEPhys_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before loadEPhys is made visible.
function loadEPhys_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to loadEPhys (see VARARGIN)

% Choose default command line output for loadEPhys
handles.output = hObject;
% establish link with the calling process
handles.parentObject =  varargin{1};
handles.parentHandles = varargin{2};

% disable all the controls that annot be operated until the EPhys data is
% loaded and properly analyzed
handles.displayLFP.Visible = 'off';
handles.logPW.Visible = 'off';
handles.displayLFP.Visible = 'off';

% prepare the frame selection slider
handles.selectFrame.Min = 1;
handles.selectFrame.Max = handles.parentHandles.metadata.frames;
handles.selectFrame.Value = 1;
handles.selectFrame.SliderStep = [0.01 0.1];    % expressed as % of entire range
handles.currentFrame = 1;
handles.currentFtxt.String = num2str(handles.currentFrame);

% set the trace duration
tlim = (handles.parentHandles.metadata.frames-1)*handles.parentHandles.metadata.meanSP;
handles.tright.String = num2str(tlim);
% Update handles structure
guidata(hObject, handles);

axes(handles.sourceMap)
try    % In here to intercept any errors from poorly used program!
    handles.imaging = imagesc(handles.parentHandles.mdImage);
    colormap('parula');
    set(gca,'fontsize',8)
end
guidata(hObject, handles);

axes(handles.outputMap)
set(gca,'fontsize',8)
set(gca,'YTick',[])

axes(handles.plotEph)
axis([0 tlim -1 1]);
handles.ref1 = line([1 1],[-1 1],'Color','r');
set(gca,'fontsize',8)
set(gca,'XTick',[])

axes(handles.RMSpower)
axis([0 tlim -1 1]);
handles.ref2 = line([1 1],[-1 1],'Color','r');
hold on;
set(gca,'fontsize',8)
set(gca,'XTick',[])

axes(handles.imagingData)
axis([0 tlim -1 1]);
handles.ref3 = line([1 1],[-1 1],'Color','r');
%hold on;
set(gca,'fontsize',8)
guidata(hObject, handles);

plotCaData (hObject, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes loadEPhys wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = loadEPhys_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function spTxt_Callback(hObject, eventdata, handles)

function spTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nChTxt_Callback(hObject, eventdata, handles)

function nChTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function workCh_Callback(hObject, eventdata, handles)

function workCh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function synchCh_Callback(hObject, eventdata, handles)

function synchCh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)

function edit5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function synchCK_Callback(hObject, eventdata, handles)

function nPoles_Callback(hObject, eventdata, handles)

function nPoles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fromHz_Callback(hObject, eventdata, handles)

function fromHz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function toHz_Callback(hObject, eventdata, handles)

function toHz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function openEPh_Callback(hObject, eventdata, handles)

display('Open Electrophysiology Data');
% Open a file holding electrophysiology data

if handles.AsciiRB.Value, fmtSt = '*.dat';
end
if handles.axopatchRB.Value, fmtSt = '*.abf';
end
if handles.meaRB.Value, fmtSt = '*.mat';
end
if handles.matlabRB.Value, fmtSt = '*.mat';
end
handles.dir_in = '';
[selection,path_in,filterI] = uigetfile ([handles.dir_in fmtSt]);
fileInfo = {'name' 'date' 'bytes' 'isdir' 'datenum'};
guidata(hObject, handles);
if selection ~= 0
    fileInfo = dir([path_in selection]);
    timeNum = datestr(fileInfo.datenum);
    handles.fileTime = timeNum (13:20);
    guidata(hObject, handles);
    openFile(hObject,handles,selection,path_in)
    handles = guidata(hObject);
end
display('End of Electrophysiology Data Input');

guidata(hObject, handles);

function openFile(hObject,handles,selection,path_in)

% this function opens the selected files using the code appropriate for the
% specific data type

handles.file_in = selection;
handles.dir_in = path_in;
handles.dtaComplete=[path_in selection];
% create an handle for the file opened in read mode
fid=fopen(handles.dtaComplete,'r');

% initialize all arrays
LFPin = [];
handles.time = [];
handles.LFP = [];
handles.bpLFP = [];
handles.BPpower = [];
handles.logBP = [];

handles.sp = str2double(handles.spTxt.String);
handles.nCh = str2double(handles.nChTxt.String);
handles.wCh = str2double(handles.workCh.String);
sCh = handles.synchCh.String;
guidata(hObject, handles);      % refresh the content of handle

% Now the processing pipelines differentiate depending on the input mode.

if handles.AsciiRB.Value, fmtSt = '*.dat';
    % This function opens a specified file in ASCII binary format.
    % This is the default mode. Data from the NI board or from the save button of Zebra.
    % The number of channels must be defined at priori only.
    
    % Time axis is included in the data exported from Zebra as ASCII. Time it is
    % not included in the data generated from the NI board. If time axis is
    % included the sampling period is obtained directly from the input data.
    nc = handles.nCh;
    % Create the string describing the data format
    fmtSt='';
    for ii=1:handles.nCh
       fmtSt=[fmtSt '%f'];
    end
    if handles.timeIncluded.Value,
        % add one column for the time axis
        fmtSt=[fmtSt '%f'];
        nc = nc + 1;
    end    
    LFPin = fscanf(fid,fmtSt,[nc,inf]);
    l = size(LFPin);
    handles.dtaLen = l(2);
    if handles.timeIncluded.Value,
        handles.spTxt.String = LFPin(1,2);      % display the sampling period 
        handles.sp = LFPin(1,2);                % set the sampling period
        handles.time = LFPin(1,:);
        for i=1:handles.nCh
            % now move to the left all channels to remove the column that
            % was occupied by the time
            LFPin(i,:) = LFPin(i+1,:);
        end
    else
        handles.time = (0:handles.sp:handles.sp*(handles.dtaLen-1));
    end
end

if handles.axopatchRB.Value,
    % in this modality data are generated by the Axopatch program and
    % are saved as binary data. sp and Ch number obtained from data header
    disp([num2str(toc) ' Begin loading of Axopatch data'])
    [LFPin,si,h] = abfload(handles.dtaComplete);
    % si    scalar           the sampling interval in micro s
    % h     struct           information on file (selected header parameters)
    % number of channels
    LFPin = LFPin';    % abfloat return the channels as columns. Traspose to set them as rows
    setGain (hObject, handles, 1000);
    handles = guidata(hObject);
    setSP (hObject, handles, 1/si*1000000);
    handles = guidata(hObject);
    setCh (hObject, handles, h.nADCNumChannels);
    handles = guidata(hObject);
end
handles.LFP = LFPin;
handles.displayLFP.Visible = 'on';

% plot the EPh
handles.tright.String = num2str(handles.sp * (handles.dtaLen-1));   % be precise, please!
plotEphys (hObject, handles)

guidata(hObject, handles);

%%

function processEPh_Callback(hObject, eventdata, handles)
% Compute BD data and RMS power

display('BP filtering and computation of RMS power');

handles.BPpwrTime = [];
handles.BPpower = [];
poles = str2double(handles.nPoles.String);
FB(1) = str2double(handles.fromHz.String);
FB(2) = str2double(handles.toHz.String);
[b,a] = butter(poles,2*FB*handles.sp);  % band pass filter
for i = 1:handles.nCh
    tmp1 = [];
    tmp2 = [];
    tmp1(1:handles.dtaLen) = handles.LFP(i,1:handles.dtaLen);
    tmp2 = filtfilt(b,a,tmp1.');
    handles.bpLFP(i,:)=tmp2(:);
end
guidata(hObject,handles); 

% compute the width of the window for power computation
winT = str2double(handles.RMSwidth.String);
pntWin = floor(winT/handles.sp);            % Width in sampling points  
ovl = str2double(handles.shiftWnd.String);
pntOvl = floor(ovl/handles.sp);
pntWin2 = floor(pntWin/2);                  % Half width in sampling points
centerBegin = floor(pntWin/2+1);

% compute the indexes of the ORIGINAL data to use for the power computation
BPindexes = centerBegin:pntOvl:(handles.dtaLen-pntWin2);  % this array contains the indexes of the 
                                                          % BP data around which are computed the RMS powers
xBPpowerEnd = length(BPindexes);                          % number of segments for power computation  
BPindLeft(1:xBPpowerEnd) = BPindexes(1:xBPpowerEnd)-pntWin2;    % left and right limits (in sampling points)
BPindRight(1:xBPpowerEnd) = BPindexes(1:xBPpowerEnd)+pntWin2;   % of the integration window
if (BPindLeft(1)<1), BPindLeft(1) = 1;
end
if (BPindRight(xBPpowerEnd)>handles.dtaLen), BPindRight(xBPpowerEnd)=handles.dtaLen;
end     % not certain about the last line...
% compute X axis
handles.BPpwrTime (1:xBPpowerEnd) = handles.sp * (BPindexes(1:xBPpowerEnd)-1);  % time coordinate
handles.BPpwrToffset = handles.BPpwrTime (1);
if (xBPpowerEnd>1),         % again, more controls to protect bad bandwidths
    handles.BPpowerSP = handles.sp * (BPindexes(2)-BPindexes(1));   % sampling period of the power. It should be equal to the Shift value
else    
    handles.BPpowerSP(BPsel) = 0;
end
% preallocate the power
handles.BPpower = zeros(handles.nCh, xBPpowerEnd);

% now compute the short time RMS power
% Lets cycle on the filtered data
for iCh = 1:handles.nCh
    % next loop computes the power, finally!
    for timI = 1:xBPpowerEnd
        handles.BPpower (iCh,timI) = rms...
           (handles.bpLFP(iCh,BPindLeft(timI):BPindRight(timI)));
    end
end
handles.logBP = log10(handles.BPpower); 
plotEphys (hObject, handles)
display('End of BP filtering and computation of RMS power');
handles.logPW.Visible = 'on';
handles.displayLFP.Visible = 'on';

guidata(hObject,handles);   

%%
function xcorr_Callback(hObject, eventdata, handles)
 
% compute the cross correlation map between the interpolated log of gamma power
% and the deltaF image set

handles.CCmap = [];

corrWin = str2double(handles.corrWinTxt.String);
intWin = str2double(handles.intWndTxt.String);
% these data must be trasformed in number of frames
corrWin = floor(corrWin/handles.parentHandles.metadata.meanSP);
intWin = floor(intWin/(2*handles.parentHandles.metadata.meanSP));
intFrom = corrWin-intWin;
intTo = corrWin+intWin;

handles.CCmap = computeXcorr(hObject, handles, corrWin, intFrom, intTo);
minEmap = min(min(handles.CCmap));
maxEmap = max(max(handles.CCmap));
% now normalize the map. At this stage the map is normalized between -1 and
% +1 with the following criteria
absMax = max([abs(minEmap) maxEmap]);
% the largest CC (be a positive or negative correlation) is mapped to 1
handles.CCmap = handles.CCmap/absMax;

handles.statsElectroMap.String = [' Min: ' num2str(minEmap/absMax, '%2.3g')...
    ' Max: ' num2str(maxEmap/absMax, '%2.3g')];
guidata(hObject, handles)

plotEmap(hObject, handles);
guidata(hObject, handles)

%% now compute limits of the map

function map = computeXcorr(hObject, handles, corrWin, intFrom, intTo)
% First, grab the time coordinate from the imaging set and interpolate the
% E-Phys data to the same vector

if handles.logPW.Value,     % use either the BP power or its log depending on selection
    ephData = handles.logBP (handles.wCh,:);
else
    ephData = handles.BPpower (handles.wCh,:);
end

vq = interp1(handles.BPpwrTime,ephData,handles.parentHandles.timeI,'linear','extrap');

% verify that there are no NaN data in the interpolation

% Next, the computation of the CC spectra requires that the input data have mean=0
% Subtract the mean to the interpolated EF data

meanEF = mean(vq);
vq = vq - meanEF;

% now loops on all pixels of the image stack to compute the CC matrix
tic
map = zeros(handles.parentHandles.WIrows,handles.parentHandles.WIcolumns);
for row = 1:handles.parentHandles.WIrows
    for col = 1:handles.parentHandles.WIcolumns
        % extract the array from the df image.
        dF = handles.parentHandles.deltaF(row,col,1,:);
        dF = permute(dF,[4 3 2 1]);
        [corr, lag] = xcorr(vq,dF,corrWin);
        % please integrate in the corr window. corr is  an array of corrWin+1 elements. Lag
        % 0 correspond to the central point of the array.
        % The XCorr metrics is computed by integrating in the intFrom:intTo window
        intCC = sum(corr(intFrom:intTo));
        % we also compute an assimmetry map that computed the Ca response
        % assimmetry with respect of the EF peak
        %intRightLobe = sum(corr(121:151));
        map (row,col) = intCC;
        %Delay (row,col) = intRightLobe;
    end
end
toc
%%

function plotEphys (hObject, handles)

% Prepare the default values for time serie display

% find time limits of the plotted data
tmin = str2double(handles.tleft.String);
tmax = str2double(handles.tright.String);
imin = tmin/handles.sp + 1;
imax = tmax/handles.sp + 1;
if imax > handles.dtaLen, imax = handles.dtaLen;
end
if imin <1, imin = 1;
end

axes(handles.plotEph)
% hide the handle of the reference line. In this way it will not be deleted
% with the subsequent plot callbacks
handles.ref1.HandleVisibility  = 'off';
if length(handles.LFP)>1 && handles.displayLFP.Value,
    plot(handles.time(imin:imax),handles.LFP(handles.wCh,imin:imax), 'Color','b');
    hold on;        % in case we get to see also the band passed data
end

if length(handles.bpLFP)>1 && handles.displayBP.Value,
    plot(handles.time(imin:imax),handles.bpLFP(handles.wCh,imin:imax), 'Color','g');
end
axis([tmin tmax -inf inf]);
set(gca,'fontsize',8)
set(gca,'XTick',[])
handles.ref1.HandleVisibility  = 'on';
hold off;
handles.plotEph.NextPlot = 'replacechildren'; 

% plot the RMS power
if length(handles.BPpower)>1
    axes(handles.RMSpower)
    % remember that there is a small time offset in the power
    imin = int32((tmin - handles.BPpwrToffset)/handles.BPpowerSP + 1);
    imax = int32((tmax - handles.BPpwrToffset)/handles.BPpowerSP + 1);
    
    % quality control of the index limits
    if imin<1, imin=1;      
    end
    if imax>length(handles.BPpower), imax=length(handles.BPpower);
    end
    
    handles.ref2.HandleVisibility  = 'off';     % hide the reference line
    if length(handles.BPpower)>1,
        if handles.logPW.Value,
            Ymin = min(handles.logBP(handles.wCh,imin:imax));
            Ymax = max(handles.logBP(handles.wCh,imin:imax));
            plot(handles.BPpwrTime(imin:imax),handles.logBP(handles.wCh,imin:imax))
        else
            Ymin = min(handles.BPpower(handles.wCh,imin:imax));
            Ymax = max(handles.BPpower(handles.wCh,imin:imax));
            plot(handles.BPpwrTime(imin:imax),handles.BPpower(handles.wCh,imin:imax))
        end    
    end
    
    axis([tmin tmax Ymin Ymax]);
    set(gca,'fontsize',8)
    set(gca,'XTick',[])
    handles.ref2.HandleVisibility  = 'on';
end

guidata (hObject, handles)
plotRefLine (hObject, handles)

guidata (hObject, handles)

function plotCaData (hObject, handles)

if length(handles.parentHandles.binarySum)>1
    axes(handles.imagingData)
    % Prepare the default values for time serie display

    % find time window
    tmin = str2double(handles.tleft.String);
    tmax = str2double(handles.tright.String);
    imin = int16(1 + tmin/handles.parentHandles.metadata.meanSP);
    imax = int16(1 + tmax/handles.parentHandles.metadata.meanSP);
    if imax > handles.parentHandles.metadata.frames,imax = handles.parentHandles.metadata.frames;
    end
    handles.ref3.HandleVisibility  = 'off';
    if handles.displaySelectedCa.Value,
        % display Ca time course of pixels with selected XCorr
        cmin = 1;
        cmax = max(handles.selBinarySum(imin:imax));
        plot (handles.parentHandles.timeI, handles.selBinarySum, 'Color','c','LineWidth',1)        
    else
        cmin = 1;
        cmax = max(handles.parentHandles.binarySum(imin:imax));
        plot (handles.parentHandles.timeI, handles.parentHandles.binarySum, 'Color','c','LineWidth',1)
    end
    handles.ref3.HandleVisibility  = 'on';
    if handles.logCa.Value
        set(gca, 'YScale', 'log')
    else
        set(gca, 'YScale', 'linear')
    end
    axis([tmin tmax cmin cmax]);
    set(gca,'fontsize',8)
end
guidata(hObject, handles)

function plotRefLine (hObject, handles)
% this function draw a line in correspondence of the selected frame in the 3
% time serie windows

% find time window for quality control
tmin = str2double(handles.tleft.String);
tmax = str2double(handles.tright.String);
lineTime = handles.parentHandles.metadata.meanSP * (handles.currentFrame-1);

imin = tmin/handles.sp + 1;
imax = tmax/handles.sp + 1;
if imax > handles.dtaLen,imax = handles.dtaLen;
    % compute the highest legal value for tmax
    tmax = (imax-1)*handles.sp;
end

% lineTime is the timing corresponding to the selected frame. If lineTime is 
% outside of the drawing area, the refeerence line is placed near the edge
% closest to the frame.
if (lineTime>tmax), lineTime = tmax;
end
if (lineTime<tmin), lineTime = tmin;
end

% ------------------------------
% Plot LFP frame reference line
axes(handles.plotEph)

% It is impossible to use the Ylim property to read the plot limits. It
% would return -Inf and +Inf which would not allow to plot the frame
% reference line. So, the range must be found for each window separately.
Ymin = 0;
Ymax = 0;
YminBP = 0;
YmaxBP = 0;

if handles.displayLFP.Value,
    Ymin = min(handles.LFP(handles.wCh,imin:imax));
    Ymax = max(handles.LFP(handles.wCh,imin:imax));
end
if handles.displayBP.Value,
    YminBP = min(handles.bpLFP(handles.wCh,imin:imax));
    YmaxBP = max(handles.bpLFP(handles.wCh,imin:imax));
end
Ymin = min([Ymin YminBP]);
Ymax = max([Ymax YmaxBP]);
if (Ymin==0 && Ymax==0),
    Ymin = -1;
    Ymax = 1;
end    
set(handles.ref1,'XData',[lineTime lineTime],'YData',[Ymin Ymax],'Color','m');

%-----------------------------------------
% Plot POWER frame reference line
if isfield(handles,'BPpwrToffset'),
    axes(handles.RMSpower)
    imin = int32((tmin - handles.BPpwrToffset)/handles.BPpowerSP + 1);
    imax = int32((tmax - handles.BPpwrToffset)/handles.BPpowerSP + 1);

    % quality control
    if imin<1, imin=1;      
    end
    if imax>length(handles.BPpower), imax=length(handles.BPpower);
    end
    if imax > length(handles.BPpower),imax = length(handles.BPpower);
        % compute the highest legal value for tmax
        tmax = (imax-1)*handles.sp;
    end
    if handles.logPW.Value,
        Ymin = min(handles.logBP(handles.wCh,imin:imax));
        Ymax = max(handles.logBP(handles.wCh,imin:imax));
    else
        Ymin = min(handles.BPpower(handles.wCh,imin:imax));
        Ymax = max(handles.BPpower(handles.wCh,imin:imax));
    end

    set(handles.ref2,'XData',[lineTime lineTime],'YData',[Ymin Ymax],'Color','m');
end

% -------------------------------
% Plot Ca frame reference line

axes(handles.imagingData)
imin = uint32(tmin/handles.parentHandles.metadata.meanSP + 1);
imax = uint32(floor(tmax/handles.parentHandles.metadata.meanSP + 1));
if imax > handles.parentHandles.metadata.frames,imax = handles.parentHandles.metadata.frames;
end
Ymin = min(handles.parentHandles.binarySum(imin:imax));
Ymax = max(handles.parentHandles.binarySum(imin:imax));

set(handles.ref3,'XData',[lineTime lineTime],'YData',[Ymin Ymax],'Color','m');

%%
function plotImages (hObject, handles)
    
% Now plot the mean image
axes(handles.sourceMap)
switch handles.sourceRB.SelectedObject.Tag
    case 'fluoRB'
        % read the axis limits
        Flim = [str2double(handles.FminTxt.String) str2double(handles.FmaxTxt.String)];
        if handles.maxPr.Value
            %handles.images.Cdata = handles.parentHandles.mdImage;
            imagesc(handles.parentHandles.mdImage, Flim);
        else
            % display single frame. Initial frame = 1
            imagesc(handles.parentHandles.workImages(:,:,:,handles.currentFrame), Flim);
        end
        
    case 'deltaFRB'
        % read the axis limits
        Flim = [str2double(handles.minDfTxt.String) str2double(handles.maxDfTxt.String)];
        if handles.maxPr.Value
            imagesc(handles.parentHandles.maxProjDF, Flim);
        else
            % display single frame. Initial frame = 1
            imagesc(handles.parentHandles.deltaF(:,:,:,handles.currentFrame), Flim);
        end        
    case 'binaryRB'
        if handles.maxPr.Value
            imagesc(handles.parentHandles.binaryProjection);
        else
            % display single frame. Initial frame = 1
            frame = squeeze(handles.parentHandles.erodedBinaryImages(:,:,1,handles.currentFrame));
            imagesc(frame);
        end        
end    
set(gca,'fontsize',8)

% try    % In here to intercept any errors from poorly used program!
%     imagesc(handles.parentHandles.mdImage)
%     colormap('parula');
%     set(gca,'fontsize',8)
% end

function plotEmap (hObject, handles)

% plot the correlation maps between Ephysiol and Ca imaging
axes(handles.outputMap)
% read the plot limits:
lim = [str2double(handles.minEmap.String) str2double(handles.maxEmap.String)];
if handles.displayThrCCmap.Value
    imagesc(handles.thresholdedCCmap, lim)
else
    imagesc(handles.CCmap, lim)
end
set(gca,'fontsize',8)
set(gca,'YTick',[])

function timeIncluded_Callback(hObject, eventdata, handles)

function RMSwidth_Callback(hObject, eventdata, handles)

function RMSwidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function shiftWnd_Callback(hObject, eventdata, handles)

function shiftWnd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%

function tleft_Callback(hObject, eventdata, handles)

plotRefLine (hObject, handles)
plotEphys (hObject, handles)
plotCaData (hObject, handles)

function tleft_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function tright_Callback(hObject, eventdata, handles)

plotRefLine (hObject, handles)
plotEphys (hObject, handles)
plotCaData (hObject, handles)

function tright_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function selectFrame_Callback(hObject, eventdata, handles)
handles.currentFrame = int16(hObject.Value);
handles.currentTime = handles.parentHandles.metadata.meanSP*(handles.currentFrame-1);
guidata(hObject, handles)
setTimeOrFrame (hObject,handles)
plotImages (hObject, handles)
plotRefLine (hObject, handles);

function selectFrame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function frame1_Callback(hObject, eventdata, handles)
handles.currentFrame = 1;
handles.currentTime = handles.parentHandles.metadata.meanSP*(handles.currentFrame-1);
setTimeOrFrame (hObject,handles)
handles.selectFrame.Value = handles.currentFrame;

guidata(hObject, handles);
plotRefLine (hObject, handles)
plotImages (hObject, handles)

function lastFrame_Callback(hObject, eventdata, handles)
handles.currentFrame = handles.parentHandles.metadata.frames;
handles.currentTime = handles.parentHandles.metadata.meanSP*(handles.currentFrame-1);
setTimeOrFrame (hObject,handles)
handles.selectFrame.Value = handles.currentFrame;

guidata(hObject, handles);
plotRefLine (hObject, handles);
plotImages (hObject, handles)

function back1F_Callback(hObject, eventdata, handles)
if handles.currentFrame == 1,
    handles.currentFrame = handles.parentHandles.metadata.frames;
else
    handles.currentFrame = handles.currentFrame - 1;
end
handles.currentTime = handles.parentHandles.metadata.meanSP*(handles.currentFrame-1);
setTimeOrFrame (hObject,handles)
handles.selectFrame.Value = handles.currentFrame;    % Align the slider

guidata(hObject, handles);
plotImages(hObject, handles)
plotRefLine (hObject, handles)

function setTimeOrFrame (hObject, handles)
if handles.frameORsec.Value
    handles.currentFtxt.String = num2str(handles.currentFrame,'%i');
else
    handles.currentFtxt.String = num2str(handles.currentTime,'f*.1');    
end

function adv1F_Callback(hObject, eventdata, handles)
if handles.currentFrame == handles.parentHandles.metadata.frames,
    handles.currentFrame = 1;
else
    handles.currentFrame = handles.currentFrame + 1;
end
handles.currentTime = handles.parentHandles.metadata.meanSP*(handles.currentFrame-1);
setTimeOrFrame (hObject,handles)
handles.selectFrame.Value = handles.currentFrame;       % align the slider
guidata(hObject, handles);
plotImages(hObject, handles)
plotRefLine (hObject, handles)

function maxPr_Callback(hObject, eventdata, handles)
plotImages (hObject, handles)

function displayLFP_Callback(hObject, eventdata, handles)
plotEphys (hObject, handles)

function displayBP_Callback(hObject, eventdata, handles)
plotEphys (hObject, handles)

function logPW_Callback(hObject, eventdata, handles)
plotEphys (hObject, handles)

function corrWinTxt_Callback(hObject, eventdata, handles)

function corrWinTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function intWndTxt_Callback(hObject, eventdata, handles)

function intWndTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minDfTxt_Callback(hObject, eventdata, handles)
plotImages (hObject, handles)

function minDfTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxDfTxt_Callback(hObject, eventdata, handles)
plotImages (hObject, handles)

function maxDfTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FminTxt_Callback(hObject, eventdata, handles)
plotImages (hObject, handles)

function FminTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FmaxTxt_Callback(hObject, eventdata, handles)
plotImages (hObject, handles)

function FmaxTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sourceRB_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in sourceRB 
plotImages (hObject, handles)

% --- Executes on button press in emapFreq.
function emapFreq_Callback(hObject, eventdata, handles)

function minEmap_Callback(hObject, eventdata, handles)
plotEmap (hObject, handles)

function minEmap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxEmap_Callback(hObject, eventdata, handles)
plotEmap (hObject, handles)

function maxEmap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fromCC_Callback(hObject, eventdata, handles)

function fromCC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function toCC_Callback(hObject, eventdata, handles)

function toCC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function extractCC_Callback(hObject, eventdata, handles)
% extract the pixels associated to a given value of the XCorr

if hObject.Value,
    CCmin = str2double(handles.fromCC.String);
    CCmax = str2double(handles.toCC.String);
    % the goodPoint array hold NaN at all coordinates that do not satisfy the
    % condition
    goodPoints = handles.CCmap;
    % set to NaN all points that do not satisfy the overall XCorr condition
    goodPoints(goodPoints>CCmax | handles.CCmap<CCmin) = NaN;
    handles.thresholdedCCmap = goodPoints;
    % set to 1 all points that do satisfy the condition
    goodPoints(~isnan(goodPoints)) = 1;
    for i=1:handles.parentHandles.metadata.frames
        frame = squeeze(handles.parentHandles.erodedBinaryImages(:,:,1,i));
        filteredFrame = goodPoints .* frame;
        numGood = sum(sum(filteredFrame(isfinite(filteredFrame))));
        handles.selBinarySum (i) = numGood; 
    end
end
guidata(hObject, handles)

function logCa_Callback(hObject, eventdata, handles)
plotCaData (hObject, handles)

function displaySelectedCa_Callback(hObject, eventdata, handles)
plotCaData (hObject, handles)

function displayThrCCmap_Callback(hObject, eventdata, handles)
plotEmap (hObject, handles)

function loadEvents_Callback(hObject, eventdata, handles)
% Open a ASCII files containing events identified by Zebra

fmtSt = '*.dat';

[selection,path_in,filterI] = uigetfile ([handles.dir_in fmtSt]);
fileInfo = {'name' 'date' 'bytes' 'isdir' 'datenum'};
eventDta = [path_in selection];
fid=fopen(eventDta,'r');
guidata(hObject, handles);
if selection ~= 0
    handles.events = [];        % clear old data
    fmtSt='%f %f %f %f';
    % skip the string descriptors
    tmp = fscanf(fid,'%s',[1,1]);
    tmp = fscanf(fid,'%s',[1,1]);
    tmp = fscanf(fid,'%s',[1,1]);
    tmp = fscanf(fid,'%s',[1,1]);
    tmp = fscanf(fid,'%s',[1,1]);
    handles.events = fscanf(fid,fmtSt,[4,inf]);
    handles.events = handles.events.';
    % field description n time-start time-end power 
end
[handles.nEvents, ~] = size(handles.events);
guidata(hObject, handles);
computeEventMaps (hObject, handles)

function computeEventMaps (hObject, handles)
% This function compute the correlation map for each event entered in the
% event file
baseline = str2double(handles.eventBaseline.String);
delay = str2double(handles.delayCaRecovery.String);
if handles.nEvents>0, handles.currentMap=1;
else
    handles.currentMap=0;
end    
oldCa2 = -1;    % used to control for collisions of baselines of successive events
oldCa3 = -1;

% Refresh the handles from the eventDetection GUI
handles.extractData = guidata(handles.extract);
for i=1:handles.nEvents
    % compute the indexes of the baseline and of the event in Ca time
    % series
    tbegin = handles.events(i,2);
    tend = handles.events(i,3);
    iCa1 = int16(1+(tbegin-baseline)/handles.parentHandles.metadata.meanSP);
    if iCa1 <1, iCa1=1;
    end    
    iCa2 = int16(1+tbegin/handles.parentHandles.metadata.meanSP);
    iCa31 = iCa2 + 1;
    iCa3 = int16(1+(tend + delay)/handles.parentHandles.metadata.meanSP);
    % now check for collision of the baseline with the previous event
    if iCa1 <= oldCa3,
        % collision!
        iCa1 = oldCa1;
        iCa2 = oldCa2;
    end
    % indexes for baseline and event
    oldCa1 = iCa1;
    oldCa2 = iCa2;
    oldCa3 = iCa3;
    % now compute the event maps
    tmp = squeeze(handles.parentHandles.deltaF(:,:,1,iCa1:iCa2));
    tmpBinBsl = squeeze(handles.parentHandles.deltaF(:,:,1,iCa1:iCa2) .* ...
        handles.parentHandles.erodedBinaryImages(:,:,1,iCa1:iCa2));
    tmpBinBsl(isnan(tmpBinBsl)) = 0;
    tmpMd1 = mean(tmp,3);
    tmp = squeeze(handles.parentHandles.deltaF(:,:,1,iCa31:iCa3));
    tmpBinEv = squeeze(handles.parentHandles.deltaF(:,:,1,iCa31:iCa3) .* ...
        handles.parentHandles.erodedBinaryImages(:,:,1,iCa31:iCa3));
    tmpBinEv(isnan(tmpBinEv)) = 0;

    tmpMd2 = mean(tmp,3);
    handles.eventMap(:,:,i) = tmpMd2-tmpMd1;
        
    % Now compute metrics from the binarized and filtered Ca sequence
    % 2 different metrics: 1) the total number of NEWLY recruited pixels minus the mean number of 
    % NOVEL pixels recruited during the baseline (and normalize for the
    % relative length). WRONG!!!!!
    
    % 2) the total number of pixel/event over the threshold during the event.
    
    bsl = squeeze(handles.parentHandles.erodedBinaryImages(:,:,1,iCa1:iCa2));
    evnt = squeeze(handles.parentHandles.erodedBinaryImages(:,:,1,iCa31:iCa3));
    
    caMeanBsl = mean(bsl,3);
    caCountBsl = sum(bsl,3);
    newPxBsl = length(find(caCountBsl));
    totalPxCountBsl = sum(sum(caCountBsl, 1));
    
    caMeanEvn = mean(evnt,3);
    caCountEvn = sum(evnt,3);
    newPxEvn = length(find(sum(evnt,3)));           % Pixels recruited during the event. Recruited area 
                                                    % > proportional to power
    totalPxCountEvt = sum(sum(caCountEvn, 1));      % Integral of the pixels recruited during the event 
                                                    % > proportional to energy 
    
    % compute total number of newly recruited pixels
    
    %sum
    handles.events(i,6) = newPxEvn; 
    handles.events(i,7) = totalPxCountEvt;
    handles.events(i,8) = totalPxCountBsl;
    % normalization of the baseline to the event duration. This factor has
    % to be computed on the frame number belonging to the baseline and event window.
    normForBaseline = (iCa31 - iCa3 + 1)/(iCa2 - iCa1 +1);
    handles.events(i,9) = totalPxCountEvt - (normForBaseline * totalPxCountBsl);
    handles.events(i,10) = sum(sum(sum(tmpBinEv,3),2)) - normForBaseline * sum(sum(sum(tmpBinBsl,3),2));
    % copy the data to the eventDetection table
    handles.extractData.USlist (i,7) = num2cell(newPxEvn);
    handles.extractData.USlist (i,8) = num2cell(handles.events(i,9));
    handles.extractData.USlist (i,9) = num2cell(handles.events(i,10));
end
set(handles.extractData.UStable,'Data',handles.extractData.USlist);
% save data in the GUI handles structure
%guidata(handles.extract,handles.extractData);
guidata(hObject, handles);

function delayCaRecovery_Callback(hObject, eventdata, handles)

function delayCaRecovery_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function eventBaseline_Callback(hObject, eventdata, handles)

function eventBaseline_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plotEventMap(hObject, handles)
minDf = str2double(handles.minDfTxt.String);
maxDf = str2double(handles.maxDfTxt.String);
axes(handles.outputMap)
imagesc(handles.eventMap(:,:,handles.currentMap),[minDf maxDf]);
set(gca,'fontsize',8)
set(gca,'YTick',[])

function backMap_Callback(hObject, eventdata, handles)
handles.currentMap = handles.currentMap-1;
if handles.currentMap==0, handles.currentMap=handles.nEvents;
end
handles.mapN.String = num2str(handles.currentMap);
plotEventMap (hObject, handles)
drawEvents (hObject, handles);
guidata(hObject,handles)

function advMap_Callback(hObject, eventdata, handles)
handles.currentMap = handles.currentMap+1;
if handles.currentMap>handles.nEvents, handles.currentMap=1;
end
handles.mapN.String = num2str(handles.currentMap);
plotEventMap (hObject, handles)
drawEvents (hObject, handles);
guidata(hObject,handles)

function drawEvents (hObject, handles)
%
if handles.nEvents>0,
    % not sure where to place the reference... try Ca=10
    axes (handles.imagingData)
    hold on
    for i=1:handles.nEvents
        if i==handles.currentMap, clr='r';
        else
            clr='g';
        end    
       axes (handles.imagingData)
       line([handles.events(i,2) handles.events(i,3)],[10 10],'Color',clr,'LineWidth',2)
    end
    hold off
end

function evMap_Callback(hObject, eventdata, handles)
if hObject.Value & handles.nEvents,
    plotEventMap (hObject, handles)
    drawEvents (hObject, handles);
end

function frameORsec_Callback(hObject, eventdata, handles)

function extractEvents_Callback(hObject, eventdata, handles)

handles = guidata(hObject);  % refresh handles. It should not be necessary, but better be safe than sorry.
if ~isfield(handles,'LFP'), errNow = errordlg('Load an Ephys data first','Monkey Error');
else
    if ~isfield(handles,'extract') % check if the child process already exists by checking the existance of the field. 
        % if it does not exist the call to Xcorr must start no matter of
        % what: open the child process
        handles.extract = eventDetectionMelies(gcf,handles);      % handle to the child figure
                                                     % arguments: parent object and
                                                     % parent handles
        handles.extractData = guidata(handles.extract);  % and its controls
        hObject.Value = 1;                           % make sure that the checkbox correctly reflects the child status
    else
        % two possibilites: it exists and then it must be closed. or it
        % exists but it has been deleted. In that case it must be restarted
        % somehow
        if ~isvalid(handles.extract) % now check whether is a valid handle.
            handles.extract = eventDetectionMelies(gcf,handles);      % handle to the child figure
                                                         % arguments: parent object and
                                                         % parent handles
            handles.extractData = guidata(handles.Xcorr);  % and its controls
            hObject.Value = 1;                           % make sure that the checkbox correctly reflects the child status        
        else
            % close the child process
            close(handles.extract);
            delete(handles.extract);      % delete the handle to the child process
            hObject.Value = 0;
        end
    end
end    

guidata(hObject, handles);

function analyzeCaTransient_Callback(hObject, eventdata, handles)
computeEventMaps (hObject, handles)

function xcorrMap_Callback(hObject, eventdata, handles)
if hObject.Value,
    plotEmap (hObject, handles)
end
