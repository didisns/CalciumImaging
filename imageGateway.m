function varargout = imageGateway(varargin)
% IMAGEGATEWAY MATLAB code for imageGateway.fig
% December 21, 2017. O&G
% January 2017.
% October 2018.
% Image Gateway has been edited to comply with MatLab 2018b vagaries.
% It includes code written by OC

%      IMAGEGATEWAY, by itself, creates a new IMAGEGATEWAY or raises the existing
%      singleton*.
%
%      H = IMAGEGATEWAY returns the handle to a new IMAGEGATEWAY or the handle to
%      the existing singleton*.
%
%      IMAGEGATEWAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEGATEWAY.M with the given input arguments.
%
%      IMAGEGATEWAY('Property','Value',...) creates a new IMAGEGATEWAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imageGateway_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imageGateway_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imageGateway

% Last Modified by GUIDE v2.5 06-Nov-2018 17:50:10

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @imageGateway_OpeningFcn, ...
                   'gui_OutputFcn',  @imageGateway_OutputFcn, ...
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


% --- Executes just before imageGateway is made visible.
function imageGateway_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imageGateway (see VARARGIN)

% Choose default command line output for imageGateway
handles.output = hObject;

% Add here all the initialization you need to run at program launch
handles.dirIn = [pwd '\'];

% initialize structure for AB-map limits
handles.ABmapLimits = struct('KSRB',[0 1],'KSstdRB',[0 1],'areaRB',[0 1],'areaStdRB',[0 1]);
handles.ABmapLimB = struct('KSRB',[0 1],'KSstdRB',[0 1],'areaRB',[0 1],'areaStdRB',[0 1],'RGBcodin',[0 0.2]);
handles.ABmap = struct('KSRB',[],'KSstdRB',[],'areaRB',[],'areaStdRB',[]);
% inactivate controls that are illegal at opening
offABnormalMap(hObject, handles);
handles.FmaskFlag = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imageGateway wait for user response (see UIRESUME)
% uiwait(handles.imageGateway);

function offABnormalMap(hObject, handles)

handles.KSstdRB.Visible = 'off';
handles.KSRB.Visible = 'off';
handles.KSstdRBb.Visible = 'off';
handles.KSRBb.Visible = 'off';
handles.areaRB.Visible = 'off';
handles.areaRBb.Visible = 'off';
handles.areaStdRB.Visible = 'off';
handles.areaStdRBb.Visible = 'off';

guidata(hObject, handles);

function onABnormalMap (hObject, handles)

handles.KSstdRB.Visible = 'on';
handles.KSRB.Visible = 'on';
handles.KSstdRBb.Visible = 'on';
handles.KSRBb.Visible = 'on';
handles.areaRB.Visible = 'on';
handles.areaRBb.Visible = 'on';
handles.areaStdRB.Visible = 'on';
handles.areaStdRBb.Visible = 'on';

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = imageGateway_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function openBTN_Callback(hObject, eventdata, handles)

% Select the proper mode opening code
offABnormalMap(hObject, handles);

if get(handles.prairieRB,'Value'), fmtSt = '*.tif';
end
if get(handles.TIFFmultiImageRB,'Value'), fmtSt = '*.tif';
end
if get(handles.aviRB,'Value'), fmtSt = '*.avi';
end

[selection,path_in,filterI] = uigetfile ([handles.dirIn fmtSt]);
% update the current data directory
handles.dirIn = path_in;
fileInfo = {'name' 'date' 'bytes' 'isdir' 'datenum'};    % structure for storage of file info
guidata(hObject, handles);
if selection ~= 0
    % set ot invisible all illegal controls
    
    fileInfo = dir([path_in selection]);
    timeNum = datestr(fileInfo.datenum);
    handles.fileTime = timeNum (13:20);
    guidata(hObject, handles);
    % time to actually open the file
    openFile(hObject,handles,selection,path_in)
    handles = guidata(hObject);
end
guidata(hObject, handles);

function openFile(hObject,handles,selection,path_in)
    
% Open a image file: first clean the info field
e = tic;
set (handles.stackSizeTxt, 'String', '');
set (handles.timeAndSpaceTxt, 'String', '');
set (handles.moreInfoTxt, 'String', '');
refresh (gcf);
if get(handles.prairieRB,'Value'), fmtSt = '*.tif';
    pathAndFileIn = [path_in selection];
    % open a Prairie data directory
    strOut = [];
    [handles.rawImages handles.metadata] = import_PrairieTimeSequence (pathAndFileIn);
    strOut = ['Size: ' num2str(handles.metadata.xdim) 'x' num2str(handles.metadata.ydim)];
    strOut = [strOut 'x' num2str(handles.metadata.frames)];
    set (handles.stackSizeTxt, 'String', strOut);
    set (handles.stackSizeTxt, 'Visible', 'on');  
    strOut = [];
    strOut = ['Mean SP (s): ' num2str(handles.metadata.meanSP, '%1.4g')];
    strOut = [strOut '    Pixel size (um): ' num2str(handles.metadata.dx_um, '%2.4g')];
    set (handles.timeAndSpaceTxt, 'String', strOut);
    set (handles.timeAndSpaceTxt, 'Visible', 'on');  
    strOut = [];
    strOut = ['Objective: ' strjoin(handles.metadata.objective) '  NA: ' num2str(handles.metadata.objNA, '%1.3g')];
    strOut = [strOut '  Zoom factor: ' num2str(handles.metadata.zoom, '%2.3g')];
    set (handles.moreInfoTxt, 'String', strOut);  
    set (handles.moreInfoTxt, 'Visible', 'on');  
end

% Open TIF multiimage
if get(handles.TIFFmultiImageRB,'Value'), fmtSt = '*.tif';
    pathAndFileIn = [path_in selection];
    handles.tifInfo = imfinfo(pathAndFileIn);
    nframe = length(handles.tifInfo);
    handles.metadata.frames = nframe;
    handles.metadata.xdim = handles.tifInfo(1).Width;
    handles.metadata.ydim = handles.tifInfo(1).Height;
    handles.metadata.tpoints = handles.tifInfo(1).StripOffsets;
    handles.rawImages=zeros( handles.metadata.xdim,handles.metadata.ydim,1,nframe);
    tmp=zeros(handles.metadata.xdim,handles.metadata.ydim,1);   
    for i=1:nframe
        %offsettif=handles.tifInfo(i).Offset;        
        tmp = imread(pathAndFileIn,'Index',i);
        %tmp=tmp-offsettif;
        handles.rawImages(:,:,:,i)=tmp;
    end    
    
    %handles.metadata.ch = 
end

% create the x axis

handles.timeI = (0:handles.metadata.meanSP:handles.metadata.meanSP*(handles.metadata.frames-1));

% create the workImages file that hold the images as they move along the
% processing pipeline. workImage holds only the target channel
if get(handles.ch1RB,'Value'), handles.workImages = handles.rawImages(:,:,1,:);
    handles.metadata.ch = 1;
end
if get(handles.ch2RB,'Value'), handles.workImages = handles.rawImages(:,:,2,:);
    handles.metadata.ch = 1;
end
if get(handles.Ch1Ch2RB,'Value'), handles.workImages = handles.rawImages;
    handles.metadata.ch = 2;
end
% Processing is performed on channel 1 unless there are 2 channels.
% Eventually it will be possible to process either ch 1 or 2 in case of a ratiometric measurement.
% if get(handles.TIFFmultiImageRB,'Value'), fmtSt = '*.tif';
%     handles.workImages = handles.rawImages;
% end
handles.processCh = 1;
guidata(hObject, handles);

% set all initial values
setAllInitialValues (hObject, handles);
e = toc(e);
disp(['End of load file: ' num2str(e) ' sec'])
handles = guidata(hObject);

function setAllInitialValues(hObject,handles)

% Set the empty data point to NaN
% This happens during spiral scan: the points outside of the spiral circle
% are set to zero, but they should not contribute to any of the statistics.
handles.workImages = cast(handles.workImages,'double');     % type recast necessary since the 
if get(handles.prairieRB,'Value')
handles.workImages(handles.workImages==0) = NaN;
end
    % original class of the data is uint16 that does not support NaN
offset = str2double(get(handles.offsetTxt,'String'));
handles.workImages = handles.workImages - offset;

handles.currentFrame = 1;
set(handles.selectFrame,'Max',handles.metadata.frames);
handles.deltaFflag = 0;
handles.statsFlag = 0;
handles.binaryFlag = 0;
handles.FmaskFlag = 0;
handles.FstatsFlag = 0;
handles.DFstatsFlag = 0;

guidata(hObject, handles);

function kernelTxt_Callback(hObject, eventdata, handles)

function kernelTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filterBtn_Callback(hObject, eventdata, handles)
% compute the low pass filtered image set

fltSize = str2double(get(handles.kernelTxt,'string'));
fltMatrix = [];
fltMatrix = ones(fltSize,fltSize) / fltSize^2;    % normalization
handles.lpImages = imfilter(handles.workImages(:,:,handles.processCh,:),fltMatrix);

if get(handles.SpLPfilterCk,'Value')
    % the filter is included in the pipeline analysis
    handles.workImages = handles.lpImages;
else
    % do nothing... pretty sure this will change.
end

guidata(hObject, handles);

function binTxt_Callback(hObject, eventdata, handles)

function binTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function computeBinning_Callback(hObject, eventdata, handles)
% This function compute the binned down image stack
% 
% This function reduces the resolution of the image by the given factor by
% creating larger pixels. The image size must be evenly disible by the bin size 

tic
display ('Compute binning');
handles.binnedImages = [];
if get(handles.binRB,'Value')
    handles.binF = str2double(get(handles.binTxt,'String'));
    % compute dimension of the binned images
    
    % how large is one image?
    % reserve space for the binned array
    handles.binnedImages = ones(handles.metadata.xdim/handles.binF, handles.metadata.ydim/handles.binF,...
        handles.metadata.ch, handles.metadata.frames, 'double');
    % Operate on one frame at the time
    % It processes both channels if they are both active, otherwise process
    % only the singleton channel of the workImages array
    for chJ = 1:handles.metadata.ch
        for i=1:handles.metadata.frames
           % extract the current frame
           cF = handles.workImages(:,:,chJ,i);
           tmp = binImage(cF,handles.binF); 
           handles.binnedImages(:,:,chJ,i) = tmp(:,:);
        end
    end
    handles.workImages = handles.binnedImages;
else
    % do nothing... pretty sure this will change.
end

[handles.WIrows, handles.WIcolumns, ~] = size(handles.workImages);

axes(handles.imageExplore);
cla;
colormap('hot');
imagesc(mean(handles.workImages,4),'ButtonDownFcn',@mdImageClick);
set(gca,'fontsize',8)

[handles.wy, handles.wx, f] = size(handles.workImages);
e = toc;
display (['End of binning: ' num2str(e) ' sec']);

guidata(hObject, handles);

function computeDf_Callback(hObject, eventdata, handles)
% Compute the DF/F matrix for the source data
% This can be computed according to three algorithms:
% Simplified: DF/F = (F(t)-mF)/mF  where mF indicates the median of F
% Adaptive: computed according to konnerth algorithm
% Bleach correction: an exponential decline is fitted to the a model of the
% data baseline and it is used to correct F before of computation of DF/F

% First: preallocate memory space. Keep in mind that the size of the stack
% changes according to the operation of the upstream pipeline.
% DF/F ic computed on 1 channel only

tic
display ('Computation of DF/F');

[xw yw ch f] = size(handles.workImages);
handles.deltaF = zeros(xw,yw,1,f,'double');

temp = handles.workImages(:,:,handles.processCh,:);

% Compute and display stats F-offset
handles.minDF = min(min(min(temp)));
handles.maxDF = max(max(max(temp)));
strOut = ['Min F: ' num2str(handles.minDF,'%2.4g')];
strOut = [strOut ' max: ' num2str(handles.maxDF,'%2.4g')];
set(handles.FstatsTxt,'String',strOut);

% third: compute DF/F according to the selected algorithm
if get(handles.simpleDfcK,'Value'),
   md = median (temp,4);
   handles.mdImage = md;
   handles.maxF = max(temp,[],4);
   
   for i=1:f
       handles.deltaF(:,:,:,i) = (temp(:,:,:,i)-md(:,:,:))./md(:,:,:);      % compute element by element division
   end
end

handles.minDF = min(min(min(handles.deltaF(isfinite(handles.deltaF)))));
handles.maxDF = max(max(max(handles.deltaF(isfinite(handles.deltaF)))));
strOut = ['Min DF/F: ' num2str(handles.minDF,'%2.4g')];
strOut = [strOut ' max: ' num2str(handles.maxDF,'%2.4g')];
set(handles.minMaxDfTxt,'String',strOut);

% compute the max projection and the average projection
handles.maxProjDF = max(handles.deltaF,[],4); 
handles.medianDF = median(handles.deltaF,4);
handles.deltaFflag = 1;

e = toc;
display (['End of DF/F computation: ' num2str(e) ' sec']);
guidata(hObject, handles);

function displayImageExplore (hObject, handles)
% February 27, 2018. Modified for alpha masking.
% April 14, 2018. Modified to display of the AB normal maps

Dflim = [str2double(get(handles.DFminTxt,'String')) ...
    str2double(get(handles.DFmaxTxt,'String'))];
Flim = [str2double(get(handles.FminTxt,'String')) ...
    str2double(get(handles.FmaxTxt,'String'))];

axes(handles.imageExplore);
cla;
%handles.mdImageW = imagesc(md,'ButtonDownFcn',@mdImageClick);

if get(handles.framesRB,'Value')
    % plot the current frame
    if get(handles.fluoRB,'Value')
        % plot fluorescence
        tmp = handles.workImages(:,:,:,handles.currentFrame);
       
        handles.mdImageW = imagesc(tmp,'ButtonDownFcn',@mdImageClick);
        set(handles.imageExplore,'Clim',Flim);  % impose the new axes limits
    else
        if handles.deltaFflag,      % check if aready computed
            % plot deltaF/F
            tmp = handles.deltaF(:,:,:,handles.currentFrame);
            handles.mdImageW = imagesc(tmp,'ButtonDownFcn',@mdImageClick);
            set(handles.imageExplore,'Clim',Dflim);
        end
    end
end

if get(handles.maxRB,'Value')
    % plot maximum projection
    if get(handles.fluoRB,'Value')
        % plot fluorescence
        handles.mdImageW = imagesc(handles.maxF,'ButtonDownFcn',@mdImageClick);
       
        set(handles.imageExplore,'Clim',Flim);
    else
        if handles.deltaFflag,      % check if aready computed
            % plot deltaF/F
            handles.mdImageW = imagesc(handles.maxProjDF,'ButtonDownFcn',@mdImageClick);
           
            set(handles.imageExplore,'Clim',Dflim);
        end
    end
end

if get(handles.averageRB,'Value')
    % plot maximum projection
    if get(handles.fluoRB,'Value')
        % plot fluorescence
        handles.mdImageW = imagesc(handles.mdImage,'ButtonDownFcn',@mdImageClick);
     
        set(handles.imageExplore,'Clim',Flim);
    else
        if handles.deltaFflag,      % check if aready computed
            % plot deltaF/F
            handles.mdImageW = imagesc(handles.medianDF,'ButtonDownFcn',@mdImageClick);
           
            set(handles.imageExplore,'Clim',Dflim);
        end
    end
end
% check if the mask preview flag is on
if handles.FmaskFlag,
    hold on
    h = imshow(handles.black); % Save the handle; we'll need it later
    hold off
    set(h, 'AlphaData', handles.Fmask);
end
% LUT and plot limits must be set at the end of drawing operations
if get(handles.fluoRB,'Value')
    set(handles.imageExplore,'Clim',Flim);
else 
    set(handles.imageExplore,'Clim',Dflim);
end
colormap('hot');

% Display the ABmaps
wndSel = handles.imageFrameSource.SelectedObject.Tag;
if isempty(strfind('framesRBmaxRBaverageRB',wndSel))
    % if the selected tag is not frames, max, or average, AB-map is selected
    %handles.mdImageW = imagesc(handles.ABmap.(wndSel),'ButtonDownFcn',@mdImageClick);
    imagesc(handles.ABmap.(wndSel));
    set(handles.imageExplore,'Clim',handles.ABmapLimits.(wndSel));
    colormap('parula');
end

set(gca,'fontsize',8)
guidata(hObject, handles);

function displayBinaryImage (hObject, handles)

axes(handles.binaryExplore);
if handles.binaryFlag,
    if handles.RGBcoding.Value
        % display the time-RGB encoded projection
        % use the ABmapLimB data to normalize the RGB image
        tmp = (handles.meanRGB-handles.ABmapLimB.RGBcodin(1));
        tmp = tmp/(handles.ABmapLimB.RGBcodin(2)-handles.ABmapLimB.RGBcodin(1));
        imagesc(tmp);
        set(handles.binaryExplore,'Clim',handles.ABmapLimB.RGBcodin);
    end
    if get(handles.FbinRB,'Value')
        % plot the current frame
        if get(handles.binaryRawRB,'Value')
            % plot raw binary frame
    %         handles.binImage = imagesc(tmp,'ButtonDownFcn',@mdImageClick);
    %         set(handles.imageExplore,'Clim',Flim);  % impose the new axes limits
            frame = squeeze(handles.binaryThrImages(:,:,1,handles.currentFrame));                     
            [row, column] = find(frame);
            scatter(column,row,1);
            axis([0 handles.WIrows 0 handles.WIcolumns]);
            set(gca,'Ydir','reverse')
        else
            frame = squeeze(handles.erodedBinaryImages(:,:,1,handles.currentFrame));
            [row, column] = find(frame);
            scatter(column,row,1);
            axis([0 handles.WIrows 0 handles.WIcolumns]);
            set(gca,'Ydir','reverse')
        end
    end
    if get(handles.binMaxRB,'Value')
        imagesc(handles.binaryProjection);
    end
end

% Display the ABmaps
wndSel = handles.binaryImageSource.SelectedObject.Tag;
if isempty(strfind('FbinRBbinMaxRB',wndSel))&& ~handles.RGBcoding.Value
    % AB-map selected
    imagesc(handles.ABmap.(wndSel(1:end-1)));
    set(handles.binaryExplore,'Clim',handles.ABmapLimB.(wndSel(1:end-1)));
    colormap('parula');
end

set(gca,'fontsize',8)
guidata(hObject, handles);

function mdImageClick(src,~)
% This function responds to the click event on the image displayed on the
% imageExplorer axes.
% src points to the graphic object that generated the call which is the
% image displayed in the window

[x,y,button] = ginput(1)           %  coordinates of the point in the axes units

% next: identify the size of the source image for quality control of the
% selected point

parentHandles = guidata(gcf);       % handles of the imageGateway figure
parentHandles.currentX = int16(x);  % let the parent figure know what we do here  
parentHandles.currentY = int16(y);  
fieldX = get(src,'xData');          % coordinates of the src image
fieldY = get(src,'yData');

% quality control
if inpolygon(x,y,fieldX,fieldY)
    strOut = ['x: ' num2str(parentHandles.currentX) ' y: ' num2str(parentHandles.currentY)];
    mf = parentHandles.mdImage(parentHandles.currentY,parentHandles.currentX);
    strOut = [strOut ' F: ' num2str(mf, '%6.3g')];
    set(parentHandles.pointerInfoTxt,'String',strOut);
    temp = [];
    % select the data to be plotted according to the radio button selection in
    % the timeExplore axes and to the active channel.
    
    if get(parentHandles.plotFcK,'Value')
        % This plot the raw fluorescence (not corrected for offset)
        % correct coordinates for the bin factor
        bf = str2double(get(parentHandles.binTxt,'String'));
        xx = bf * parentHandles.currentX;
        yy = bf * parentHandles.currentY;
        temp(1:parentHandles.metadata.frames) = parentHandles.rawImages(xx...
            ,yy,1,1:parentHandles.metadata.frames); 
    end
    if get(parentHandles.plotBinnedFcK,'Value')
        temp(1:parentHandles.metadata.frames) = parentHandles.binnedImages(parentHandles.currentY...
            ,parentHandles.currentX,1,1:parentHandles.metadata.frames); 
    end
    if get(parentHandles.plotDeltaFcK,'Value')
        temp(1:parentHandles.metadata.frames) = parentHandles.deltaF(parentHandles.currentY...
            ,parentHandles.currentX,1,1:parentHandles.metadata.frames); 
    end
    axes(parentHandles.timeExplore);
    plot (temp);
end

guidata(gco, parentHandles);

function imageExplore_ButtonDownFcn(hObject, eventdata, handles)

[x,y] = ginput(1);

function computeAll_Callback(hObject, eventdata, handles)

function SpLPfilterCk_Callback(hObject, eventdata, handles)

function computeBinaryStack_Callback(hObject, eventdata, handles)

tic
disp('Begin compute binary stack')
% compute the distribution of all pixel's DF/F
dSigma = 2;

% put all the DF data into a single vector to compute its global stats

computeDFstats_Callback(hObject, eventdata, handles)
handles = guidata(hObject);

newMd = handles.medianDF;
newStd = handles.stdDF;

allDF = handles.deltaF(:);
temp = allDF;
%newMd = median(temp(~isnan(temp)));       % stats of the entire sample
%handles.medianDF = newMd;
%newStd = std(temp(~isnan(temp)));
handles.stdDF = newStd;

oldMd = 10;                  % ficticious value used to define the variable
newtemp = [];

% the following procedure refine the median value by an interative
% procedure;
% 1) compute the median and the STD
% 2) take off all the data that are +dSigma away from the median
% 3) compute again the median
% 4) reiterate until the median stabilizes

while (abs(newMd-oldMd)>0.04*newStd)
    oldMd = newMd;
    oldStd = newStd;
    % extract good elements by removing outliers
    newtemp = temp(temp <= (oldMd+dSigma*oldStd));
    newMd = median(newtemp(isfinite(newtemp)));
    newStd = std(newtemp(isfinite(newtemp)));
    temp = newtemp;
end

handles.globalMd = newMd;
handles.meanMainMode = mean(temp(isfinite(temp)));
handles.stdMainMode = newStd;
handles.numMainMode = length(temp);
handles.numTotal = length(allDF);
guidata(hObject, handles);
disp('End of statistics computation')

% Score each pixel according to the DF metrics
% pixels are scored with the following metrics:
% 1) number of frames falling outside of the normality threshold
% 2) total power of the frames outside of the normality threshold
% 3) entropy of the sequence
% 4) number of events lasting longer than n frames

handles.nstd = str2double(get(handles.thrStdTxt,'String'));
handles.minLen = str2double(get(handles.minDurationTxt,'String'));
structEl = 1:handles.minLen+2;
structEl = structEl ./ structEl;
% now zero the edges of the structuring element
structEl(1) = 0;
structEl(handles.minLen+2) = 0;

handles.minSize = str2double(get(handles.minSizeTxt,'String'));
handles.minsizefiltering=str2double(get(handles.minsizefilt,'String'));
handles.maxsizefiltering=str2double(get(handles.maxsizefilt,'String'));
handles.rangefilt=[handles.minsizefiltering handles.maxsizefiltering];
%handles.rangefilt=handles.rangefilt';
structElxy = ones(handles.minSize);

thrIs = handles.meanMainMode + handles.nstd*handles.stdMainMode;
% the image stack binaryThrImage contains a binary representation of the DF stack with 1 in the 
% pixels where DF/F is larger than the 2*std limit
handles.binaryThrImages = handles.deltaF>thrIs;
% now it is time to remove pixels that are filtered away because of low S/N ratio
if get(handles.thresholdFck,'Value'),
    for i=1:handles.metadata.frames
        tmp = handles.binaryThrImages (:,:,1,i); 
        tmp = tmp .* handles.binMask;
        handles.binaryThrImages (:,:,1,i) = tmp(:,:);
    end
end
e = toc;
disp(['End of computation of the binary stack. Begin of filtering: ' num2str(e) ' sec'])

[r c f] = size(handles.deltaF);

% now remove all events shorter than the minimal duration
for sheetN=1:r
    % extract one sheet for the morphological erosion
    tmp = handles.binaryThrImages(sheetN,:,1,:);
    sheet = squeeze(tmp);
    erodedSheet = imerode(sheet,structEl);
    % and now apply the complementary dilation
    erodedSheet = imdilate(erodedSheet,structEl);
    handles.erodedBinaryImages(sheetN,1:c,1,1:f) = erodedSheet(1:c,1:f);
end

% now remove events smaller than a specified size
for i=1:handles.metadata.frames
    tmp = handles.erodedBinaryImages (:,:,1,i); 
    filteredFrame = imerode(tmp,structElxy);
    % 20-12-18 by Didi: added the imdilate function back in, with the
    % correct figure to dilate: the previously eroded binary stack (before
    % the figure input was 'temp' which therefore effectively undid all
    % filtering done by the imerode function
    filteredFrame = imdilate(filteredFrame,structElxy);
    %filteredFrame = bwareafilt(tmp,handles.rangefilt);
    handles.erodedBinaryImages (:,:,1,i) = filteredFrame;
end

e = toc;
disp(['End of computation of the binary stack: ' num2str(e) ' sec'])

% compute the projection of all images
handles.binaryProjection = sum(handles.erodedBinaryImages,4);

% compute the RGB coded binary time series. A standard color map has only
% 64 RGB entries. We will expand the LUT later on. The jet colormap plot
% the rainbow starting from blue to green, yellow and red. Time axis starts
% at the blue end of the spectra.

clMap = colormap('jet');        % clMap contains the RGB triplets
% for j=1:63
%     newj = 1 + 10 * (j-1)
%     newCmap (newj:newj+10,1) = linspace(clmap(j,1),clmap(j+1,1),10);         % R component
%     newCmap (newj:newj+10,2) = linspace(clmap(j,2),clmap(j+1,2),10);         % G component
%     newCmap (newj:newj+10,3) = linspace(clmap(j,3),clmap(j+1,3),10);         % B component
% end
colN = length(clMap);   % dimension of the colormap

% initialize RGB stack
handles.binRGB = zeros(r,c,3,handles.metadata.frames);
for i=1:handles.metadata.frames
    % compute the RGB code belonging to the time slice 
    tim = int16(1+(colN-1)*i/handles.metadata.frames);
    handles.binRGB (:,:,1,i) = clMap(tim,1) * handles.erodedBinaryImages (:,:,1,i);
    handles.binRGB (:,:,2,i) = clMap(tim,2) * handles.erodedBinaryImages (:,:,1,i);
    handles.binRGB (:,:,3,i) = clMap(tim,3) * handles.erodedBinaryImages (:,:,1,i);
    % for proper display of RGB images, the RGB chanels must be the third
    % dimension of the array
end
% project the RGB stack
handles.meanRGB = mean(handles.binRGB,4);

% compute the projection x, y of each frame
% this vector contains the time course of the number of px over threshold
%Emanuele 19-10-2018 aggiunto '.' al plot
handles.binarySum = squeeze(sum(sum(handles.erodedBinaryImages,2),1));
% plot the temporal profile of the number of 'meaningful' pixels
axes(handles.timeExplore);
plot(handles.binarySum,'.');
set(gca,'fontsize',8)

handles.binaryFlag = 1;
displayBinaryImage (hObject, handles)

guidata(hObject, handles);

function plotDistribution(hObject,handles)
% Ricompute and plot the distribution of the elected entity
% The recomputation allows to change the bin size on the fly

axes(handles.distributionW);
cla;
if (get(handles.binnedFck,'Value') && handles.FstatsFlag),
    bar(handles.Fcent,handles.Fnel);
    axis([handles.Fleft handles.Fright -inf inf])
end
if (get(handles.dFck,'Value') && handles.DFstatsFlag),
    bar(handles.DFcent,handles.DFnel);
    axis([handles.DFleft handles.DFright -inf inf])
end

% finally check whether we should apply a log scale 
if (get(handles.distributionLogCk,'Value')), set(handles.distributionW, 'YScale', 'log')
end
set(gca,'fontsize',8)

function minDistTxt_Callback(hObject, eventdata, handles)

function minDistTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxDistTxt_Callback(hObject, eventdata, handles)

function maxDistTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function binDistTxt_Callback(hObject, eventdata, handles)

function binDistTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function distributionLogCk_Callback(hObject, eventdata, handles)
plotDistribution(hObject,handles);

function thrStdTxt_Callback(hObject, eventdata, handles)

function thrStdTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minDurationTxt_Callback(hObject, eventdata, handles)

function minDurationTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function backMultiFrames_Callback(hObject, eventdata, handles)
if handles.currentFrame <= 10,
    handles.currentFrame = handles.metadata.frames;
else
    handles.currentFrame = handles.currentFrame - 10;
end
set(handles.currentFrTxt,'String',handles.currentFrame);
guidata(hObject, handles);
displayImageExplore (hObject, handles);
displayBinaryImage (hObject, handles);

function back1Frame_Callback(hObject, eventdata, handles)
if handles.currentFrame == 1,
    handles.currentFrame = handles.metadata.frames;
else
    handles.currentFrame = handles.currentFrame - 1;
end
set(handles.currentFrTxt,'String',handles.currentFrame);
guidata(hObject, handles);
displayImageExplore (hObject, handles);
displayBinaryImage (hObject, handles);

function advance1frame_Callback(hObject, eventdata, handles)
if handles.currentFrame == handles.metadata.frames,
    handles.currentFrame = 1;
else
    handles.currentFrame = handles.currentFrame + 1;
end
set(handles.currentFrTxt,'String',handles.currentFrame);
guidata(hObject, handles);
displayImageExplore (hObject, handles);
displayBinaryImage (hObject, handles);

function advanceMultiFrames_Callback(hObject, eventdata, handles)
if handles.currentFrame >= (handles.metadata.frames-9),
    handles.currentFrame = 1;
else
    handles.currentFrame = handles.currentFrame + 10;
end
set(handles.currentFrTxt,'String',handles.currentFrame);
guidata(hObject, handles);
displayImageExplore (hObject, handles);
displayBinaryImage (hObject, handles);

function globalStatsCk_Callback(hObject, eventdata, handles)

function offsetTxt_Callback(hObject, eventdata, handles)

function offsetTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function resetProcessing_Callback(hObject, eventdata, handles)
% reset the processing flow by reloading the rawImages to the workImage array
if get(handles.ch1RB,'Value'), handles.workImages = handles.rawImages(:,:,1,:);
    handles.metadata.ch = 1;
end
if get(handles.ch2RB,'Value'), handles.workImages = handles.rawImages(:,:,2,:);
    handles.metadata.ch = 1;
end
if get(handles.Ch1Ch2RB,'Value'), handles.workImages = handles.rawImages;
    handles.metadata.ch = 2;
end

% Set the empty data point to NaN
% This happens during spiral scan: the points outside of the spiral circle
% are set to zero, but they should not contribute to any of the statistics.
handles.workImages = cast(handles.workImages,'double');     % type recast necessary since the 
handles.workImages(handles.workImages==0) = NaN;            % original class of the data is uint16 that does not support NaN
offset = str2double(get(handles.offsetTxt,'String'));
handles.workImages = handles.workImages - offset;

guidata(hObject, handles);

function selectFrame_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.currentFrame = uint16(get(hObject,'Value'));
set(handles.currentFrTxt,'String',handles.currentFrame);
% display the frame pointed by handles.currentFrame
guidata(hObject, handles);

displayImageExplore (hObject, handles);
displayBinaryImage (hObject, handles);

function selectFrame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function animateFCK_Callback(hObject, eventdata, handles)

if hObject.Value,       % Activate time lapse viewer
    if get(handles.framesRB,'Value')
        % play the binned frames
        handles.Fmovie = implay(handles.workImages(:,:,1,:));
        handles.Fmovie.Visual.setPropertyValue('UseDataRange',true);
        handles.Fmovie.Visual.setPropertyValue('DataRangeMin',-Inf);
        handles.Fmovie.Visual.setPropertyValue('DataRangeMax',Inf);
    else
        handles.Fmovie = implay(handles.deltaF(:,:,1,:));        
    end
else    
    % Close already open animation window
    if ishandle(handles.Fmovie)
        close(handles.Fmovie);
    end    
end    
guidata(hObject, handles)

function DFminTxt_Callback(hObject, eventdata, handles)
displayImageExplore (hObject, handles)

function DFminTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function DFmaxTxt_Callback(hObject, eventdata, handles)
displayImageExplore (hObject, handles)

function DFmaxTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FminTxt_Callback(hObject, eventdata, handles)
displayImageExplore (hObject, handles)

function FminTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FmaxTxt_Callback(hObject, eventdata, handles)
displayImageExplore (hObject, handles)

function FmaxTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes when selected object is changed in imageFrame.
function imageFrame_SelectionChangedFcn(hObject, eventdata, handles)
displayImageExplore (hObject, handles)

function imageFrameSource_SelectionChangedFcn(hObject, eventdata, handles)
% set the scale limits in the 'Display Control' panel if a AB-map is being displayied
wndSel = handles.imageFrameSource.SelectedObject.Tag;
if isempty(strfind('framesRBmaxRBaverageRB', wndSel))
   % AB-map selected
   mn = handles.ABmapLimits.(wndSel)(1);
   mx = handles.ABmapLimits.(wndSel)(2);
   handles.ABmin.String = num2str(mn);
   handles.ABmax.String = num2str(mx);
end 
guidata(hObject, handles)
displayImageExplore (hObject, handles)

function createAVIbtn_Callback(hObject, eventdata, handles)
% Feb. 18th, 2018. 
% This function create a non compressed AVI movie of the binarized data

% create the output object
namevideo = [handles.dirIn 'BinaryVideo.avi'];
AVIout = VideoWriter (namevideo,'Uncompressed AVI');
open(AVIout);

% prepare the figure for the video out.
% remember to set the dimension according with the image stack size.
hFig = figure('units','pixels','position',[200 200 256 256]);
hMovie = axes('units','pixels','position',[0 0 256 256],'Visible','off');
axis ([0 handles.WIrows 0 handles.WIcolumns]);
set(hMovie, 'XTick', [], 'YTick', []);

%Emanuele 19-10-2018 cambiatto scatter con plot'.'
axes(hMovie);
for i=1:handles.metadata.frames   
    frame = squeeze(handles.erodedBinaryImages(:,:,1,i));    
    [row, column] = find(frame);    
    scatter(column,row,1);       
    axis ([0 handles.WIrows 0 handles.WIcolumns]);
    set(gca,'Ydir','reverse');
    %set(hMovie, 'XTick', [], 'YTick', []);
    refresh;
%    myavi = addframe(myavi, hFig);
    im = frame2im(getframe(hFig));
    im = rgb2gray(im);
    writeVideo(AVIout,im);
end
close(AVIout);

function uibuttongroup7_SelectionChangedFcn(hObject, eventdata, handles)
displayBinaryImage (hObject, handles)

function binaryImageSource_SelectionChangedFcn(hObject, eventdata, handles)
wndSel = handles.binaryImageSource.SelectedObject.Tag;
if isempty(strfind('FbinRBbinMaxRB',wndSel))
   % AB-map selected
   mn = handles.ABmapLimB.(wndSel(1:end-1))(1);
   mx = handles.ABmapLimB.(wndSel(1:end-1))(2);
   handles.ABminb.String = num2str(mn);
   handles.ABmaxb.String = num2str(mx);
end 
guidata(hObject, handles)

displayBinaryImage (hObject, handles)

function animateBin_Callback(hObject, eventdata, handles)

if hObject.Value,       % Activate time lapse viewer
    if get(handles.binaryRawRB,'Value')
        % play the raw (no binary flters) binned frames 
        handles.Bmovie = implay(handles.workImages(:,:,1,:));
        handles.Bmovie.Visual.setPropertyValue('UseDataRange',true);
        handles.Bmovie.Visual.setPropertyValue('DataRangeMin',-Inf);
        handles.Bmovie.Visual.setPropertyValue('DataRangeMax',Inf);
    else
        handles.Fmovie = implay(handles.deltaF(:,:,1,:));        
        handles.Bmovie.Visual.setPropertyValue('UseDataRange',true);
        handles.Bmovie.Visual.setPropertyValue('DataRangeMin',-Inf);
        handles.Bmovie.Visual.setPropertyValue('DataRangeMax',Inf);
    end
else    
    % Close already open animation window
    if ishandle(handles.Bmovie)
        close(handles.Bmovie);
    end    
end    
guidata(hObject, handles)

function minFdistTxt_Callback(hObject, eventdata, handles)

function minFdistTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxFdistTxt_Callback(hObject, eventdata, handles)

function maxFdistTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function distFbin_Callback(hObject, eventdata, handles)

function distFbin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function computeFstats_Callback(hObject, eventdata, handles)
% compute statistics and distribution of binned fluorescence

allF = handles.workImages(:);
temp = allF;
handles.medianF = median(temp(~isnan(temp)));       % stats of the entire sample
handles.stdF = std(temp(~isnan(temp)));

bins = str2double(get(handles.distFbin,'String'));
handles.Fleft = str2double(get(handles.minFdistTxt,'String'));
handles.Fright = str2double(get(handles.maxFdistTxt,'String'));
xValues = handles.Fleft:(handles.Fright-handles.Fleft)/(bins-1):handles.Fright;

[handles.Fnel,handles.Fcent] = hist(allF,xValues);
handles.FstatsFlag = 1;
guidata(hObject, handles);

plotDistribution(hObject,handles)

function computeDFstats_Callback(hObject, eventdata, handles)
% Compute the statistics and distribution of the DF/F signal
if get(handles.previewFthr,'Value'),
    % if the preview mode is active the statistics are computed only on the
    % NON masked data
    temp = handles.deltaF;
    binmask = handles.Fmask + 1;
    binmask(binmask==2) = NaN;    
    for i=1:handles.metadata.frames
        temp(:,:,1,i) = temp(:,:,1,i) .* binmask;
    end
    allDF = temp(:);
else
    allDF = handles.deltaF(:);
end    
temp = allDF;
handles.medianDF = median(temp(isfinite(temp)));       % stats of the entire sample
handles.stdDF = std(temp(isfinite(temp)));

strOut = ['DF/F median ' num2str(handles.medianDF,'%2.3g')];
strOut = [strOut ' std ' num2str(handles.stdDF,'%2.3g')];
strOut = [strOut '  Median, STD: ' num2str(handles.medianDF,'%1.2g') ' ' num2str(handles.stdDF,'%1.2g')];
set(handles.statDFtxt,'String',strOut);

bins = str2double(get(handles.binDistTxt,'String'));
handles.DFleft = str2double(get(handles.minDistTxt,'String'));
handles.DFright = str2double(get(handles.maxDistTxt,'String'));
xValues = handles.DFleft:(handles.DFright-handles.DFleft)/(bins-1):handles.DFright;

[handles.DFnel,handles.DFcent] = hist(allDF,xValues);
guidata(hObject, handles);
handles.DFstatsFlag = 1;
plotDistribution(hObject,handles)

guidata(hObject, handles);

function uibuttongroupDistribution_SelectionChangedFcn(hObject, eventdata, handles)
plotDistribution(hObject,handles);

function thrFtxt_Callback(hObject, eventdata, handles)

function thrFtxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function thresholdFck_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of thresholdFck
if get(hObject,'Value'),
    % Compute the thresholded image stack
    % First, save the workImages stack in a buffer for rapid recovery
    handles.buffer = handles.workImages;
    % Compute the mask
    thr = str2double(get(handles.thrFtxt,'String'));
    handles.binMask = handles.mdImage;
    handles.binMask(handles.binMask<=thr) = 0;
    % all masked elements NOT = 0 are made equal to 1
    handles.binMask(handles.binMask>0) = 1;
    % This mask is multiplied with the binarized stack BEFORE spatial and
    % temporal filtering.
    handles.binMask(isnan(handles.binMask)) = 0;
else
    % restore the original data
    handles.workImages = handles.buffer;    
end

guidata(hObject,handles)

function previewFthr_Callback(hObject, eventdata, handles)
% Here we compute a alpha channel by comparing the median fluorescence with the threshold

if (get(hObject,'Value')==1),
    thr = str2double(get(handles.thrFtxt,'String'));
    handles.Fmask = handles.mdImage;
    handles.Fmask(handles.Fmask<=thr) = 0;
    handles.Fmask(handles.Fmask>0) = 1;
    handles.Fmask = -(handles.Fmask-1);
    % the mask holds 1 where the fdata is invalid and he mask is opaque and 0
    % were the data is valid and the mask is transparent
    handles.FmaskFlag = 1;

    handles.black = handles.mdImage;
    handles.black = handles.black * 0;
    guidata(hObject);
    displayImageExplore (hObject, handles);
else
    handles.FmaskFlag = 0;
end    

guidata(hObject, handles);

function minSizeTxt_Callback(hObject, eventdata, handles)

function minSizeTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function computeGaussMap_Callback(hObject, eventdata, handles)
% This function computes the Gauss map of the imaging stack. The map is
% obtained by computing a 'distance' between the unitary gaussian distribution and
% the cdf of each pixel
% The distance is computed using two different metrics:
% mode = 1, 2 the probability p that each pixel process belongs to
% a normal distribution. The larger is p the darker is the map.
% mode = 3, 4 the area between the pixel cdf and the unitary gaussian
% distribution

% Computationally speaking the longest task is the computation of the stats
% of each process Once that is done, the computation of the metrics is
% straighforward. Therefore, I compute the four metrics at once.

% GMR April, 2018.

tic
% preallocate output maps. The CDF matrix makes use of the a priori
% understanding that the bin number is 1401. To be fixed if we extend the
% range. The present range is -6 +8 standard deviation units that, quite
% frankly, should be more than enough.

handles.ABmap.KSRB = zeros(handles.WIrows, handles.WIcolumns);
handles.ABmap.KSstdRB = zeros(handles.WIrows, handles.WIcolumns);
handles.ABmap.areaRB = zeros(handles.WIrows, handles.WIcolumns);
handles.ABmap.areaStdRB = zeros(handles.WIrows, handles.WIcolumns);
handles.CDFmap =  zeros(handles.WIrows, handles.WIcolumns, 1401);

% first, read the degree of freedom to be used for the KS test (mode = 1 or 2)
nFreedom = str2double(handles.KSnTxt.String);
% second, initialize the ktestGauss function. This is performed by calling
% the function with the first field different from 0. 
[ks, ksstd, area, areastd, cdf] = kstestGauss(1, [],-6,8,nFreedom);
source = handles.ABsourceF.Value;
for column = 1:handles.WIcolumns
    display(['Process column ' num2str(column)]);
    for row = 1:handles.WIrows
        % extract the process to analyze
        if source
            process = handles.workImages(row,column,:);
        else
            process = handles.deltaF(row,column,:);
        end            
        process = permute(process,[3 2 1]);     % transform process in a row vector
        [ks, ksstd, area, areastd, cdf] = kstestGauss(0, process,-6,8,nFreedom);
        handles.ABmap.KSRB (row,column) = ks;
        handles.ABmap.KSstdRB (row,column) = ksstd;
        handles.ABmap.areaRB (row,column) = area;
        handles.ABmap.areaStdRB (row,column) = areastd;
        handles.CDFmap(row,column,:) = cdf(:);

    end
end
% Normalisation of the maps (except ks) between 0 and 1
handles.ABmap.KSstdRB = handles.ABmap.KSstdRB/max(max(handles.ABmap.KSstdRB));
handles.ABmap.areaRB = handles.ABmap.areaRB/max(max(handles.ABmap.areaRB));
handles.ABmap.areaStdRB = handles.ABmap.areaStdRB/max(max(handles.ABmap.areaStdRB));
% now set the default values for the display limits
handles.ABmapLimits.KSRB = [0 1];
handles.ABmapLimB.KSRB = [0 1];
handles.ABmapLimits.KSstdRB = [0 1];
handles.ABmapLimB.KSstdRB = [0 1];
handles.ABmapLimits.areaRB = [0 1];
handles.ABmapLimB.areaRB = [0 1];
handles.ABmapLimits.areaStdRB = [0 1];
handles.ABmapLimB.areaStdRB = [0 1];

onABnormalMap(hObject, handles);
e = toc;
display(['End of AB-map computation: ' num2str(e) ' sec']);

function KSnTxt_Callback(hObject, eventdata, handles)

function KSnTxt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ABmin_Callback(hObject, eventdata, handles)
% make sure that an AB-map is selected for display in the imageExplore window
mn = str2double(get(hObject, 'String'));
wndSel = handles.imageFrameSource.SelectedObject.Tag;
if isempty(strfind('framesRBmaxRBaverageRB',wndSel))
   % AB-map selected
   handles.ABmapLimits.(wndSel)(1) = mn;
end 

guidata(hObject, handles)
displayImageExplore (hObject, handles)

function ABmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ABmax_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ABmax as text
%        str2double(get(hObject,'String')) returns contents of ABmax as a double

% make sure that an AB-map is selected for display in the imageExplore
% window
mx = str2double(get(hObject, 'String'));
wndSel = handles.imageFrameSource.SelectedObject.Tag;
if isempty(strfind('framesRBmaxRBaverageRB', wndSel))
   % AB-map selected
   handles.ABmapLimits.(wndSel)(2) = mx;
end 

guidata(hObject, handles)
displayImageExplore (hObject, handles)

function ABmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ABminb_Callback(hObject, eventdata, handles)
% make sure that an AB-map is selected for display in the binaryExplore window
mn = str2double(get(hObject, 'String'));
wndSel = handles.binaryImageSource.SelectedObject.Tag;
if isempty(strfind('FbinRBbinMaxRB',wndSel))
    % AB-map selected
    wndSel=wndSel(1:end-1);         % remove last char
    handles.ABmapLimB.(wndSel)(1) = mn;
end 

guidata(hObject, handles)
displayBinaryImage (hObject, handles)

function ABminb_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ABmaxb_Callback(hObject, eventdata, handles)
% make sure that an AB-map is selected for display in the binaryExplore window
mx = str2double(get(hObject, 'String'));
wndSel = handles.binaryImageSource.SelectedObject.Tag;
if isempty(strfind('FbinRBbinMaxRB',wndSel))
    % AB-map selected
    wndSel=wndSel(1:end-1);         % remove last char
    handles.ABmapLimB.(wndSel)(2) = mx;
end 

guidata(hObject, handles)
displayBinaryImage (hObject, handles)

function ABmaxb_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LFPcorrelate_Callback(hObject, eventdata, handles)

handles = guidata(hObject);  % refresh handles. It should not be necessary, but better be safe than sorry.
if ~isfield(handles,'metadata'), errNow = errordlg('Load an imaging data first','Monkey Error');
else
    if ~isfield(handles,'Xcorr') % check if the child process already exists by checking the existance of the field. 
        % if it does not exist the call to Xcorr must start no matter of
        % what
        %    if ishghandle(handles.Xcorr)==[]
        % open the child process
        handles.Xcorr = loadEPhys(gcf,handles);      % handle to the child figure
                                                     % arguments: parent object and
                                                     % parent handles
        handles.XcorrData = guidata(handles.Xcorr);  % and its controls
        hObject.Value = 1;                           % make sure that the checkbox correctly reflects the child status
    else
        % two possibilites: it exist and then it must be closed. or it
        % exists but it has been deleted. In that case it must be restarted
        % somehow
        if ~isvalid(handles.Xcorr) % now check whether is a valid handle.
            handles.Xcorr = loadEPhys(gcf,handles);      % handle to the child figure
                                                         % arguments: parent object and
                                                         % parent handles
            handles.XcorrData = guidata(handles.Xcorr);  % and its controls
            hObject.Value = 1;                           % make sure that the checkbox correctly reflects the child status        
        else
            % close the child process
            close(handles.Xcorr);
            delete(handles.Xcorr);      % delete the handle to the child process
            hObject.Value = 0;
        end
    end
end    
guidata(hObject, handles);

function RGBcoding_Callback(hObject, eventdata, handles)

% Display 


% --- Executes on button press in BMapp_pb.
function BMapp_pb_Callback(hObject, eventdata, handles)
% hObject    handle to BMapp_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
applicationBinaryMap(handles)


% --- Executes on button press in stamp1.
%Emanuele 24-10-2018 open a new figure with the current image displayed on
%imageExplore and saves automatically as 'image1' (jpg file)
function stamp1_Callback(hObject, eventdata, handles)
% hObject    handle to stamp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
currentimage=getimage(handles.imageExplore);
imshow(currentimage);
saveas(gcf,'image1.jpg', 'jpg');

%displayImageExplore (hObject, handles)
%dislay(ImageExplore(hObject, handles));


% --- Executes on button press in stamp2.
%Emanuele 24-10-2018 open a new figure with the current image displayed on
%binaryExplore and saves automatically as 'image2' (jpg file)
function stamp2_Callback(hObject, eventdata, handles)
% hObject    handle to stamp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
currentimage=getimage(handles.binaryExplore);
imshow(currentimage); 
saveas(gcf,'image2.jpg', 'jpg');


% --- Executes on button press in bintomat.
function bintomat_Callback(hObject, eventdata, handles)
% hObject    handle to bintomat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
binaryguide=handles.erodedBinaryImages;
save('binaryguide.mat','binaryguide');



function minsizefilt_Callback(hObject, eventdata, handles)
% hObject    handle to minsizefilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minsizefilt as text
%        str2double(get(hObject,'String')) returns contents of minsizefilt as a double


% --- Executes during object creation, after setting all properties.
function minsizefilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minsizefilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxsizefilt_Callback(hObject, eventdata, handles)
% hObject    handle to maxsizefilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxsizefilt as text
%        str2double(get(hObject,'String')) returns contents of maxsizefilt as a double


% --- Executes during object creation, after setting all properties.
function maxsizefilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxsizefilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
