function varargout = applicationBinaryMap(varargin)
% APPLICATIONBINARYMAP MATLAB code for applicationBinaryMap.fig
%      APPLICATIONBINARYMAP, by itself, creates a new APPLICATIONBINARYMAP or raises the existing
%      singleton*.
%
%      H = APPLICATIONBINARYMAP returns the handle to a new APPLICATIONBINARYMAP or the handle to
%      the existing singleton*.
%
%      APPLICATIONBINARYMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APPLICATIONBINARYMAP.M with the given input arguments.
%
%      APPLICATIONBINARYMAP('Property','Value',...) creates a new APPLICATIONBINARYMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before applicationBinaryMap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to applicationBinaryMap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help applicationBinaryMap

% Last Modified by GUIDE v2.5 29-Oct-2018 13:38:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @applicationBinaryMap_OpeningFcn, ...
                   'gui_OutputFcn',  @applicationBinaryMap_OutputFcn, ...
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


% --- Executes just before applicationBinaryMap is made visible.
function applicationBinaryMap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to applicationBinaryMap (see VARARGIN)

% Choose default command line output for applicationBinaryMap
handles.output = hObject;

% establish link with the calling process
handles.parentHandles = varargin{1};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes applicationBinaryMap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = applicationBinaryMap_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.output = hObject;


% --- Executes on button press in applyBinMask_pb.
function applyBinMask_pb_Callback(hObject, eventdata, handles)
% hObject    handle to applyBinMask_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%ima=getimage(handles.parentHandles.erodedBinaryImages(:,:,1,:));
%ima=imbinarize(ima);
%handles.parentHandles.erodedBinaryImages(:,:,1,:)=logical(handles.parentHandles.erodedBinaryImages(:,:,1,:));
%handles.parentHandles.erodedBinaryImages(:,:,:,:)=double(handles.parentHandles.erodedBinaryImages(:,:,:,:));
if get(handles.parentHandles.prairieRB,'Value')
handles.BinRoiDF(:,:,1,:)= handles.parentHandles.erodedBinaryImages(:,:,1,:) .* handles.parentHandles.deltaF(:,:,1,:);
    handles.currentFrame = 1;
end

if get(handles.parentHandles.TIFFmultiImageRB,'Value')
handles.BinRoiDF(:,:,1,:)= handles.parentHandles.erodedBinaryImages(:,:,1,:) .* handles.parentHandles.workImages(:,:,1,:);
    handles.currentFrame = 1;
end

% compute the limits of the image range
handles.BinRoiMin = 0;
handles.BinRoiMax = max(max(max(squeeze(handles.BinRoiDF(:,:,1,:)))));
handles.meanBinRoi = (handles.BinRoiMin + handles.BinRoiMax)/2;
         
tmp = handles.BinRoiDF(:,:,:,handles.currentFrame);
axes(handles.axes5)
imagesc(tmp,[0 handles.BinRoiMax]);
%handles.mdImageW = imagesc(handles.ABmap.(wndSel),'ButtonDownFcn',@mdImageClick);
colormap hot;       
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in back1frameBM.
function back1frameBM_Callback(hObject, eventdata, handles)
  
if handles.currentFrame == 1
    handles.currentFrame = handles.parentHandles.metadata.frames;
else
    handles.currentFrame = handles.currentFrame - 1;
end
set(handles.parentHandles.currentFrTxt,'String',handles.currentFrame);
guidata(hObject, handles);
handles.axes5 = imagesc(handles.BinRoiDF(:,:,:,handles.currentFrame),[0 handles.BinRoiMax]);


% --- Executes on button press in advace1frameBM.
function advace1frameBM_Callback(hObject, eventdata, handles)
if handles.currentFrame == handles.parentHandles.metadata.frames
    handles.currentFrame = 1;
else
    handles.currentFrame = handles.currentFrame + 1;
end
set(handles.parentHandles.currentFrTxt,'String',handles.currentFrame);
guidata(hObject, handles);
handles.axes5 = imagesc(handles.BinRoiDF(:,:,:,handles.currentFrame),[0 handles.BinRoiMax]);


% --- Executes on button press in SaveBM_dF_pb.
function SaveBM_dF_pb_Callback(hObject, eventdata, handles)
% hObject    handle to SaveBM_dF_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%actualmap = hot;
%imwrite(handles.BinRoiDF(:,:,:,1),actualmap,'BM_dF.tif','tiff')
%Emanuele 25-10-2018 save the file in AVI format
AVIoutdF=VideoWriter ('BM_dF.avi','Uncompressed AVI');
open(AVIoutdF);

hFig = figure('units','pixels','position',[200 200 256 256]);
hMovie = axes('units','pixels','position',[0 0 256 256],'Visible','off');
%axis ([0 handles.WIrows 0 handles.WIcolumns]);
set(hMovie, 'XTick', [], 'YTick', []);

%video=handles.BinRoiDF(:,:,:,1);
axes(hMovie);
if get(handles.parentHandles.prairieRB,'Value')
    for  cf=1:handles.parentHandles.metadata.frames
        %BinRoiDFsave=zeros(200,200,256,256);
        handles.BinRoiMax = max(max(max(squeeze(handles.BinRoiDF(:,:,1,:)))));
        BinRoiDFsave(:,:,1,cf)= handles.parentHandles.erodedBinaryImages(:,:,1,cf) .* handles.parentHandles.deltaF(:,:,1,cf);
tmp = BinRoiDFsave(:,:,:,cf);
           imagesc(tmp,[0 handles.BinRoiMax]);
           colormap hot;
        %cf=2:handles.parentHandles.metadata.frames
        %imwrite(handles.BinRoiDF(:,:,:,cf),'BM_dF.tif','WriteMode','append')
        %video=(handles.BinRoiDF(:,:,:,cf));
        %writeVideo(AVIoutdF,cf);
        
        frame =(handles.BinRoiDF(:,:,:,cf));
    [row, column] = find(frame);
   
    %axis ([0 handles.WIrows 0 handles.WIcolumns]);
    %set(hMovie, 'XTick', [], 'YTick', []);
    refresh;
%    myavi = addframe(myavi, hFig);
    im = frame2im(getframe(hFig));
    %im = rgb2gray(im);
    writeVideo(AVIoutdF,im);
    end
end
if get(handles.parentHandles.TIFFmultiImageRB,'Value')
     for  cf=1:handles.parentHandles.metadata.frames
        %BinRoiDFsave=zeros(200,200,256,256);
        handles.BinRoiMax = max(max(max(squeeze(handles.BinRoiDF(:,:,1,:)))));
        BinRoiDFsave(:,:,1,cf)= handles.parentHandles.erodedBinaryImages(:,:,1,cf) .* handles.parentHandles.workImages(:,:,1,cf);
tmp = BinRoiDFsave(:,:,:,cf);
           imagesc(tmp,[0 handles.BinRoiMax]);
           colormap hot;
        %cf=2:handles.parentHandles.metadata.frames
        %imwrite(handles.BinRoiDF(:,:,:,cf),'BM_dF.tif','WriteMode','append')
        %video=(handles.BinRoiDF(:,:,:,cf));
        %writeVideo(AVIoutdF,cf);
        
        frame =(handles.BinRoiDF(:,:,:,cf));
    [row, column] = find(frame);
   
    %axis ([0 handles.WIrows 0 handles.WIcolumns]);
    %set(hMovie, 'XTick', [], 'YTick', []);
    refresh;
%    myavi = addframe(myavi, hFig);
    im = frame2im(getframe(hFig));
    %im = rgb2gray(im);
    writeVideo(AVIoutdF,im);
     end
end
    display('end of BM_dF');
    close(AVIoutdF);
    

% --- Executes on button press in loadROI.
function loadROI_Callback(hObject, eventdata, handles)
% hObject    handle to loadROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[selection,path_in,filterI] = uigetfile ('*.*','Aprila!!!','C:\Users\SNS\Desktop');
file=xlsread(selection,'B:B');
ROI=file';
%meanROI=mean(ROI);
%varROI=var(ROI);
handles.ROI=normalize(ROI);
display('ROI loaded');
guidata(hObject, handles);


% --- Executes on button press in CrossCROI.
function CrossCROI_Callback(hObject, eventdata, handles)
% hObject    handle to CrossCROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.parentHandles.deltaF(isnan(handles.parentHandles.deltaF(:,:,1,:)))=0;
for i=1:size(handles.parentHandles.erodedBinaryImages(:,:,1,:),4)
    
handles.BinRoiDFcorr(:,:,1,i)= handles.parentHandles.erodedBinaryImages(:,:,1,i) .* handles.parentHandles.deltaF(:,:,1,i);

i=i+1;
end
corrWin = str2double(handles.corrwindow.String); %Window for the computation of Xcorrelation spectra
intWin = str2double(handles.intwindow.String); %Integration window for the computation of Xcorrelation measure
% these data must be trasformed in number of frames
corrWin = floor(corrWin/handles.meanBinRoi);
intWin = floor(intWin/(2*handles.meanBinRoi));
intFrom = corrWin-intWin;
intTo = corrWin+intWin;
handles.correlationarray=zeros(128,128);
for i=1:size(handles.BinRoiDFcorr(:,:,:,:),1)
    for j=1:size(handles.BinRoiDFcorr(:,:,:,:),2)
        pixelts=squeeze(handles.BinRoiDFcorr(i,j,1,:));
        pixelts=pixelts';
     corr=xcorr(handles.ROI,pixelts,corrWin);
        corrsum=sum(corr(intFrom:intTo));
        handles.correlationarray(i,j,1)=corrsum;
       
        j=j+1;
    end
     i=i+1;
end
axes(handles.axes4)
imagesc(handles.correlationarray);
colormap hot;
set(handles.axes4, 'Visible', 'off');
guidata(hObject, handles);
%handles.plotCC=


% --- Executes on button press in ccROI.
function ccROI_Callback(hObject, eventdata, handles)
% hObject    handle to ccROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure;
currentimage=getimage(handles.axes4);
imagesc(currentimage);
colormap hot;
saveas(gcf,'ccROI.jpg', 'jpg');

function intwindow_CreateFcn(hObject, eventdata, handles)



function windowdisplacement_Callback(hObject, eventdata, handles)
% hObject    handle to windowdisplacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of windowdisplacement as text
%        str2double(get(hObject,'String')) returns contents of windowdisplacement as a double


% --- Executes during object creation, after setting all properties.
function windowdisplacement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowdisplacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slidingwindow_Callback(hObject, eventdata, handles)
% hObject    handle to slidingwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slidingwindow as text
%        str2double(get(hObject,'String')) returns contents of slidingwindow as a double


% --- Executes during object creation, after setting all properties.
function slidingwindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidingwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stCCROI.
function stCCROI_Callback(hObject, eventdata, handles)
% hObject    handle to stCCROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.parentHandles.deltaF(isnan(handles.parentHandles.deltaF(:,:,1,:)))=0;
for i=1:size(handles.parentHandles.erodedBinaryImages(:,:,1,:),4)
    
handles.BinRoiDFcorr(:,:,1,i)= handles.parentHandles.erodedBinaryImages(:,:,1,i) .* handles.parentHandles.deltaF(:,:,1,i);

i=i+1;
end
handles.totalframes=i;
slidewindow=str2double(handles.slidingwindow.String);%time window for the computation of the short time CC
dispwindow=str2double(handles.windowdisplacement.String); %window displacement (in terms of frames) for the short time CC
corrWin = str2double(handles.corrwindow.String); %Window for the computation of Xcorrelation spectra
intWin = str2double(handles.intwindow.String); %Integration window for the computation of Xcorrelation measure
% these data must be trasformed in number of frames
corrWin = floor(corrWin/handles.meanBinRoi);
intWin = floor(intWin/(2*handles.meanBinRoi));
intFrom = corrWin-intWin;
intTo = corrWin+intWin;
t=size(handles.parentHandles.erodedBinaryImages(:,:,1,:),4);
handles.stcorrelationarray=zeros(256,256,1,t);
if corrWin<slidewindow
    for z=1:(t-slidewindow)
        if slidewindow<=t
            for i=1:size(handles.BinRoiDFcorr(:,:,:,:),1)
                for j=1:size(handles.BinRoiDFcorr(:,:,:,:),2)
                    stCCpixelts=squeeze(handles.BinRoiDFcorr(i,j,1,z:slidewindow));
                    stCCpixelts=stCCpixelts';
                    actualROI=handles.ROI(1,z:slidewindow);
                    corr=xcorr(actualROI,stCCpixelts,corrWin);
                    corrsum=sum(corr(intFrom:intTo));
                    handles.stcorrelationarray(i,j,1,z)=corrsum;
                    
                    j=j+1;
                end
                i=i+1;
            end
        end
        z=z+dispwindow %#ok<NOPRT>
        slidewindow=slidewindow+dispwindow;
    end
else
    display('The time window must be smaller than the sliding window. Try again.');
end
handles.stcorrarrayMin = 0;
handles.stcorrarrayMax = max(max(max(squeeze(handles.BinRoiDF(:,:,1,:)))));
handles.meanstcorrarray = (handles.stcorrarrayMin + handles.stcorrarrayMax)/2;
load handel; %#ok<LOAD>
player = audioplayer(y, Fs);
play(player);
handles.actualframe=1;
guidata(hObject, handles);
showstCC (handles);


function showstCC (handles)
axes(handles.axes6);
imagesc(handles.stcorrelationarray(:,:,:,handles.actualframe),[0 handles.stcorrarrayMax]);
set(handles.axes6, 'Visible', 'off');
colormap jet;



% --- Executes on button press in backward.
function backward_Callback(hObject, eventdata, handles)
% hObject    handle to backward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.actualframe>=1
handles.actualframe=handles.actualframe+1;
guidata(hObject,handles);
end
showstCC (handles);



% --- Executes on button press in forward.
function forward_Callback(hObject, eventdata, handles)
% hObject    handle to forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.actualframe<=handles.totalframes
handles.actualframe=handles.actualframe+1;
guidata(hObject,handles);
end
showstCC (handles);

function corrwindow_CreateFcn(hObject, eventdata, handles)

function axes5_CreateFcn(hObject, eventdata, handles)

function axes4_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in savestccROI.
function savestccROI_Callback(hObject, eventdata, handles)
% hObject    handle to savestccROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AVIoutstcc=VideoWriter ('st-ccROI.avi','Uncompressed AVI');
open(AVIoutstcc);

hFig = figure('units','pixels','position',[200 200 256 256]);
hMovie = axes('units','pixels','position',[0 0 256 256],'Visible','off');
%axis ([0 handles.WIrows 0 handles.WIcolumns]);
set(hMovie, 'XTick', [], 'YTick', []);
axes(hMovie);

    for  cf=1:size(handles.stcorrelationarray,4)
        %BinRoiDFsave=zeros(200,200,256,256);
        tmp=handles.stcorrelationarray(:,:,:,cf);
           imagesc(tmp,[0 handles.stcorrarrayMax]);
           colormap hot;
        %cf=2:handles.parentHandles.metadata.frames
        %imwrite(handles.BinRoiDF(:,:,:,cf),'BM_dF.tif','WriteMode','append')
        %video=(handles.BinRoiDF(:,:,:,cf));
        %writeVideo(AVIoutdF,cf);
        
        frame =(handles.stcorrelationarray(:,:,:,cf));
    [row, column] = find(frame);
   
    %axis ([0 handles.WIrows 0 handles.WIcolumns]);
    %set(hMovie, 'XTick', [], 'YTick', []);
    refresh;
%    myavi = addframe(myavi, hFig);
    im = frame2im(getframe(hFig));
    %im = rgb2gray(im);
    writeVideo(AVIoutstcc,im);
    end
    display('end of BM_dF');
    close(AVIoutdF);
    
