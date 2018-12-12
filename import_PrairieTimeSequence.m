function [xyct_img metadata] = import_PrairieTimeSequence(varargin)
% SYNTAX: [xyct_img metadata] = import_PrairieTimeSequence(varargin)
%   import_PrairieTimeSequence('C:\Images\SingleScan001-1')
%
% METHOD: Imports an image taken from a prairie CLSM into MATLAB memory
% with an axis convention of xyct.  By default, removes any singleton 
% dimensions (see ALGORITHM). Also parses some metadata from the xml file.
%
% INPUT: 
% img_full_path [o, string]: path to the folder containing prairie tif 
%   files, xml file, and config file.  If not provided, prompts user to 
%   browse.
% channel_index [p, 1xnchan numeric]: index which reorder channels, the
%   value of each element is the index for each channel.  Extra channels 
%   are zero padded (i.e. can be used to force 2 channel image to be RGB by
%   specifying [1 2 3]).  Defaults to not reordering index or padding to
%   RGB.
% DataType [p, class]: the image will be scaled and casted to this
%   class denoted by this string. Default: uint16.
% SQUEEZE_DIMS [p, bool]: removes singleton dimensions of imported image.  
%   Default = true.
%
% OUTPUT:
% xyct_img [img]: imported image, singleton dimensions may be
%   removed.
% metadata [struct]: a MATLAB struct containing a subset of the metadata
%   found in the xml file.  These are gleaned using regular expressions.  
%   See below for explanation.
%
% EXAMPLE:
%  [xyct_img, metadata] = import_PrairieTimeSequence('C:\Images\ZSerie-015');
%  [xyct_img, metadata] = import_PrairieTimeSequence();
%
% ALGORITHM: 
% The img metadata is parsed from the xml file included with each prairie
% image.  Instead of parsing the entire xml file (which takes WAY TOO LONG 
% in MATLAB, even for small xml files), relevant metadata is gleaned using 
% regular expressions in token mode.  Code can easily be added to get more 
% metadata if desired.  
%   After loading the metadata, the individual prairie tif files are loaded
% into a 5D image of the form xyct.  For our system, Prairie images are 12
% bit.  Since there is no native 12 bit dataype in matlab the images are 
% cast and scaled to 16 bit.  This behavior can be changed with the INPUT
% parameters.
%   By default (SQUEEZE_DIMS), singleton dimensions of the image are 
% removed, so a single image would be xyc  and a time series single image would be xyct. 
% Disabling this will return a xyct image regardless of whether those 
% dimensions are singleton or not.
%
% ASSUMPTIONS:
% The code assumes the bit depth and pixel dimensions remain constant
% across images.
%
% NOTES:
%  MATLAB Regular Expressions: mostly follows the standard regular expression
% syntax.  One exception found so far is '.' matches to any character,
% including newline.  So '.*' goes ACROSS newline characters, which is not
% default behavior. To make '.*' behave in correct way the easiest
% modification is '.*?' which means 'repeat matching of any char in a 
% RELUCTANT fashion ACROSS newline characters'.  So it matches to a single
% line if possible.
% 
%  Fields of 'metadata' structure.
%   .xdim
%   .ydim
%   .x_stg: x axis stage coordinate in um
%   .y_stg: y axis stage coordinate in um
%   .middle_z_stg: z axis stage coordinate of the middle z slice
%   .dx_um: x axis um/pixel ratio
%   .dy_um: y axis um/pixel ratio
%   .dz_um: um distance between z slices in um, NAN if single image.
%
%
% DEPENDENCIES: none
%   Tested MATLAB R2011a win32
%
% AUTHOR: Bruce A. Corliss, Scientist III
% REVISIONS: GMR, February-June 2017.
% Adapted to new Prairie format. The code has been modified to disable
% input of z, zt series. Realistically, I will never need this and the
% original code was not working with the current Prairie xml syntax.
% The structure of the metadata file has been heavily modified in the
% latest Prairie releases and the code has been modified accordingly. These
% changes will certainly affect older data.

% LICENSE: BSD-3.  See Below.
%
% ORGANIZATION: GrassRoots Biotechnology
%
% CONTACT:
%  bruce.a.corliss@gmail.com
%  bruce.corliss@grassrootsbiotech.com
%
% DATE: 6/15/2011
%
% VERSION: 1.0.1

% LICENSE
% Copyright (c) 2011, GrassRoots Biotechnology
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%     Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%     Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     Neither the name of the GrassRoots Biotechnology nor the names of its 
%       contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
% AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
% USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% ARGUMENT PARSING
p = inputParser;
p.addOptional('img_full_path', '',@ischar);
p.addParamValue('channel_index', [], @isnumeric);
p.addParamValue('DataType', 'uint16', @char);
p.addParamValue('SQUEEZE_DIMS', 1, @(x)(x==0 || x==1));
p.StructExpand = false;  
p.parse(varargin{:});
% Import parsed variables into workspace
fargs = fields(p.Results);
for n=1:numel(fargs); eval([fargs{n} '=' 'p.Results.' fargs{n} ';']);  end

% Prompt for manual browse to image folder if null input
if isempty(img_full_path); 
    img_full_path  = uigetdir('Select Folder containing Prairie Image');
    if(img_full_path ==0); return; end
end

%% Make path *nix compatible, extract path and file name of image data
img_full_path = regexprep(img_full_path, '\', '/');
[pathStr, img_name] = fileparts(img_full_path);

%% If image folder or tif or XML files DNE, return empty
if isempty(dir(pathStr)) || numel(dir([pathStr '/*.tif']))==0 ...
       || numel(dir([pathStr '/*.xml']))==0; 
   xyct_img = []; return;
end

%% Read in metadata from xml file
metafile = dir ([pathStr '/*.xml']);
%text = fileread([img_full_path img_name '.xml']);
text = fileread([pathStr '/' metafile.name]);

%% Find dimensions of 5D image
% XY dimensions: it is assumed that all frames have the same scale
metadata.xdim = str2double(regexp(text, 'key="pixelsPerLine".*?value="(\d*)"', 'tokens','once'));
metadata.ydim = str2double(regexp(text, 'key="linesPerFrame".*?value="(\d*)"', 'tokens','once'));

% extract infos about objective and zoom
metadata.objective = regexp(text, 'key="objectiveLens".*?value="(.*?)"', 'tokens','once');
metadata.objMag = str2double(regexp(text, 'key="objectiveLensMag".*?value="(\d*)"', 'tokens','once'));
metadata.objNA = str2double(regexp(text, 'key="objectiveLensNA".*?value="(\d*)"', 'tokens','once'));
metadata.zoom = str2double(regexp(text, 'key="opticalZoom".*?value="(\d*)"', 'tokens','once'));

% Revision June 03, 2017 (GMR). The original code did not interpret correctly the time depth of the sequence. 
% each frame comes with two entries 'Frame relative...' but the index is
% always set to 1. 
% t_cell =  regexp(text, 'Frame relative.*?index="(\d*)"', 'tokens');
% t_cell =  regexp(text, 'Frame relative.*?absoluteTime="(\d*)"', 'tokens');
% tpoints = cellfun(@(x) str2double(x{1}),t_cell);

%% Revision September 11, 2017 (GMR). 
%  the following code loads in the array tpoint the timing of each frame
pointerToTime = strfind(text,'absoluteTime=');      
pointerToTime = pointerToTime + 13;     % Array of pointers to the beginning of the time stamp string
nframes = length (pointerToTime);

for i=1:nframes
   shortText = text(pointerToTime(i):pointerToTime(i)+25);
   shorterText = extractBetween(shortText,'"','"');
   % remove '"' from start and end of string
   shorterText = shorterText(2:end-1);
   metadata.tpoints(i) = str2double(shorterText);
end
metadata.frames = nframes;
metadata.meanSP = (metadata.tpoints(end)-metadata.tpoints(1))/(nframes-1);

% % Z axis dimension. Original code commented out.
% z_cell =  regexp(text, 'Frame relative.*?index="(\d*)"', 'tokens');
% zslice = cellfun(@(x) str2double(x{1}),z_cell);
% T axis dimension
% t_cell =  regexp(text, 'Sequence type=.*?cycle="(\d*)"', 'tokens');
% tpoints = cellfun(@(x) str2double(x{1}),t_cell);
% if tpoints == 0; tpoints=1; end

% Ch_img contains the file name of each image for each channel as an ordered table
ch_img_names_cell =  regexp(text, 'File channel="(\d*)".*?filename="(.*?)"', 'tokens');
% ch specify the channel of each file in the table above
ch = cellfun(@(x) str2double(x{1}),ch_img_names_cell);
% and the next is a structure containing the file names
img_names = cellfun(@(x) x{2},ch_img_names_cell, 'UniformOutput', 0);
metadata.ch = max(ch);      % number of channels

% Bit depth of img in filesystem. This line has been modified to reflect
% changes of the metafile sintax.
tif_bit_depth = str2double(regexp(text, 'key="bitDepth".*?value="(\d*)"', 'tokens','once'));


%% Parse extra subset of metadata
% XYZ Coordinates of each img

% This section has been completely rewritten on account of the different
% structure of the metafile. GMR, Sep 10, 2017.

% The xml file has at least two alternative syntaxes that are identified
% during parsing. The code has been modified to extract the correct data
% regardless of the format. GMR, September 23, 2017.

% shortText = extractBetween(text,'key="positionCurrent"','</PVStateValue>');
% shorterText = extractBetween(shortText,'index="XAxis"','index="YAxis"');
% if isempty(shorterText), metadata.x_stg = 0
% else
%     metadata.x_stg = mode(cellfun(@(x) str2double(x{1}), regexp(shorterText,'<SubindexedValue subindex="0".*?value="([-\d\.]*?)"', 'tokens')));
% end
% 
% shorterText = extractBetween(shortText,'index="YAxis"','index="ZAxis"');
% metadata.y_stg = mode(cellfun(@(x) str2double(x{1}), regexp(shorterText,'<SubindexedValue subindex="0".*?value="([-\d\.]*?)"', 'tokens')));
% 
% shorterText = extractBetween(shortText,'index="ZAxis"','</SubindexedValues>');
% metadata.middle_z_stg = median(unique(cellfun(@(x) str2double(x{1}),regexp(shorterText, '<SubindexedValue subindex="0".*?value="([-\d\.]*?)"', 'tokens'))));

shortText = extractBetween(text,'key="micronsPerPixel"','key="minVoltage"');
if isempty(shortText),
    % to cover for alternative format of xml file. This is the Pisa format
    shortText = extractBetween(text,'key="micronsPerPixel','key="pmtGain');
    metadata.dx_um = mode(cellfun(@(x) str2double(x{1}), regexp(shortText,'XAxis".*?value="([-\d\.]*?)"', 'tokens'))); 
    metadata.dy_um = mode(cellfun(@(x) str2double(x{1}), regexp(shortText,'YAxis".*?value="([-\d\.]*?)"', 'tokens'))); 
else
    % Padova format
    metadata.dx_um = mode(cellfun(@(x) str2double(x{1}), regexp(shortText,'index="XAxis".*?value="([-\d\.]*?)"', 'tokens'))); 
    metadata.dy_um = mode(cellfun(@(x) str2double(x{1}), regexp(shortText,'index="YAxis".*?value="([-\d\.]*?)"', 'tokens'))); 
end

%% Create 1x1 mapping between img_name, channel, timepoint
% Channel index
ch_ind = ch;  
if ~isempty(channel_index)
    if ~isempty(setdiff(ch_ind, channel_index)); error(['Specified channel index must'...
            ' contain all channel index indices parsed from xml file']);
    end
    % Reassign channel indices
    switch_cell = arrayfun(@(x) x == ch_ind, 1:numel(unique(ch)),'UniformOutput', 0);
    for n = 1:numel(switch_cell); ch_ind(n) = channel_index(n); end
end

% timpoint index need to be repeated to match img_names elements
step1 = [1:1:nframes];
step2 = repmat(step1,2,1);
t_ind = step2(:);   % this line folds the two columns of step2 in a single vector

%% Clear xml text form memory, not needed for final img import
clear text;

%% Initialize image specified datatype
xyct_img = zeros(metadata.xdim, metadata.ydim, max(horzcat(ch, channel_index)), nframes, DataType);

%% Read individual tif files into 5D img
% modified on Jan 14, 2018 to accomodate single channel data
% First: Ineed to know how many channels are contained in the series
if ch(1)==ch(2),
    % the image sequence contains only one channel. Therefore...
    t_ind = [1:nframes];
end

for n = 1:numel(img_names)
    % The following statement has been modified on Feb 13, 2017. The RDIVIDE
    % (./) operator has been used instead of the simple (/) operator.
    % img_names is a structure holding all the file names
    tmp = imread([pathStr '/' img_names{n}]);
    xyct_img(:,:,ch_ind(n), t_ind(n)) = tmp; %* double(intmax(DataType)./2^tif_bit_depth);
end

%% Remove extra singleton dimensions
% KJan 12, 2018. The squeeze has been removed. IF the original data has
% only one channel the aray would loose one dimension which messes up the
% remain of the processing pipeline.
%  if SQUEEZE_DIMS
%      xyct_img = squeeze(xyct_img);
%  end

%ch2only(:,:,:)=xyct_img(:,:,2,:);
%implay (ch2only);

end
