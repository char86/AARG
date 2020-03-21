function varargout = graphAARG(varargin)
% GRAPHAARG MATLAB code for graphAARG.fig
%      GRAPHAARG, by itself, creates a new GRAPHAARG or raises the existing
%      singleton*.
%
%      H = GRAPHAARG returns the handle to a new GRAPHAARG or the handle to
%      the existing singleton*.
%
%      GRAPHAARG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GRAPHAARG.M with the given input arguments.
%
%      GRAPHAARG('Property','Value',...) creates a new GRAPHAARG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before graphAARG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to graphAARG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help graphAARG

% Last Modified by GUIDE v2.5 29-Aug-2018 22:24:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @graphAARG_OpeningFcn, ...
                   'gui_OutputFcn',  @graphAARG_OutputFcn, ...
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


% --- Executes just before graphAARG is made visible.
function graphAARG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to graphAARG (see VARARGIN)

% Choose default command line output for graphAARG

osType = computer;
handles.osType = osType;

screenSize = get(0, 'ScreenSize');
hObject.Visible = 'Off';

screenWidth = screenSize(3); screenHeight = screenSize(4);
if screenWidth/screenHeight > 1.6
    widthScaleExtension = 1.3; heightScaleExtension = 1.1;
else
    widthScaleExtension = 1; heightScaleExtension = 1;
end
if screenWidth > screenHeight
    [hObject] = ResizeGUIForCurrentScreen(hObject, screenWidth, screenHeight, widthScaleExtension, heightScaleExtension);
end

bottomPositionScaleFactor = 0.95;
guiPosition = hObject.Position;
loweringValue = guiPosition(2)*(1-bottomPositionScaleFactor);
hObject.Position = [guiPosition(1), guiPosition(2)-loweringValue, guiPosition(3), guiPosition(4)];

handles.output = hObject;

suffixStrs = {}; Dirs = {}; 
if nargin == 0
    keywordStrs = {}; handles.keywordStrs = keywordStrs;
else
    keywordStrs = varargin{1};
    if ~isempty(keywordStrs)
        handles.keywordStrs = keywordStrs; numberOfKeywords = size(keywordStrs,2);
        for cWordIdx = 1 : numberOfKeywords
            if cWordIdx == 1
                handles.edit_keywords.String = keywordStrs{cWordIdx};
            else
                handles.edit_keywords.String = strcat(handles.edit_keywords.String,{' '},keywordStrs{cWordIdx});
            end
        end

        handles.edit_keywords.ForegroundColor = [0 0 0];
        handles.edit_keywords.FontAngle = 'normal';
    else
        handles.keywordStrs = {};
    end
end
handles.suffixStrs = suffixStrs; handles.Dirs = Dirs;

%CG: disable customization controls. These will be enabled only if the user
%checks the 'Customize start point distances' checkbox.
handles.pushbutton_getBranchIDs.Enable = 'off';
handles.edit_experimentDirectoryIndex.Enable = 'off';
handles.edit_branchIDs.Enable = 'off';
handles.edit_customStartPointDistances.Enable = 'off';

configStruct = struct('dontDisplayHistograms',0,...
                          'displayNValues',1,...
                          'xAxisMinHistogram',0,...
                          'yAxisMinHistogram',0,...
                          'xAxisMaxHistogram','?',...
                          'yAxisMaxHistogram','?',...
                          'binWidth','default',...
                          'dontDisplayMedianPlots',0,...
                          'displayTTestResult',1,...
                          'xAxisMinMedian',0,...
                          'yAxisMinMedian',0,...
                          'xAxisMaxMedian','?',...
                          'yAxisMaxMedian','?',...
                          'hideUnnecessaryTickLabels',0,...
                          'displayPeakAmplitude',1,...
                          'displayDecayTimeConstant',0,...
                          'lengthOfDendritePerGraph','default',...
                          'maximumDistance', []);
%CG: must display something. When user clicks the decay radio button,
%which is inactive, amplitude data radio button will switch off and
%vice versa. A minimum of one must be checked. 
handles.configStruct = configStruct;

%CG: ******testing only.
% dircells = load('test_dirs.mat'); textstring = load('test_text.mat');
% handles.Dirs = dircells.Dirs; %handles.edit_selectedExperiments.String = sprintf(textstring.text);
% keywordStrs = {'a_baseline', 'b_TTAP2'}; handles.keywordStrs = keywordStrs;
% handles.edit_keywords.String = 'a_baseline b_TTAP2';
%CG: ******testing only.

set( findall( hObject, '-property', 'Units' ), 'Units', 'Normalized' )
set( hObject, 'Resize','on' ) 
% Update handles structure

guidata(hObject, handles);

% UIWAIT makes graphAARG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = graphAARG_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_customizeStartPointDistances.
function checkbox_customizeStartPointDistances_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_customizeStartPointDistances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_customizeStartPointDistances
if hObject.Value == 1
    
    handles.pushbutton_getBranchIDs.Enable = 'on';
    handles.edit_experimentDirectoryIndex.Enable = 'on';
    handles.edit_branchIDs.Enable = 'on';
    handles.edit_customStartPointDistances.Enable = 'on';
    
elseif hObject.Value == 0
%CG: if the user unchecks the 'customize start point distances' checkbox, then 
%the customization controls will be disabled and all their associated
%variables will be reset. Graphs will be plotted with all start points set
%to the default distance value, which is 0. 

    
    handles.pushbutton_getBranchIDs.Enable = 'off';
    
    handles.edit_experimentDirectoryIndex.ForegroundColor = [0.502 0.502 0.502];
    handles.edit_experimentDirectoryIndex.FontAngle = 'italic';
    handles.edit_experimentDirectoryIndex.String = 'Experiment directory index';
    handles.edit_experimentDirectoryIndex.Enable = 'off';
    
    handles.edit_branchIDs.ForegroundColor = [0.502 0.502 0.502];
    handles.edit_branchIDs.FontAngle = 'italic';
    handles.edit_branchIDs.String = 'Branch IDs';
    handles.edit_branchIDs.Enable = 'off';
    
    handles.edit_customStartPointDistances.ForegroundColor = [0.502 0.502 0.502];
    handles.edit_customStartPointDistances.FontAngle = 'italic';
    handles.edit_customStartPointDistances.String = 'Custom start point distances';
    handles.edit_customStartPointDistances.Enable = 'off';
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_getBranchIDs.
function pushbutton_getBranchIDs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_getBranchIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

keywordStrs = handles.keywordStrs;
getBranchIDs(keywordStrs)
%CG: show the dendrite outlines created with the AARGII interface.
%getBranchIDs will show the branch ID number for each branch, which the
%user should enter into the branchIDs editable text box. 

function edit_experimentDirectoryIndex_Callback(hObject, eventdata, handles)
% hObject    handle to edit_experimentDirectoryIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_experimentDirectoryIndex as text
%        str2double(get(hObject,'String')) returns contents of edit_experimentDirectoryIndex as a double
expDirIndexString = handles.edit_experimentDirectoryIndex.String;
commaSeparationIdx = strfind(expDirIndexString,',');
if isempty(commaSeparationIdx)
%CG: in this case we can expect the user to have entered only one numeric
%format.
    expDirIndex = str2double(expDirIndexString);
elseif ~isempty(commaSeparationIdx)
%CG: in this case there is more than one numeric value contained within the string.
%Because str2double returns NaN with multiple numeric values within a
%single string, it is necessary to iteratively collect each numeric value
%using commaSeparationIdx.
    numberOfCommas = numel(commaSeparationIdx); expDirIndex = [];
    for cComma = 1 : numberOfCommas

%CG: check in front and behind.
           if cComma == numberOfCommas && cComma ~= 1
               expDirIndex = [expDirIndex,str2double(expDirIndexString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
               expDirIndex = [expDirIndex,str2double(expDirIndexString(commaSeparationIdx(cComma):end))];
           elseif cComma == numberOfCommas && cComma == 1
               expDirIndex = [str2double(expDirIndexString(1:commaSeparationIdx-1)), str2double(expDirIndexString(commaSeparationIdx+1:end))];
           elseif cComma == 1
               expDirIndex = [expDirIndex,str2double(expDirIndexString(1:commaSeparationIdx(1)-1))];
           else
               expDirIndex = [expDirIndex,str2double(expDirIndexString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
           end
           nanIndices = isnan(expDirIndex); expDirIndex(nanIndices) = [];
    end
end

if ~isnan(expDirIndex)
    handles.edit_experimentDirectoryIndex.ForegroundColor = [0 0 0]; 
    handles.edit_experimentDirectoryIndex.FontAngle = 'normal';
else
    handles.edit_experimentDirectoryIndex.String = 'Experiment directory index';
    handles.edit_experimentDirectoryIndex.ForegroundColor = [0.5020 0.5020 0.5020];
    handles.edit_experimentDirectoryIndex.FontAngle = 'italic';
end


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_experimentDirectoryIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_experimentDirectoryIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_branchIDs_Callback(hObject, eventdata, handles)
% hObject    handle to edit_branchIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_branchIDs as text
%        str2double(get(hObject,'String')) returns contents of edit_branchIDs as a double

branchIDsString = handles.edit_branchIDs.String;
commaSeparationIdx = strfind(branchIDsString,',');
if isempty(commaSeparationIdx)
    branchIDs = str2double(branchIDsString);
elseif ~isempty(commaSeparationIdx)
    numberOfCommas = numel(commaSeparationIdx); branchIDs = [];
    for cComma = 1 : numberOfCommas
           if cComma == numberOfCommas && cComma ~= 1
               branchIDs = [branchIDs,str2double(branchIDsString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
               branchIDs = [branchIDs,str2double(branchIDsString(commaSeparationIdx(cComma):end))];
           elseif cComma == numberOfCommas && cComma == 1
               branchIDs = [str2double(branchIDsString(1:commaSeparationIdx-1)), str2double(branchIDsString(commaSeparationIdx+1:end))];
           elseif cComma == 1
               branchIDs = [branchIDs,str2double(branchIDsString(1:commaSeparationIdx(1)-1))];
           else
               branchIDs = [branchIDs,str2double(branchIDsString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
           end
    end
end


if ~isnan(branchIDs)
    hObject.ForegroundColor = [0 0 0]; hObject.FontAngle = 'normal';
else
    handles.edit_branchIDs.String = 'Branch IDs';
    handles.edit_branchIDs.ForegroundColor = [0.5020 0.5020 0.5020];
    handles.edit_branchIDs.FontAngle = 'italic';
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_branchIDs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_branchIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_customStartPointDistances_Callback(hObject, eventdata, handles)
% hObject    handle to edit_customStartPointDistances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_customStartPointDistances as text
%        str2double(get(hObject,'String')) returns contents of edit_customStartPointDistances as a double

customStartPointDistancesString = handles.edit_customStartPointDistances.String;
commaSeparationIdx = strfind(customStartPointDistancesString,',');
if isempty(commaSeparationIdx)
    customStartPointDistances = str2double(customStartPointDistancesString);
elseif ~isempty(commaSeparationIdx)
    numberOfCommas = numel(commaSeparationIdx); customStartPointDistances = [];
    for cComma = 1 : numberOfCommas

       if cComma == numberOfCommas && cComma ~= 1
           customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
           customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(commaSeparationIdx(cComma):end))];
       elseif cComma == numberOfCommas && cComma == 1
           customStartPointDistances = [str2double(customStartPointDistancesString(1:commaSeparationIdx-1)), str2double(customStartPointDistancesString(commaSeparationIdx+1:end))];
       elseif cComma == 1
           customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(1:commaSeparationIdx(1)-1))];
       else
           customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
       end
    end
end

if ~isnan(customStartPointDistances)
    handles.edit_customStartPointDistances.ForegroundColor = [0 0 0]; 
    handles.edit_customStartPointDistances.FontAngle = 'normal';
else
    handles.edit_customStartPointDistances.String = 'Custom start point distances';
    handles.edit_customStartPointDistances.ForegroundColor = [0.5020 0.5020 0.5020];
    handles.edit_customStartPointDistances.FontAngle = 'italic';
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_customStartPointDistances_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_customStartPointDistances (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_plotGraphs.
function pushbutton_plotGraphs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_plotGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 expDirIndexString = handles.edit_experimentDirectoryIndex.String;
 commaSeparationIdx = strfind(expDirIndexString,',');
 if isempty(commaSeparationIdx)
 %CG: in this case we can expect the user to have entered only one numeric
 %format.
 	expDirIndex = str2double(expDirIndexString);
 elseif ~isempty(commaSeparationIdx)
 %CG: in this case there is more than one numeric value contained within the string.
 %Because str2double returns NaN with multiple numeric values within a
 %single string, it is necessary to iteratively collect each numeric value
 %using commaSeparationIdx.
    numberOfCommas = numel(commaSeparationIdx); expDirIndex = [];
    for cComma = 1 : numberOfCommas
        
 %CG: check in front and behind.
           if cComma == numberOfCommas && cComma ~= 1
               expDirIndex = [expDirIndex,str2double(expDirIndexString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
               expDirIndex = [expDirIndex,str2double(expDirIndexString(commaSeparationIdx(cComma):end))];
           elseif cComma == numberOfCommas && cComma == 1
               expDirIndex = [str2double(expDirIndexString(1:commaSeparationIdx-1)), str2double(expDirIndexString(commaSeparationIdx+1:end))];    
           elseif cComma == 1
               expDirIndex = [expDirIndex,str2double(expDirIndexString(1:commaSeparationIdx(1)-1))];
           else
               expDirIndex = [expDirIndex,str2double(expDirIndexString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
           end
    end
 end
 %CG: the same must now be done for both Branch IDs and custom start
 %distances.
 branchIDsString = handles.edit_branchIDs.String;
 commaSeparationIdx = strfind(branchIDsString,',');
 if isempty(commaSeparationIdx)
 	branchIDs = str2double(branchIDsString);
 elseif ~isempty(commaSeparationIdx)
    numberOfCommas = numel(commaSeparationIdx); branchIDs = [];
    for cComma = 1 : numberOfCommas
           if cComma == numberOfCommas && cComma ~= 1
               branchIDs = [branchIDs,str2double(branchIDsString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
               branchIDs = [branchIDs,str2double(branchIDsString(commaSeparationIdx(cComma):end))];
           elseif cComma == numberOfCommas && cComma == 1
               branchIDs = [str2double(branchIDsString(1:commaSeparationIdx-1)), str2double(branchIDsString(commaSeparationIdx+1:end))];
           elseif cComma == 1
               branchIDs = [branchIDs,str2double(branchIDsString(1:commaSeparationIdx(1)-1))];
           else
               branchIDs = [branchIDs,str2double(branchIDsString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
           end
    end
 end
 
 customStartPointDistancesString = handles.edit_customStartPointDistances.String;
 commaSeparationIdx = strfind(customStartPointDistancesString,',');
 if isempty(commaSeparationIdx)
 	customStartPointDistances = str2double(customStartPointDistancesString);
 elseif ~isempty(commaSeparationIdx)
    numberOfCommas = numel(commaSeparationIdx); customStartPointDistances = [];
    for cComma = 1 : numberOfCommas
        
           if cComma == numberOfCommas && cComma ~= 1
               customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
               customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(commaSeparationIdx(cComma):end))];
           elseif cComma == numberOfCommas && cComma == 1
               customStartPointDistances = [str2double(customStartPointDistancesString(1:commaSeparationIdx-1)), str2double(customStartPointDistancesString(commaSeparationIdx+1:end))];    
           elseif cComma == 1
               customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(1:commaSeparationIdx(1)-1))];
           else
               customStartPointDistances = [customStartPointDistances,str2double(customStartPointDistancesString(commaSeparationIdx(cComma-1)+1:commaSeparationIdx(cComma)-1))];
           end
    end
 end
 
 %CG: Call BranchTracks to make final processing steps and create graphs
 %only if numel are equal for all customization arrays.  
 Dirs = handles.Dirs; keywordStrs = handles.keywordStrs; configStruct = handles.configStruct;
 if (numel(expDirIndex) == numel(branchIDs) && numel(branchIDs) == numel(customStartPointDistances)) &&...
         iscell(Dirs) && iscell(keywordStrs)
     BranchTracks(Dirs, keywordStrs, expDirIndex, branchIDs,...
     customStartPointDistances, configStruct)
%      BranchTracksForDiameterPlots(Dirs, keywordStrs, expDirIndex, branchIDs,...
%          customStartPointDistances, configStruct)
 end

function edit_selectedExperiments_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selectedExperiments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selectedExperiments as text
%        str2double(get(hObject,'String')) returns contents of edit_selectedExperiments as a double

Dirs = handles.Dirs; Sz_Dirs = size(Dirs); NumDirs = Sz_Dirs(1)*Sz_Dirs(2);
checkbox_bulkLoad = handles.checkbox_bulkLoad;

if handles.pushbutton_selectExperiments.Value == 1
    if checkbox_bulkLoad.Value == 0
        text = ''; hObject.Max = NumDirs;
        for cDirIdx = 1 : NumDirs
            cDir = Dirs{cDirIdx};
            if ~isempty(strfind(handles.osType, 'MAC')) 
                Slashes = strfind(cDir, '/'); 
            else 
                Slashes = strfind(cDir, '\');
            end 
%CG: no need to present the whole directory. Just something showing the
%folder name while illustrating that a directory is the captured string. 
            if strcmp(cDir(1),'/') 
                sh_cDir = strcat('/...',cDir(Slashes(end):end));
            else
                sh_cDir = strcat('...',cDir(Slashes(end)+1:end));
            end 
            if NumDirs > 1; text = strcat(text,sh_cDir,'\n'); elseif NumDirs == 1; text = sh_cDir; end
        end
    elseif checkbox_bulkLoad.Value == 1 
        dfcount = 0; text = ''; rDirs = {}; rDir_counter = 0;
        for cDirIdx = 1 : NumDirs
            cDir = Dirs{cDirIdx}; cd(cDir); 
            dirInfo = dir; NumItems = size(dirInfo,1); dfcount = dfcount + NumItems;

            for cItemIdx = 1 : NumItems
                cItem = dirInfo(cItemIdx).name;
                if exist(cItem) == 7 && ~strcmp(cItem,'.') && ~strcmp(cItem,'..') && isempty(strfind(cItem,'ThresholdRecords'))
                    sh_cDir = strcat('/.../',cItem); text = strcat(text,sh_cDir,'\n');
                    rDir_counter = rDir_counter + 1; rDirs{rDir_counter} = strcat(cDir,'/',cItem);
                end
            end
        end
        hObject.Max = dfcount; handles.Dirs = rDirs;
    end
    hObject.ForegroundColor = [0 0 0];
    hObject.FontAngle = 'normal';
    hObject.String = sprintf(text);
elseif handles.pushbutton_selectExperiments.Value == 0
    
end

if isempty(hObject.String)
    handles.cDir = ''; hObject.String = 'Select experiment';
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
else
    if strcmp(hObject.String,' ')
        handles.cDir = ''; hObject.String = 'Select experiment';
        hObject.ForegroundColor = [0.502 0.502 0.502]; 
        hObject.FontAngle = 'italic';
    end
end
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_selectedExperiments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selectedExperiments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_selectExperiments.
function pushbutton_selectExperiments_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectExperiments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Dirs = uipickfiles('prompt','Pick experiments for AARGI');
% Dirs{1}
if iscell(Dirs)
    handles.Dirs = Dirs;
    edit_selectedExperiments_Callback(handles.edit_selectedExperiments, eventdata, handles);
end
guidata(hObject, handles);

% --- Executes on button press in checkbox_bulkLoad.
function checkbox_bulkLoad_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_bulkLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_bulkLoad



function edit_keywords_Callback(hObject, eventdata, handles)
% hObject    handle to edit_keywords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_keywords as text
%        str2double(get(hObject,'String')) returns contents of edit_keywords as a double

keywordStrs = handles.keywordStrs;

kw_str = hObject.String; space_idces = strfind(kw_str, ' '); 
if iscell(kw_str); kw_str = kw_str{1}; end
if iscell(space_idces); space_idces = space_idces{1}; end
%CG: for some reason matlab can change kw_str and/or space_idces to a cell.
if ~strcmp(kw_str, 'Enter key word(s) here') && ~isempty(kw_str)

%CG: if the user puts a space after the last keyword, the number of
%keywords will be: NumStrs = numel(space_idces);
    if strcmp(kw_str(end),' '); NumStrs = numel(space_idces);
    elseif ~strcmp(kw_str(end),' '); NumStrs = numel(space_idces) + 1;
    end
    
%CG: any spaces before the first keyword should be cleared
    if strcmp(kw_str(1),' ')
        aa = diff(space_idces);
        if ~isempty(aa)
            aa = [1,aa]; g1Idx = find(aa > 1); 
            if isempty(g1Idx) && sum(g1Idx) == numel(g1Idx)
                kw_str = kw_str(1:numel(g1Idx));
%CG: if isempty(g1Idx) && sum(g1Idx) == numel(g1Idx), then this probably means 
%only 1 keyword was entered with at least 1 space before it. 
            else
                kw_str = kw_str(1:g1Idx);
            end
        end
    end

    if ~isempty(space_idces);
        for cstrIdx = 1 : NumStrs
            if cstrIdx == 1
                cKW = kw_str(1:space_idces(cstrIdx)-1);
            else
                if numel(space_idces) < cstrIdx
                    cKW = kw_str(space_idces(cstrIdx-1)+1:end);
                else
                    cKW = kw_str(space_idces(cstrIdx-1)+1:space_idces(cstrIdx)-1);
                end
            end
            keywordStrs{cstrIdx} = cKW;
%CG: the contents of suffixStrs should mirror the contents of the key
%word(s) editable text box.
        end
    elseif isempty(space_idces);
%CG: if isempty(space_idces), this probably means that only a single keyword 
%has been entered and no space has been placed after it.
            keywordStrs = {};
            keywordStrs{1} = kw_str;
    end
    if size(keywordStrs,2) == 1
        disp(strcat('"',keywordStrs{1},'"'))
    elseif size(keywordStrs,2) == 2
        disp(strcat('"',keywordStrs{1}, '" and "', keywordStrs{2}, '"'))
    end
    
    hObject.ForegroundColor = [0 0 0];
    hObject.FontAngle = 'normal';
    
elseif strcmp(kw_str, 'Enter keyword(s) here')
    hObject.ForegroundColor = [0.502 0.502 0.502];
    hObject.FontAngle = 'italic';
    
end

if isempty(kw_str)
    hObject.String = 'Enter keyword(s) here';
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
end
handles.keywordStrs = keywordStrs;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_keywords_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_keywords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_configureGraphs.
function pushbutton_configureGraphs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_configureGraphs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
configStruct = handles.configStruct;

[configGraphHandles] = configureGraphs(configStruct);
handles.configStruct = configGraphHandles.configStruct;
delete(configGraphHandles.figure1)

guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

delete(hObject);

