function varargout = AARGII(varargin)

% Outline
% AARGII uses Mathworks' gui-building feature 'guide' to generate a GUI also
% named AARGII. The AARGII gui enables dendrites to be traced so that all
% ROIs are spatially referenced to a single user-defined start point. The
% next step is to verify each detection manually in the Review Data stage.
% A preferred playlist or audiobook is recommended for this stage.

% Author: Charlie J Gilbride
% version: 1.0.0

% AARGII MATLAB code for AARGII.fig
%      AARGII, by itself, creates a new AARGII or raises the existing
%      singleton*.
%
%      H = AARGII returns the handle to a new AARGII or the handle to
%      the existing singleton*.
%
%      AARGII('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AARGII.M with the given input arguments.
%
%      AARGII('Property','Value',...) creates a new AARGII or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AARGII_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AARGII_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AARGII

% Last Modified by GUIDE v2.5 26-Sep-2018 20:46:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AARGII_OpeningFcn, ...
                   'gui_OutputFcn',  @AARGII_OutputFcn, ...
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


% --- Executes just before AARGII is made visible.
function AARGII_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AARGII (see VARARGIN)

% Choose default command line output for AARGII
osType = computer;
handles.osType = osType;

functionName = 'AARGII.m'; functionDir = which(functionName);
functionDir = functionDir(1:end-length(functionName));
addpath(genpath(functionDir)); 

if ~isempty(strfind(handles.osType, 'MAC')) 
    Slashes = strfind(functionDir, '/'); 
else 
    Slashes = strfind(functionDir, '\');
end

functionDir = functionDir(1:Slashes(end-1)); 
functionDir = strcat(functionDir, 'commonFunctions');
addpath(genpath(functionDir));

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

suffixStrs = {}; handles.suffixStrs = suffixStrs;
keywordStrs = {}; handles.keywordStrs = keywordStrs;
cDir = ''; handles.cDir = cDir; handles.frameRate = '';
%CG: default status of the resize checkbox should be that it is checked.
%Currently, data files are resized by default to reduce computational load.
handles.checkbox_NewBuild.Value = 1;
handles.currentExpSetByTextMod = '';

%CG: Setup default state for the modify sub-section
handles.radiobutton_AddStartPntsorLines.Enable = 'off';
handles.radiobutton_AddWhiteLines.Enable = 'off';
handles.radiobutton_RemROIConx.Enable = 'off';
handles.radiobutton_Skip2ROISelection.Enable = 'off';
handles.uipanel_modify.ForegroundColor = [0.502 0.502 0.502];


%*************testing only
% handles.edit_keywords.String = 'a_baseline b_CTX';
% cDir = '/Volumes/No2/MATLAB/AARG_Data/ReAnalysis/conotoxin/testing/generalTesting/CG0706173';
% handles.cDir = cDir; Slashes = strfind(cDir, '/');
% text = strcat('/...',cDir(Slashes(end):end));
% [~,CellName,~] = fileparts(cDir);
% handles.edit_cExp.String = CellName;
% handles.edit_cExp.ForegroundColor = [0 0 0];
% handles.edit_selexpDisp.String = sprintf(text);
% pause(1)
% handles.edit_MaxDist.String = '5';
% handles.edit_shROIsFrCon.String = 'b_CTX';
% handles.edit_Resolution.String = '6.25';
% handles.keywordStrs = {'a_baseline','b_CTX'};
%*************testing only

% Update handles structure
set( findall( hObject, '-property', 'Units' ), 'Units', 'Normalized' )
set( hObject, 'Resize','on' ) 
guidata(hObject, handles);

% UIWAIT makes AARGII wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AARGII_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_keywords_Callback(hObject, eventdata, handles)
% hObject    handle to edit_keywords (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_keywords as text
%        str2double(get(hObject,'String')) returns contents of edit_keywords as a double
keywordStrs = handles.keywordStrs;

kw_str = hObject.String; space_idces = strfind(kw_str, ' '); 

%CG: contents of the editable text box must be something else other than
%the default contents.
if ~strcmp(kw_str, 'Enter keyword(s) here') && ~isempty(kw_str)

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
        keywordStrs = {};
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
%CG: the contents of suffixStrs should mirror the contents of the
%keyword(s) editable text box.
        end
    elseif isempty(space_idces);
%CG: if isempty(space_idces), this probably means that only a single keyword 
%has been entered and no space has been placed after it.
        keywordStrs{1} = kw_str;
    end
    
    hObject.ForegroundColor = [0 0 0];
    hObject.FontAngle = 'normal';
%     if size(keywordStrs,2) == 1
%         disp(strcat('"',keywordStrs{1},'"'))
%     elseif size(keywordStrs,2) == 2
%         disp(strcat('"',keywordStrs{1}, '" and "', keywordStrs{2}, '"'))
%     end
    
    if ~isempty(hObject.String) && ~strcmp(hObject.String,'Enter keyword(s) here')
        handles.edit_shROIsFrCon.String = keywordStrs{end};
    end
elseif strcmp(kw_str, 'Enter keyword(s) here')
    hObject.ForegroundColor = [0.502 0.502 0.502];
    hObject.FontAngle = 'italic';
    
end
if isempty(kw_str)
    hObject.String = 'Enter keyword(s) here';
    keywordStrs = {};
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
    handles.edit_shROIsFrCon.String = 'no keyword(s) specified';
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



function edit_selexpDisp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_selexpDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_selexpDisp as text
%        str2double(get(hObject,'String')) returns contents of edit_selexpDisp as a double

cDir = handles.cDir; 

if ~isempty(strfind(handles.osType, 'MAC')) 
    Slashes = strfind(cDir, '/'); 
else 
    Slashes = strfind(cDir, '\');
end 

oldCellName = handles.edit_cExp.String; newText = hObject.String; 
if isempty(handles.currentExpSetByTextMod)
    handles.currentExpSetByTextMod = 'SetByText';
end

if exist(cDir) == 7
    if strcmp(cDir(1),'/') 
        cDirText = strcat('/...',cDir(Slashes(end):end));
    else
        cDirText = strcat('...',cDir(Slashes(end)+1:end));
    end
    if isempty(oldCellName) || strcmp(newText,cDirText) ||...
            strcmp(handles.currentExpSetByTextMod, 'SetByPushButton')
        [~,CellName,~] = fileparts(cDir);
        handles.edit_cExp.String = CellName;
        hObject.ForegroundColor = [0 0 0];
        hObject.FontAngle = 'normal';
        hObject.String = sprintf(cDirText);
    elseif ~strcmp(newText,cDirText)
        %newText = get(hObject,'String');
        if strcmp(cDir(1),'/') 
            moreSlashes = strfind(get(hObject,'String'),'/');
        else 
            moreSlashes = strfind(get(hObject,'String'),'\');
        end
        if ~isempty(moreSlashes)
            newCellName = newText(moreSlashes(end)+1:end);
            newDir = strcat(cDir(1:Slashes(end)),newCellName);
        
            if exist(newDir) == 7
%CG: if the modification produces a recognizible directory, then update
%variables.
                handles.cDir = newDir; handles.edit_cExp.String = newCellName;
                hObject.String = sprintf(newText);
            else
%CG: if no directory is recognized, use callback function: pushbutton_selexp_Callback
                handles.cDir = ''; handles.edit_cExp.String = 'none selected';
                pushbutton_selexp_Callback(handles.pushbutton_selexp, eventdata, handles);
            end
        else
            handles.cDir = ''; handles.edit_cExp.String = 'none selected';
%             pushbutton_selexp_Callback(handles.pushbutton_selexp, eventdata, handles);
        end
    end
    if strcmp(handles.currentExpSetByTextMod, 'SetByPushButton')
        handles.currentExpSetByTextMod = '';
    end
    handles.currentExpSetByTextMod = 'No';
    
end

if isempty(hObject.String)
    handles.cDir = ''; hObject.String = 'Select experiment for analysis';
    pause(0.25)
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
end
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_selexpDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selexpDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton_selexp.
function pushbutton_selexp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cDir = uipickfiles('prompt','Pick ONE experiment for AARGII');
Sz_cDir = size(cDir); 
if iscell(cDir) 
    if ~isempty(cDir) && Sz_cDir(2) == 1
        handles.cDir = cDir{1}; handles.currentExpSetByTextMod = 'SetByPushButton';
        edit_selexpDisp_Callback(handles.edit_selexpDisp, eventdata, handles);
        text_status.String = 'Pick only ONE experiment at a time';
        text_status.ForegroundColor = [0.502 0.502 0.502];
        text_status.FontSize = 12;
    elseif ~isempty(cDir) && Sz_cDir(2) > 1 
        text_status = handles.text_status;
        text_status.String = 'Pick only ONE experiment at a time';
        text_status.ForegroundColor = [1 0 0];
        text_status.FontSize = 16;
    else
        text_status = handles.text_status;
        text_status.String = 'Idle';
        text_status.ForegroundColor = [0.502 0.502 0.502];
        text_status.FontSize = 12;
    end

    handles.tex_status = text_status;
    guidata(hObject, handles);
end



function edit_cExp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_cExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_cExp as text
%        str2double(get(hObject,'String')) returns contents of edit_cExp as a double


% --- Executes during object creation, after setting all properties.
function edit_cExp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_cExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_shROIsFrCon_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shROIsFrCon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shROIsFrCon as text
%        str2double(get(hObject,'String')) returns contents of edit_shROIsFrCon as a double


% --- Executes during object creation, after setting all properties.
function edit_shROIsFrCon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shROIsFrCon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Resolution_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Resolution as text
%        str2double(get(hObject,'String')) returns contents of edit_Resolution as a double

if isempty(hObject.String)
    handles.cDir = ''; hObject.String = 'Enter value';
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
else
    if strcmp(hObject.String,' ')
        handles.cDir = ''; hObject.String = 'Enter value';
        hObject.ForegroundColor = [0.502 0.502 0.502]; 
        hObject.FontAngle = 'italic';
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_Resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MaxDist_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxDist as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxDist as a double

if isempty(hObject.String)
    handles.cDir = ''; hObject.String = 'Enter value';
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
else
    if strcmp(hObject.String,' ')
        handles.cDir = ''; hObject.String = 'Enter value';
        hObject.ForegroundColor = [0.502 0.502 0.502]; 
        hObject.FontAngle = 'italic';
    end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_MaxDist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_NewBuild.
function checkbox_NewBuild_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_NewBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_NewBuild
 
if handles.checkbox_NewBuild.Value == 1
    handles.checkbox_ModOldBuild.Value = 0;
%     handles.checkbox_traceDataOnly.Value = 0;
    handles.radiobutton_AddStartPntsorLines.Value = 0;
    handles.radiobutton_AddStartPntsorLines.Enable = 'off';
    handles.radiobutton_AddWhiteLines.Value = 0;
    handles.radiobutton_AddWhiteLines.Enable = 'off';
    handles.radiobutton_RemROIConx.Value = 0;
    handles.radiobutton_RemROIConx.Enable = 'off';
    handles.radiobutton_Skip2ROISelection.Value = 0;
    handles.radiobutton_Skip2ROISelection.Enable = 'off';
    handles.uipanel_modify.ForegroundColor = [0.502 0.502 0.502];
% elseif handles.checkbox_NewBuild.Value == 0
%     handles.radiobutton_AddStartPntsorLines.Enable = 'on';
%     handles.radiobutton_AddWhiteLines.Enable = 'on';
%     handles.radiobutton_RemROIConx.Enable = 'on';
%     handles.radiobutton_Skip2ROISelection.Enable = 'on';
%     handles.uipanel_modify.ForegroundColor = [0 0 0];
end
guidata(hObject, handles);    

% --- Executes on button press in checkbox_ModOldBuild.
function checkbox_ModOldBuild_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ModOldBuild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ModOldBuild
if handles.checkbox_ModOldBuild.Value == 1
    handles.checkbox_NewBuild.Value = 0;
%     handles.checkbox_traceDataOnly.Value = 0;
    handles.radiobutton_AddStartPntsorLines.Enable = 'on';
    handles.radiobutton_AddWhiteLines.Enable = 'on';
    handles.radiobutton_RemROIConx.Enable = 'on';
    handles.radiobutton_Skip2ROISelection.Enable = 'on';
    handles.uipanel_modify.ForegroundColor = [0 0 0];
elseif handles.checkbox_ModOldBuild.Value == 0
    handles.radiobutton_AddStartPntsorLines.Value = 0;
    handles.radiobutton_AddStartPntsorLines.Enable = 'off';
    handles.radiobutton_AddWhiteLines.Value = 0;
    handles.radiobutton_AddWhiteLines.Enable = 'off';
    handles.radiobutton_RemROIConx.Value = 0;
    handles.radiobutton_RemROIConx.Enable = 'off';
    handles.radiobutton_Skip2ROISelection.Value = 0;
    handles.radiobutton_Skip2ROISelection.Enable = 'off';
    handles.uipanel_modify.ForegroundColor = [0.502 0.502 0.502];
end
guidata(hObject, handles);    

% --- Executes on button press in pushbutton_Launch.
function pushbutton_Launch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Launch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
AbortLaunch = 0; cDir = handles.cDir;
if strcmp(handles.edit_cExp.String, 'none selected')
    handles.text_status.String = 'Select experiment to analyse';
    handles.text_status.ForegroundColor = [1 0 0];
    handles.text_status.FontSize = 16; AbortLaunch = 1;
end
if strcmp(handles.edit_shROIsFrCon.String, 'no keyword(s) specified')
    handles.text_status.String = 'Specify key word(s)';
    handles.text_status.ForegroundColor = [1 0 0];
    handles.text_status.FontSize = 16; AbortLaunch = 1;
end
if strcmp(handles.edit_Resolution.String, 'Enter value')
    handles.text_status.String = 'Set resolution';
    handles.text_status.ForegroundColor = [1 0 0];
    handles.text_status.FontSize = 16; AbortLaunch = 1;
end
if strcmp(handles.edit_MaxDist.String, 'Enter value')
    handles.text_status.String = 'Set maximum distance ROIs can be from any dendrite';
    handles.text_status.ForegroundColor = [1 0 0];
    handles.text_status.FontSize = 16; AbortLaunch = 1;
end
if AbortLaunch == 0
    
    pause(1)
    handles.text_status.String = 'Initiating...';
    handles.text_status.ForegroundColor = [0.502 0.502 0.502];
    handles.text_status.FontSize = 12; 
    
    if handles.checkbox_NewBuild.Value == 1
        cd(cDir); dirInfo = dir; NumItems = size(dirInfo,1); cItemIdx = 0; TgtFound = 0;
        GoStr = 'Go!';
        while cItemIdx < NumItems && TgtFound == 0
            cItemIdx = cItemIdx + 1;
            cItem = dirInfo(cItemIdx).name;
            if ~isempty(strfind(cItem,'CellBones'))                
                choice = questdlg('Overwrite existing cell skeleton data for current experiment? ', ...
                    'Overwriting',...
                    'Overwrite', 'Cancel','Cancel');
                switch choice                        
                    case 'Cancel'
                        GoStr = 'NoGo';
                end
                TgtFound = 1;
            end
        end
        
        if strcmp(GoStr,'Go!')
            md = str2double(handles.edit_MaxDist.String);
            keywordStrs = handles.keywordStrs;
            suffixStrs = keywordStrs(end);
            rr = str2double(handles.edit_Resolution.String);
            if TgtFound == 1; overwriteStatus = 'On'; elseif TgtFound == 0; overwriteStatus = 'Off'; end
%CG: no need to overwrite if there is nothing to overwrite.
            CellBones('datadir',cDir,'suffixStrs', suffixStrs,'MaxDistance',md,...
                'Resolution',rr,'Overwrite',overwriteStatus)
        end
    elseif handles.checkbox_ModOldBuild.Value == 1
        cd(handles.cDir); dirInfo = dir; NumItems = size(dirInfo,1); cItemIdx = 0; TgtFound = 0;
        while cItemIdx < NumItems && TgtFound == 0
            cItemIdx = cItemIdx + 1;
            cItem = dirInfo(cItemIdx).name;
            if ~isempty(strfind(cItem,'CellBones'))                
                TgtFound = 1;
            end
        end
        if TgtFound == 0
            handles.text_status.String = 'No previous build detected for current experiment. Start New Build';
            handles.text_status.ForegroundColor = [1 0 0];
            handles.text_status.FontSize = 16; 
        elseif TgtFound == 1
            GoStr = 'Go!';
            if strcmp(GoStr,'Go!')
                md = str2double(handles.edit_MaxDist.String);
                keywordStrs = handles.keywordStrs;
                suffixStrs = keywordStrs(end);
                rr = str2double(handles.edit_Resolution.String);

                if handles.radiobutton_AddStartPntsorLines.Value == 1
%                     pause(2); handles.text_status.String = 'Initiating...';
%                     handles.text_status.ForegroundColor = [0.502 0.502 0.502];
%                     handles.text_status.FontSize = 12; 
                    CellBones('datadir',cDir,'suffixStrs', suffixStrs,'MaxDistance',md,...
                        'Resolution',rr, 'Overwrite','Off','AddSPorWL','On')
                elseif handles.radiobutton_RemROIConx.Value == 1
%                     pause(2); handles.text_status.String = 'Initiating...';
%                     handles.text_status.ForegroundColor = [0.502 0.502 0.502];
%                     handles.text_status.FontSize = 12; 
                    CellBones('datadir',cDir,'suffixStrs', suffixStrs,'MaxDistance',md,...
                        'Resolution',rr,'Overwrite','Off','DeleteGreenLines','On')
                elseif handles.radiobutton_Skip2ROISelection.Value == 1
%                     pause(2); handles.text_status.String = 'Initiating...';
%                     handles.text_status.ForegroundColor = [0.502 0.502 0.502];
%                     handles.text_status.FontSize = 12; 
                    CellBones('datadir',cDir,'suffixStrs', suffixStrs,'MaxDistance',md,...
                        'Resolution',rr,'Overwrite','Off','SkipDistanceMeasurement','On')
                end
                
                    
            end
                
        end
    else
%CG: don't do anything.             
        
    end
    pause(1); handles.text_status.String = 'Idle';
    handles.text_status.ForegroundColor = [0.502 0.502 0.502];
    handles.text_status.FontSize = 12;
end
guidata(hObject,handles)

% --- Executes on button press in pushbutton_Help.
function pushbutton_Help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

open('HelpAARGII.pdf')

%--------8. Cell skeleton------------

%Resize checkbox: by default, this box is checked. It should always be kept
%checked, unless the user clear understands why (s)he would have it
%unchecked. When events are detected (see AARGI), the data from the tif
%stacks is resized such that pixels are binned by a factor of three. The
%binning factor is (at the time of writing) hardcoded. If the checkbox is
%unchecked, then the function that manages dendrite tracing and distance
%measurement will not "remember" that the data has been resized. It is
%equivalent to using the wrong pixel resolution and will result in distance
%measurements being (very) wrong. If the pixel resolution before resizing
%is x, the user could enter the solution to x/3 into the Resolution text
%box and uncheck the Resize checkbox and get the correct answer (but why
%would that be a more sensible option?)


% --- Executes on button press in radiobutton_AddStartPntsorLines.
function radiobutton_AddStartPntsorLines_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_AddStartPntsorLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_AddStartPntsorLines

% --- Executes on button press in radiobutton_RemROIConx.
function radiobutton_RemROIConx_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_RemROIConx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_RemROIConx
if handles.radiobutton_RemROIConx.Value == 1
    handles.radiobutton_AddStartPntsorLines.Value = 0;
    handles.radiobutton_AddStartPntsorLines.Enable = 'Off';
elseif handles.radiobutton_RemROIConx.Value == 0
    handles.radiobutton_AddStartPntsorLines.Enable = 'On';
end

guidata(hObject,handles)

% --- Executes on button press in radiobutton_Skip2ROISelection.
function radiobutton_Skip2ROISelection_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_Skip2ROISelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_Skip2ROISelection

if handles.radiobutton_Skip2ROISelection.Value == 1
    handles.radiobutton_AddStartPntsorLines.Value = 0;
    handles.radiobutton_AddStartPntsorLines.Enable = 'Off';
    handles.radiobutton_RemROIConx.Value = 0;
    handles.radiobutton_RemROIConx.Enable = 'Off';
elseif handles.radiobutton_Skip2ROISelection.Value == 0
    handles.radiobutton_AddStartPntsorLines.Enable = 'On';
    handles.radiobutton_RemROIConx.Enable = 'On';
end

guidata(hObject,handles)


% --- Executes on button press in pushbutton_reviewLaunch.
function pushbutton_reviewLaunch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reviewLaunch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frameRate = handles.frameRate; validFrameRate = 0;

if ~isempty(frameRate) 
    if ~strcmp(frameRate, 'Enter frame rate (Hz)') && isnumeric(str2double(frameRate))
        validFrameRate = 1;
    end
end
        
text_status = handles.text_status;
if validFrameRate == 1

    text_status.String = 'Launching...';
    pause(1)
    frameRate = str2double(frameRate); suffixStrs = handles.keywordStrs;
    cDir = handles.cDir; CellName = handles.edit_cExp.String;
    try 
        numberOfConditions = size(suffixStrs,2); cConditionIdx = 0; dataIncomplete = 0;
        while cConditionIdx < numberOfConditions && dataIncomplete == 0
            cConditionIdx = cConditionIdx + 1;
%             breaktry
            peakAmplitudeData = load(strcat(CellName,suffixStrs{cConditionIdx},'_MeasuresRAW.mat'));
            undostring = peakAmplitudeData.undostring; parameters = peakAmplitudeData.parameters;
            CA_Measure = peakAmplitudeData.CA_Measure; Sz_CA_Measure = size(CA_Measure);
            totaltraces = 0;
            for cROI = 1 : Sz_CA_Measure(1);
                if ~isempty(CA_Measure{cROI,11})
                    CA_ew = CA_Measure{cROI,11}; 
                    numpotseg = size(CA_ew,2); seg_indices = [];
                    for cRow = 1 : numpotseg
                        if ~isempty(CA_ew{cRow}); totaltraces = totaltraces + 1; end
                    end
                end
            end
            if numel(undostring) < totaltraces
                dataIncomplete = 1;
            end
        end
        breakTry
        ReviewFit(cDir,suffixStrs,text_status,cConditionIdx,frameRate)
    catch
        conditionsToShow = (1:1:numberOfConditions); dataIncomplete = 1;
%CG: GetPeakAmplitudes analyses ROIs which are sorted according to event
%frequency. ROIs with the most events are listed at the top. 
%CG: conditionsToShow specifies for which conditions ShowTrace will apply.
%example: if ShowTrace = [1]; and conditionsToShow =
%(1:1:numberOfConditions), then GetPeakAmplitudes will show the trace from
%the most active ROI for each condition of the current experiment. If ShowTrace 
%is changed to ShowTrace = [2], then the trace of the second-most active ROI
%will be displayed for each condition. 
        pause(1)
        text_status.String = 'Identifying peak amplitude locations...please wait a few moments'; 
        Dirs = {}; Dirs{1} = cDir; 
        GetPeakAmplitudes('Dirs',Dirs,'suffixStrs',suffixStrs,'ShowTrace',[1],...
            'ConditionSpec',conditionsToShow);
        pause(1)
        text_status.String = 'Peaks identified...proceeding to review';
        pause(5)
        ReviewFit(cDir,suffixStrs,text_status,cConditionIdx,frameRate)
    end
    if dataIncomplete == 1
        text_status.String = 'Idle'; text_status.ForegroundColor = [0.502 0.502 0.502]; text_status.FontSize = 12;
    elseif dataIncomplete == 0
        text_status.String = strcat('"',CellName,'" has already been completely analysed');
        text_status.ForegroundColor = [0 0 0]; text_status.FontSize = 18;
        pause(5)
        text_status.String = 'Idle'; text_status.ForegroundColor = [0.502 0.502 0.502]; text_status.FontSize = 12;
    end
elseif validFrameRate == 0
    text_status.String = 'A valid frame rate must be entered before traces can be reviewed!';
    text_status.ForegroundColor = [0 0 0]; text_status.FontSize = 18;
end

% --- Executes on button press in pushbutton_graphSetup.
function pushbutton_graphSetup_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_graphSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selectExpText = handles.edit_cExp.String; keywordStrs = handles.keywordStrs;

if ~ischar(keywordStrs) && ~strcmp(selectExpText, 'Select experiment for analysis')
    graphAARG(keywordStrs)
end



function edit_frameRate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frameRate as text
%        str2double(get(hObject,'String')) returns contents of edit_frameRate as a double

frameRate = hObject.String; text_status = handles.text_status;

if ~strcmp(frameRate, 'Enter frame rate (Hz)') && ~isempty(frameRate)
    
    hObject.ForegroundColor = [0 0 0];
    hObject.FontAngle = 'normal';
    
    handles.frameRate = frameRate;
    if ~strcmp(text_status.String, 'Idle')
        text_status.String = 'Idle'; 
        text_status.ForegroundColor = [0.502 0.502 0.502]; text_status.FontSize = 12;
        pause(1)
    end
    
elseif isempty(frameRate)
    
    hObject.String = 'Enter frame rate (Hz)';
    hObject.ForegroundColor = [0.502 0.502 0.502];
    hObject.FontAngle = 'italic';
    
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_frameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

