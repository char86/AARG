function varargout = configureGraphs(varargin)

% Author: Charlie J Gilbride
% version: 1.0.0

% CONFIGUREGRAPHS MATLAB code for configureGraphs.fig
%      CONFIGUREGRAPHS, by itself, creates a new CONFIGUREGRAPHS or raises the existing
%      singleton*.
%
%      H = CONFIGUREGRAPHS returns the handle to a new CONFIGUREGRAPHS or the handle to
%      the existing singleton*.
%
%      CONFIGUREGRAPHS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONFIGUREGRAPHS.M with the given input arguments.
%
%      CONFIGUREGRAPHS('Property','Value',...) creates a new CONFIGUREGRAPHS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before configureGraphs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to configureGraphs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help configureGraphs

% Last Modified by GUIDE v2.5 08-Apr-2019 20:13:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @configureGraphs_OpeningFcn, ...
                   'gui_OutputFcn',  @configureGraphs_OutputFcn, ...
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


% --- Executes just before configureGraphs is made visible.
function configureGraphs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to configureGraphs (see VARARGIN)

% Choose default command line output for configureGraphs
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

%CG: first argument is expected to be a structure provided by the calling
%function graphAARG.
configStruct = varargin{1};
handles.configStruct = configStruct;
handles.configStructOldVersion = varargin{1};

applyButtonClicked = 0; handles.applyButtonClicked = applyButtonClicked;
cancelButtonClicked = 0; handles.cancelButtonClicked = cancelButtonClicked;

%CG: setup GUI according to configStruct. These will be default settings or
%will contain modifications set by the user on a previous configureGraphs
%call.
%CG: for the histograms...
handles.checkbox_dontDisplayHistograms.Value = configStruct.dontDisplayHistograms;
handles.checkbox_displayNValues.Value = configStruct.displayNValues;
handles.edit_xMinHistogram.String = num2str(configStruct.xAxisMinHistogram);
if isnan(configStruct.xAxisMaxHistogram)
    handles.edit_xMaxHistogram.String = '?';
elseif isnumeric(configStruct.xAxisMaxHistogram)
    handles.edit_xMaxHistogram.String = num2str(configStruct.xAxisMaxHistogram);
end
handles.edit_yMinHistogram.String = num2str(configStruct.yAxisMinHistogram);
if isnan(configStruct.yAxisMaxHistogram)
    handles.edit_yMaxHistogram.String = '?';
elseif isnumeric(configStruct.yAxisMaxHistogram)
    handles.edit_yMaxHistogram.String = num2str(configStruct.yAxisMaxHistogram);
end
handles.edit_binWidth.String = configStruct.binWidth;
%CG: for the median plots...
handles.checkbox_dontDisplayMedianPlots.Value = configStruct.dontDisplayMedianPlots;
handles.checkbox_displayTTestResult.Value = configStruct.displayTTestResult;
handles.edit_xMinMedianPlot.String = num2str(configStruct.xAxisMinMedian);
if isnan(configStruct.xAxisMaxMedian)
    handles.edit_xMaxMedianPlot.String = '?';
elseif isnumeric(configStruct.xAxisMaxMedian)
    handles.edit_xMaxMedianPlot.String = num2str(configStruct.xAxisMaxMedian);
end
handles.edit_yMinMedianPlot.String = num2str(configStruct.yAxisMinMedian);
if isnan(configStruct.yAxisMaxMedian)
    handles.edit_yMaxMedianPlot.String = '?';
elseif isnumeric(configStruct.yAxisMaxMedian)
    handles.edit_yMaxMedianPlot.String = num2str(configStruct.yAxisMaxMedian);
end
%CG: for all graphs...
handles.checkbox_hideTickLabels.Value = configStruct.hideUnnecessaryTickLabels;
handles.radiobutton_displayPeakAmplitude.Value = configStruct.displayPeakAmplitude;
handles.radiobutton_displayDecay.Value = configStruct.displayDecayTimeConstant;
handles.edit_lengthOfDendritePerGraph.String = num2str(configStruct.lengthOfDendritePerGraph);

if isnan(configStruct.maximumDistance)
    handles.edit_maximumDistance.String = 'maximum distance';
    handles.edit_maximumDistance.ForegroundColor = [0.502 0.502 0.502];
    handles.edit_maximumDistance.FontAngle = 'italic';
elseif ~isnan(configStruct.maximumDistance)
    handles.edit_maximumDistance.String = num2str(configStruct.maximumDistance);
    handles.edit_maximumDistance.ForegroundColor = [0 0 0];
    handles.edit_maximumDistance.FontAngle = 'normal';
end

set(findall(hObject, '-property', 'Units'), 'Units', 'Normalized')
set(hObject, 'Resize','on')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes configureGraphs wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = configureGraphs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

applyButtonClicked = handles.applyButtonClicked;
cancelButtonClicked = handles.cancelButtonClicked;
if applyButtonClicked == 1 || cancelButtonClicked == 1
    varargout{1} = handles;
end

guidata(hObject, handles); 
% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.configStruct = handles.configStructOldVersion;
%CG: discard all changes to configStruct and then end configureGraphs.
handles.cancelButtonClicked = 1;
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.applyButtonClicked = 1;
guidata(hObject, handles);
uiresume(handles.figure1);

% --- Executes on button press in checkbox_dontDisplayMedianPlots.
function checkbox_dontDisplayMedianPlots_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_dontDisplayMedianPlots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_dontDisplayMedianPlots

configStruct = handles.configStruct;
if hObject.Value == 1
    
%CG: don't allow no graphs to be displayed. Why would the user configure
%the settings for no graphs to be displayed?
    if handles.checkbox_dontDisplayHistograms.Value == 1
        hObject.Value = 0;
    elseif handles.checkbox_dontDisplayHistograms.Value == 0
        configStruct.dontDisplayMedianPlots = 1;
    end
elseif hObject.Value == 0
    configStruct.dontDisplayMedianPlots = 0;
end
handles.configStruct = configStruct;

guidata(hObject, handles);

% --- Executes on button press in checkbox_displayTTestResult.
function checkbox_displayTTestResult_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayTTestResult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayTTestResult

configStruct = handles.configStruct;
configStruct.displayTTestResult = hObject.Value;
handles.configStruct = configStruct;

guidata(hObject, handles);

function edit_xMinMedianPlot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xMinMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xMinMedianPlot as text
%        str2double(get(hObject,'String')) returns contents of edit_xMinMedianPlot as a double

configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.xAxisMinMedian = str2double(get(hObject,'String'));
else
    configStruct.xAxisMinMedian = 0;
end
handles.configStruct = configStruct;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_xMinMedianPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xMinMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yMinMedianPlot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yMinMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yMinMedianPlot as text
%        str2double(get(hObject,'String')) returns contents of edit_yMinMedianPlot as a double
configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.yAxisMinMedian = str2double(get(hObject,'String'));
else
    configStruct.yAxisMinMedian = 0;
end
handles.configStruct = configStruct;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_yMinMedianPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yMinMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xMaxMedianPlot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xMaxMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xMaxMedianPlot as text
%        str2double(get(hObject,'String')) returns contents of edit_xMaxMedianPlot as a double

configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.xAxisMaxMedian = str2double(get(hObject,'String'));
else
    configStruct.xAxisMaxMedian = '?';
end
handles.configStruct = configStruct;
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edit_xMaxMedianPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xMaxMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yMaxMedianPlot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yMaxMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yMaxMedianPlot as text
%        str2double(get(hObject,'String')) returns contents of edit_yMaxMedianPlot as a double
configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.yAxisMaxMedian = str2double(get(hObject,'String'));
else
    configStruct.yAxisMaxMedian = '?';
end
handles.configStruct = configStruct;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_yMaxMedianPlot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yMaxMedianPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_dontDisplayHistograms.
function checkbox_dontDisplayHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_dontDisplayHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_dontDisplayHistograms

configStruct = handles.configStruct;
if hObject.Value == 1
    
%CG: don't allow no graphs to be displayed. Why would the user configure
%the settings for no graphs to be displayed?
    if handles.checkbox_dontDisplayMedianPlots.Value == 1
        hObject.Value = 0;
    elseif handles.checkbox_dontDisplayMedianPlots.Value == 0
        configStruct.dontDisplayHistograms = 1;
    end
elseif hObject.Value == 0
    configStruct.dontDisplayHistograms = 0;
end
handles.configStruct = configStruct;

guidata(hObject, handles);

% --- Executes on button press in checkbox_displayNValues.
function checkbox_displayNValues_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayNValues (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayNValues

configStruct = handles.configStruct;
configStruct.displayNValues = hObject.Value;
handles.configStruct = configStruct;

guidata(hObject, handles);

function edit_xMinHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xMinHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xMinHistogram as text
%        str2double(get(hObject,'String')) returns contents of edit_xMinHistogram as a double

configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.xAxisMinHistogram = str2double(get(hObject,'String'));
else
    configStruct.xAxisMinHistogram = 0;
end
handles.configStruct = configStruct;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_xMinHistogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xMinHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yMinHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yMinHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yMinHistogram as text
%        str2double(get(hObject,'String')) returns contents of edit_yMinHistogram as a double

configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.yAxisMinHistogram = str2double(get(hObject,'String'));
else
    configStruct.yAxisMinHistogram = '?';
end
handles.configStruct = configStruct;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_yMinHistogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yMinHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_xMaxHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to edit_xMaxHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_xMaxHistogram as text
%        str2double(get(hObject,'String')) returns contents of edit_xMaxHistogram as a double

configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.xAxisMaxHistogram = str2double(get(hObject,'String'));
else
    configStruct.xAxisMaxHistogram = '?';
end
handles.configStruct = configStruct;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_xMaxHistogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_xMaxHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_yMaxHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to edit_yMaxHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_yMaxHistogram as text
%        str2double(get(hObject,'String')) returns contents of edit_yMaxHistogram as a double

configStruct = handles.configStruct;
if ~isnan(str2double(get(hObject,'String')))
    configStruct.yAxisMaxHistogram = str2double(get(hObject,'String'));
else
    configStruct.yAxisMaxHistogram = '?';
end
handles.configStruct = configStruct;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_yMaxHistogram_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_yMaxHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_hideTickLabels.
function checkbox_hideTickLabels_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_hideTickLabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_hideTickLabels

configStruct = handles.configStruct;
configStruct.hideUnnecessaryTickLabels = hObject.Value;
handles.configStruct = configStruct;

guidata(hObject, handles);

% --- Executes on button press in radiobutton_displayPeakAmplitude.
function radiobutton_displayPeakAmplitude_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_displayPeakAmplitude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_displayPeakAmplitude
configStruct = handles.configStruct;

if hObject.Value == 1
    configStruct.displayPeakAmplitude = 1;
    configStruct.displayDecayTimeConstant = 0; 
    handles.radiobutton_displayDecay.Value = 0;
elseif hObject.Value == 0
%CG: don't allow all radio buttons to be off. Why would the user be
%configuring the settings for no data to be displayed?
    if handles.radiobutton_displayDecay.Value == 0
        hObject.Value = 1;
    elseif handles.radiobutton_displayDecay.Value == 1
        configStruct.displayPeakAmplitude = 0;
    end
end
handles.configStruct = configStruct;

guidata(hObject, handles);

% --- Executes on button press in radiobutton_displayDecay.
function radiobutton_displayDecay_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_displayDecay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_displayDecay

configStruct = handles.configStruct;

if hObject.Value == 1
    configStruct.displayDecayTimeConstant = 1;
    configStruct.displayPeakAmplitude = 0; 
    handles.radiobutton_displayPeakAmplitude.Value = 0;
elseif hObject.Value == 0
%CG: don't allow all radio buttons to be off. Why would the user be
%configuring the settings for no data to be displayed?
    if handles.radiobutton_displayPeakAmplitude.Value == 0
        hObject.Value = 1;
    elseif handles.radiobutton_displayPeakAmplitude.Value == 1
        configStruct.displayDecayTimeConstant = 0;
    end
end
handles.configStruct = configStruct;

guidata(hObject, handles);

function edit_binWidth_Callback(hObject, eventdata, handles)
% hObject    handle to edit_binWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_binWidth as text
%        str2double(get(hObject,'String')) returns contents of edit_binWidth as a double

if isnan(str2double(hObject.String))
    handles.edit_binWidth = 'default';
    handles.configStruct.binWidth = 'default';
elseif isnumeric(str2double(hObject.String))
    handles.edit_binWidth = str2double(hObject.String);
    handles.configStruct.binWidth = str2double(hObject.String);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_binWidth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_binWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_lengthOfDendritePerGraph_Callback(hObject, eventdata, handles)
% hObject    handle to edit_lengthOfDendritePerGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_lengthOfDendritePerGraph as text
%        str2double(get(hObject,'String')) returns contents of edit_lengthOfDendritePerGraph as a double

if isnan(str2double(hObject.String))
    handles.edit_lengthOfDendritePerGraph = 'default';
    handles.configStruct.lengthOfDendritePerGraph = 'default';
elseif isnumeric(str2double(hObject.String))
    handles.edit_lengthOfDendritePerGraph = str2double(hObject.String);
    handles.configStruct.lengthOfDendritePerGraph = str2double(hObject.String);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_lengthOfDendritePerGraph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_lengthOfDendritePerGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

handles.configStruct = handles.configStructOldVersion;
%CG: discard all changes to configStruct and then end configureGraphs.
handles.cancelButtonClicked = 1;
guidata(hObject, handles);
uiresume(handles.figure1);



function edit_maximumDistance_Callback(hObject, eventdata, handles)
% hObject    handle to edit_maximumDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_maximumDistance as text
%        str2double(get(hObject,'String')) returns contents of edit_maximumDistance as a double
configStruct = handles.configStruct;
if ~isnan(hObject.String)
    if ~isnan(hObject.String)
        hObject.ForegroundColor = [0 0 0]; 
        hObject.FontAngle = 'normal';
        configStruct.maximumDistance = str2double(hObject.String);
    elseif isnan(hObject.String) 
        hObject.String = 'invalid max distance value';
        hObject.ForegroundColor = [1 0 0];
        hObject.FontAngle = 'normal';
        configStruct.maximumDistance = [];
    end
elseif isnan(hObject.String)
    hObject.String = 'maximum distance';
    hObject.ForegroundColor = [0.502 0.502 0.502];
    hObject.FontAngle = 'italic';
    configStruct.maximumDistance = [];
end
handles.configStruct = configStruct;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_maximumDistance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_maximumDistance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

