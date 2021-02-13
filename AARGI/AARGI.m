function varargout = AARGI(varargin)

% Outline
% AARGI uses Mathworks' gui-building feature 'guide' to generate a GUI also
% named AARGI. AARGI gui contains all controls the user will need to
% segment his or her data and find the fluorescence intensity of the raw
% pixel data in each ROI across all frames. 

% Author: Charlie J Gilbride
% version: 1.0.0


% AARGI MATLAB code for AARGI.fig
%      AARGI, by itself, creates a new AARGI or raises the existing
%      singleton*.
%
%      H = AARGI returns the handle to a new AARGI or the handle to
%      the existing singleton*.
%
%      AARGI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AARGI.M with the given input arguments.
%
%      AARGI('Property','Value',...) creates a new AARGI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AARGI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AARGI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AARGI

% Last Modified by GUIDE v2.5 20-Apr-2019 11:12:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AARGI_OpeningFcn, ...
                   'gui_OutputFcn',  @AARGI_OutputFcn, ...
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


% --- Executes just before AARGI is made visible.
function AARGI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AARGI (see VARARGIN)

% Choose default command line output for AARGI

%CG: get operating system type for correct handling of directory strings.
osType = computer;
handles.osType = osType;

functionName = 'AARGI.m'; functionDir = which(functionName);
functionDir = functionDir(1:end-length(functionName));
addpath(genpath(functionDir)); 

if ~isempty(strfind(handles.osType, 'MAC')) 
    slashType = '/'; 
else 
    slashType = '\';
end
handles.slashType = slashType;

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

suffixStrs = {};
handles.suffixStrs = suffixStrs;
keywordStrs = {};
handles.keywordStrs = keywordStrs;
Dirs = {}; handles.Dirs = Dirs;

%CG: disable threshold intervals text boxes:
edit_XCorr = handles.edit_Incr;
edit_XCorr.Enable = 'off'; handles.edit_Incr = edit_XCorr;
edit_RawAmp = handles.edit_numDataChunks;
edit_RawAmp.Enable = 'off'; handles.edit_numDataChunks = edit_RawAmp;

% %************Development only
% handles.edit_MaxLat.String = '3'; handles.edit_Eval_thres.String = '5'; 
% handles.edit_RawAmpthres.String = '3.5'; handles.edit_PDE.String = '70';
% handles.edit_ROIsl.String = '3';
%************
handles.edit_MaxLat.String = 'Enter value'; handles.edit_Eval_thres.String = 'Enter value'; 
handles.edit_RawAmpthres.String = 'Enter value'; handles.edit_PDE.String = 'Enter value';
handles.edit_ROIsl.String = 'Enter value';

handles.edit_MaxLat.ForegroundColor = [0.502 0.502 0.502]; handles.edit_MaxLat.FontAngle = 'italic';
handles.edit_Eval_thres.ForegroundColor = [0.502 0.502 0.502]; handles.edit_Eval_thres.FontAngle = 'italic';  
handles.edit_RawAmpthres.ForegroundColor = [0.502 0.502 0.502]; handles.edit_RawAmpthres.FontAngle = 'italic'; 
handles.edit_PDE.ForegroundColor = [0.502 0.502 0.502]; handles.edit_PDE.FontAngle = 'italic'; 
handles.edit_ROIsl.ForegroundColor = [0.502 0.502 0.502]; handles.edit_ROIsl.FontAngle = 'italic'; 

%CG: default check for lateral shift
handles.checkbox_CheckLateralShift.Value = 1;


handles.checkbox_closeFigures.Value = 1;
%CG: By default, the figures created by the AARG function should be closed
%after ROIs have been established. 
% Update handles structure
handles.output.Visible = 'On';
set(findall(hObject, '-property', 'Units'), 'Units', 'Normalized')
set(hObject, 'Resize','on') 

guidata(hObject, handles);

% UIWAIT makes AARGI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AARGI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_selexps.
function pushbutton_selexps_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selexps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Dirs = uipickfiles('prompt','Pick experiments for AARGI');
% Dirs{1}
if iscell(Dirs)
    handles.Dirs = Dirs;
    edit_selexpDisp_Callback(handles.edit_selexpDisp, eventdata, handles);
end
%guidata(hObject, handles);

% --- Executes on button press in pushbutton_remexps.
function pushbutton_remexps_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remexps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_Help.
function pushbutton_Help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


open('HelpAARGI.pdf')
%--------1. Key word(s)------------

%CG: Keyword(s): Delete default text and enter
%your keywords*. Press 'enter' when ready. Don't put more than one space
%between each keyword. If there are some keywords being applicable to some
%data folders in experiment list but not all, AARGI will match data folders
%to the relevant key words before proceeding.
%*Keywords are any combination of keyboard strokes you have used to
%distinguish data files belonging to different conditions within the same
%experiment.

%CG: IMPORTANT: the order in which the keywords are entered into the
%keyword(s) text box determines the order in which they will be analysed.
%It is important that keywords are entered into the dialogue box in the
%same order as the data in the corresponding condition was acquired. Simple
%example: if condition 1 (keyword: 'baseline'), condition 2 (keyword:
%'drugA') and condition 3 (keyword: 'washout') compose Experiment 1, then
%the following keywords should be entered into the text box: 'baseline
%drugA washout' in that order. More complex example: if the user wishes to
%analyse data from Experiment 2 in the same batch as Experiment 1, but
%there is different keyword for condition 2 (namely, drugB), then an
%appropriate string of keywords would be: 'baseline drugA drugB washout'.
%When a keyword does not exist within a given experiment, that keyword will
%be skipped by AARGI for that experiment. 

%--------2. Experiment List------------

%CG: Experiment List:
%Select each data folder to be included in analysis. This will appear in
%the form of an abbreviated directory showing only the name of lowest
%folder in the directory (i.e. the selected data folder). Click the 'Select
%experiment' button to begin selecting data folders - these should be
%folders directly containing the raw data files (i.e. the data files should
%not be contained within subfolders within the selected folder). But see
%note on 'Bulk load' checkbox for an exception to this rule. Remove
%selected experiments from the list simply by highlighting them in the
%editable text box and deleting.

%CG: Bulk load: when checked, the user can
%select a parent folder containing multiple data folders rather than each
%individual data folder.

%--------3. Optimize thresholds------------

%CG: Multi thresholds (E-value): If checked, enables the user to set an
%interval value for the E-value threshold (see section 4 for more info on
%E-value threshold). In this case, the objective should be to find an
%appropriate or more appropriate E-value threshold that will maximise the
%number of events detected without introducing noise. When this checkbox is
%checked interval value should be specified in the 'Incr.' editable text
%box and instead of a single value entered in the value box for the E-value
%threshold (see section 4. Setup), the user should set a range of values in
%the form: 'i-j', where i is the smallest threshold value the user wishes
%to test and j is the largest. In addition, the user can limit the amount
%of data analysed. This is a very useful feature, because the user might
%want to test up to 10 different thresholds. Including all the available
%data in such a case might be very time consuming! By default, data is
%chunked into 500-frame segments. A sufficient number of chunks depends on
%the frame rate and frequency of events. The threshold value applied will
%increase from i in increments determined by the interval ('Incr.') value
%up to j. If this box is checked AARGI will detect events using the
%thresholds specified, while keeping the raw amplitude threshold fixed at a
%single value. ROIs will not be identified if the user is optimizing thresholds.

%CG: Multi thresholds (raw amp):
%same principle as for the multi threshold checkbox for E-values. A
%increment value should be given and the number of data chunks to be analysed
%can also be specified.

%CG: Specify keyword: threshold optimization will be applied to only one
%condition and the user should specify which condition that will be by
%entering one of the keywords from section 1. 

%--------4. Setup------------

%CG: Post-Detection Exclusion (PDE): The user should select the
%number of frames equal to the duration of the longest event. If an
%appropriate value is not known, the frequency distribution
%for event duration plots generated after the event detection algorithm is
%complete can be used to judge what an appropriate value should be. These
%graphs compile all events across all experiments, in each condition, into
%a single distribution.

%CG: Max. lateral shift permitted: lateral shift refers to steady small
%drift along the x/y axis during the experiment. If lateral shift along
%either axis is detected above the value entered here, then AARG will
%automatically correct it before detecting events and applying ROIs.

%CG: E-value threshold: Multiplication of
%cross-covariance and cross-correlation matrices produces the Evaluation or
%E matrix. The cross-covariance matrix has the same dimensions as the
%re-sized data chunk in which events are being detected. For each pixel in
%the data chunk, the covariance is calculated. This is the covariance
%between 37 pixel values ahead of this pixel and its corresponding value on
%each event template. The event templates are a range of event shapes that
%could plausibly occur within the data. If any given segment of data fits
%well with any one of these templates, then the covariance value for the
%relevant pixel will be high. This is the principle of the cross-covariance
%calculation. In practice, this calculation is computationally too
%demanding for conventional computers. To reduce computational demands,
%fast fourier transform is used. The cross-correlation is simply the
%cross-covariance matrix with each value divided by the standard deviation.
%of the data in each pixel across time. The cross-covariance and
%cross-correlation calculations have complementary strengths and weaknesses
%if used as two thresholds. To reduce the number of thresholds the user
%must set, the cross-covariance and cross-correlation matrices have been
%crossed with each other to produce the E matrix. 

%CG: Raw amplitude threshold: 2nd threshold is applied
%to the amplitude of the raw fluorescence intensity data to further reduce
%noise contamination. As with the E-value threshold, the user can set up
%multiple raw amplitude thresholds to be tested. This requires that E-value
%threshold is fixed at a single value. 

%CG: ROI side length: In AARG analysis, every ROI is a square
%with the same side length. This must be specified by the user. An
%appropriate value depends on the pixel resolution, magnification and
%objectives of the experiment. 

%CG: Close figures: Box is checked when gui opens. This will result in
%figures displaying placement of ROIs after AARG analysis to be closed
%(i.e. deleted). This is preferable when many experiments have been loaded
%in one sitting because it prevents too much RAM being used up by these
%figure windows. The user might want to uncheck this box and use the
%figures generated during analysis for illustration.

%--------5. Run AARGI------------
%CG: Check lateral shift: when active, the user will be prompted to select
%a region of the field of view for each experiment. This will be the
%reference image by which AARG measures the level of shift (comparing an
%average of the first three frames to the last frame). Checked by default.

%CG: Detect events: calls the functions which threshold the raw data to
%find events. Unless the user simply wishes to measure the level of lateral
%shift, this one should be checked. 

%CG: Check PDE: when ticked, events will be detected, but
%ROIs will not be established. This is intended to be used once an
%appropriate threshold combination is identified and the user wishes to
%find an appropriate inter-event interval value to add to the editable text
%box in section 4.

%CG: Find ROIs: this will prompt AARG to carry out its primary function,
%which is to use detected events to find ROIs. 

%CG: Launch: the functions, or some combination thereof (depending
%on selections made by the user), will be called. These functions do the
%following: test lateral shift, correct lateral shift, detect events with
%one or more threshold values, establish ROIs and gather events into
%common populations so that an appropriate PDE value can be selected. 

% --- Executes on button press in pushbutton_Launch.
function pushbutton_Launch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Launch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

keywordStrs = handles.keywordStrs; 
Dirs = handles.Dirs; text_status = handles.text_status;
checkbox_checkPDE = handles.checkbox_checkPDE;
ROI_sl = str2double(handles.edit_ROIsl.String);
PDE = str2double(handles.edit_PDE.String);

checkbox_MultiThresEval = handles.checkbox_MultiThresEval;
checkbox_MultiThresRawAmp = handles.checkbox_MultiThresRawAmp;

checkbox_closeFigures = handles.checkbox_closeFigures;

selectedExpDispStr = handles.edit_selexpDisp.String;
txtdim = size(selectedExpDispStr); pureStringSelexpDisp = '';
for cRowIdx = 1 : txtdim(1)
    if ~isempty(selectedExpDispStr)
        if ~strcmp(selectedExpDispStr(cRowIdx,1),'')
            pureStringSelexpDisp = strcat(pureStringSelexpDisp,selectedExpDispStr(cRowIdx,:));
        end
    end
end
%CG: make sure what is displayed in the select experiment display is what
%AARGI will analyse. 
NumDirs = size(Dirs,2); newDirs = {}; dirCounter = 0;
for cDirIdx = 1 : NumDirs
    [~,cCellName,~] = fileparts(Dirs{cDirIdx});
    if ~isempty(strfind(pureStringSelexpDisp,cCellName));
        dirCounter = dirCounter + 1; newDirs{dirCounter} = Dirs{cDirIdx};
    end
end
Dirs = newDirs; NumDirs = size(Dirs,2);
checkbox_bulkload = handles.checkbox_bulkload;

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
handles.edit_selexpDisp.String = sprintf(text);

AbortLaunch = 0;
%CG: make sure the user has specified all the necessary inputs for each
%section, starting with sections 1 and 2:
if isempty(keywordStrs) && ~isempty(Dirs)
    text_status.String = 'Enter keyword(s) to begin analysis';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
elseif isempty(Dirs) && ~isempty(keywordStrs)
    text_status.String = 'Select experiment(s) to begin analysis';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
elseif isempty(keywordStrs) && isempty(Dirs)
    text_status.String = 'Enter keyword(s) and select experiment(s) to begin analysis';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end

%CG: section 3:
edit_Eval_thres = handles.edit_Eval_thres; edit_RawAmpthres = handles.edit_RawAmpthres;

if (~isempty(strfind(edit_Eval_thres.String, '-')) || ~isempty(strfind(edit_RawAmpthres.String, '-'))) && ...
        (checkbox_MultiThresRawAmp.Value == 0 && checkbox_MultiThresEval.Value == 0)
%CG: in case the user has put a '-' in the 'E-value threshold' or 'Fluo.
%amplitude threshold' text box without having checked one of the checkboxes
%in section 3. This could mean they have found an appropriate threshold
%combination, but forgot to enter this properly in section 4. (Setup).
    text_status.String = 'Check one of the boxes in "3. Optimize thresholds" if you want to test different thresholds';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end
if checkbox_MultiThresEval.Value == 1 || checkbox_MultiThresRawAmp.Value == 1
    specifiedKeyword = handles.edit_enterKeyword.String;
    if strcmp(specifiedKeyword, '') || strcmp(specifiedKeyword, 'Enter keyword') 
        text_status.String = 'Enter ONE keyword in "Specify keyword" text box';
        text_status.ForegroundColor = [1 0 0];
        text_status.FontSize = 16; AbortLaunch = 1;
    end
    
    if checkbox_MultiThresEval.Value == 1 && isempty(strfind(edit_Eval_thres.String, '-'))
        text_status.String = 'Enter RANGE of E threshold values';
        text_status.ForegroundColor = [1 0 0];
        text_status.FontSize = 16; AbortLaunch = 1;
    end
    if checkbox_MultiThresRawAmp.Value == 1 && isempty(strfind(edit_RawAmpthres.String, '-'))
        text_status.String = 'Enter RANGE of raw amplitude threshold values';
        text_status.ForegroundColor = [1 0 0];
        text_status.FontSize = 16; AbortLaunch = 1;
    end
    
else
%CG: make specifiedKeyword empty to prompt AARGI to search for appropriate
%keyword in suffixStrs.
    specifiedKeyword = '';
end

%CG: section 4:
edit_MaxLat = handles.edit_MaxLat;
if isempty(edit_MaxLat.String) || strcmp(edit_MaxLat.String, 'Enter value')
    text_status.String = 'Enter maximum permitted lateral shift';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end

if isempty(edit_Eval_thres.String) || strcmp(edit_Eval_thres.String, 'Enter value')
    text_status.String = 'Enter E-value threshold';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end

if isempty(edit_RawAmpthres.String) || strcmp(edit_RawAmpthres.String, 'Enter value')
    text_status.String = 'Enter amplitude threshold';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end

edit_PDE = handles.edit_PDE;
if checkbox_checkPDE.Value == 0 && strcmp(edit_PDE.String,'Enter value')
    text_status.String = 'Enter PDE or use PDE checkbox';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end

edit_pBF = handles.edit_pixelBinningFactor;
if isempty(edit_pBF.String) || strcmp(edit_pBF.String, 'Enter value')
    text_status.String = 'Enter a pixel binning factor. This should be "1" if you do not want any binning';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end

edit_ROIsl = handles.edit_ROIsl;
if isempty(edit_ROIsl.String) || strcmp(edit_ROIsl.String,'Enter value') ||...
        ~mod(str2double(edit_ROIsl.String),2)
%CG: check a roi size has been specified unless the user wants to measure
%event sizes only.
    text_status.String = 'Enter a valid side length for ROIs (odd numbers only)';
    text_status.ForegroundColor = [1 0 0];
    text_status.FontSize = 16; AbortLaunch = 1;
end
if AbortLaunch == 0
%CG: first task is to match up the keywords and data directories.

%CG: go through each data folder and find the key words that apply to each
%data file contained therein. First, identify which combinations of key
%words apply.
    NumDirs = size(Dirs,2); NumKWStrs = size(keywordStrs,2); kwCombos = cell(NumDirs,NumKWStrs);
%CG: fill up kwCombos with detected combinations and then remove the
%duplicates.

    latshift = handles.checkbox_CheckLateralShift.Value;
    detectons = handles.checkbox_detection.Value;
    CheckPDE = handles.checkbox_checkPDE.Value;
    findROIs = handles.checkbox_findROIs.Value;

    for ckwIdx = 1 : NumKWStrs
        ckw = keywordStrs{ckwIdx};
        for cDirIdx = 1 : NumDirs
            cDir = Dirs{cDirIdx}; cd(cDir)
            dirInfo = dir; NumItems = size(dirInfo,1);
            for cItemIdx = 1 : NumItems
                cItem = dirInfo(cItemIdx).name;
                if ~isempty(strfind(cItem,ckw))
                    kwCombos{cDirIdx,ckwIdx} = ckw;
                end
            end
        end
    end
    kwCombos_copy = kwCombos; numberOfRows = size(kwCombos_copy,1); kwCombos = {};
%CG: in case the user enters: a b a c, instead of: a b c, where a is the string denoting the baseline
%condition and b and c are some other conditions, the for-loop for cRowIdx = 1 : numberOfRows
%renders kwCombos what it would be had the user written a b c
    for cRowIdx = 1 : numberOfRows
        row_kwCombos_copy = kwCombos_copy(cRowIdx,:);
        numberOfCells = size(row_kwCombos_copy, 2);
        for cCellIdx = 1 : numberOfCells
            if isempty(kwCombos) && ~isempty(row_kwCombos_copy{cCellIdx})
                kwCombos{1} = row_kwCombos_copy{cCellIdx};
            elseif ~isempty(kwCombos) && ~isempty(row_kwCombos_copy{cCellIdx})
                
                if cCellIdx == 1
                    kwCombos{cRowIdx,1} = row_kwCombos_copy{cCellIdx};
                elseif cCellIdx > 1 
                    cCellIdx2 = 1; addIt = 1;
                    while cCellIdx2 <= size(kwCombos,2)
                        if strcmp(row_kwCombos_copy{cCellIdx}, kwCombos{cRowIdx,cCellIdx2})
                            addIt = 0;
                        end
                        cCellIdx2 = cCellIdx2 + 1;
                    end
                    if addIt == 1; kwCombos{cRowIdx,end+1} = row_kwCombos_copy{cCellIdx}; end
                end
                
            end
        end
    end
    
%CG: each row in kwCombos corresponds to a directory in Dirs, which
%provides access to the data file. If there are two non-empty cells in row
%1, this means the first directory leads to an experiment with two
%conditions (assuming each condition is designated by a keyword). 
    NumKWCols = size(kwCombos,2); CellNameList = {}; allsuffixStrs = {};
    medianVals = []; CorrectShiftSwt_list = {};
    for cDirIdx = 1 : NumDirs 
        cDir = Dirs{cDirIdx}; strcount = 0; suffixStrs = {}; specifiedConditionIdx = [];
        genpath(cDir)
%         [~,CellName,~] = fileparts(cDir);
        for cColIdx = 1 : NumKWCols
            if ~isempty(kwCombos{cDirIdx,cColIdx})
                strcount = strcount + 1; suffixStrs{strcount} = kwCombos{cDirIdx,cColIdx};
                if ~isempty(specifiedKeyword)
                    if strcmp(suffixStrs{strcount},specifiedKeyword)
                        specifiedConditionIdx = strcount;
                    end
                end
            end
        end
        if ~isempty(specifiedKeyword) && isempty(specifiedConditionIdx)
            text_status.String = strcat('Specified keyword not found. Ensure specified keyword has been listed',...
                '\n in the Keyword(s) section (section 1).');
            text_status.ForegroundColor = [1 0 0];
            text_status.FontSize = 16; AbortLaunch = 1;
        end
            
        allsuffixStrs{cDirIdx,1} = suffixStrs;
        if AbortLaunch == 0 
            
            if latshift == 1 && detectons == 0 && CheckPDE == 0 && findROIs == 0
                text_status.String = 'Checking lateral shift...(see command line for output)';
                text_status.ForegroundColor = [0.502 0.502 0.502]; text_status.FontSize = 12;
                [CellNameList, medianVals] = LatShiftCheck(cDir,suffixStrs,...
                    CellNameList,medianVals,checkbox_closeFigures.Value);

            elseif latshift == 1 && detectons == 1
%                 edit_MaxLat = edit_MaxLat.handles; edit_Eval_thres = edit_Eval_thres.handles;
%                 edit_RawAmpThres = edit_RawAmpThres.handles;
                if strcmp(edit_MaxLat.String,'Enter value') || strcmp(edit_MaxLat.String,'')
                    text_status.String = 'Enter maximum acceptable lateral shift';
                    text_status.ForegroundColor = [1 0 0];
                    text_status.FontSize = 16; AbortLaunch = 1;
                elseif strcmp(edit_Eval_thres.String,'Enter value') || strcmp(edit_Eval_thres.String,'')
                    text_status.String = 'Enter E-value threshold';
                    text_status.ForegroundColor = [1 0 0];
                    text_status.FontSize = 16; AbortLaunch = 1;
                elseif strcmp(edit_RawAmpthres.String,'Enter value') || strcmp(edit_RawAmpthres.String,'')
                    text_status.String = 'Enter raw amplitude threshold';
                    text_status.ForegroundColor = [1 0 0];
                    text_status.FontSize = 16; AbortLaunch = 1;
                else
                    if ~isempty(strfind(edit_Eval_thres.String, '-')) && ~isempty(strfind(edit_RawAmpthres.String, '-'))
%CG: one threshold must be kept constant. 
                        text_status.String = 'Re-assign one threshold to a single value';
                        text_status.ForegroundColor = [1 0 0];
                        text_status.FontSize = 16; AbortLaunch = 1;
                    else
                        text_status.String = 'Checking lateral shift...';
                        [CellNameList, medianVals] = LatShiftCheck(cDir,suffixStrs,...
                            CellNameList,medianVals,checkbox_closeFigures.Value);
                        MaxLatShift = str2double(edit_MaxLat.String);
                        if ~isempty(find(medianVals==MaxLatShift,1))
                            text_status.String = strcat('Lateral shift > "',edit_MaxLat.String,...
                                '" detect in "',CellNameList{end},'". This shift will be corrected.');
                            text_status.ForegroundColor = [1 0 0];
                            text_status.FontSize = 16; AbortLaunch = 1;
                            CorrectShiftSwt_list{cDirIdx,1} = 'On';
                        elseif isempty(find(medianVals==MaxLatShift,1))
                            CorrectShiftSwt_list{cDirIdx,1} = 'Off';
                        end
                    end
                end
            elseif latshift == 0 
                CorrectShiftSwt_list{cDirIdx,1} = 'Off';
            end
        
        end
    end
    if AbortLaunch == 0 
        text_status.String = 'Launching...';
        text_status.ForegroundColor = [0.502 0.502 0.502]; text_status.FontSize = 12;
    end
    
    if AbortLaunch == 0 && detectons == 1
        ThresholdDataFolderList = {}; strlists = struct;
        for cDirIdx = 1 : NumDirs
            cDir = Dirs{cDirIdx}; suffixStrs = allsuffixStrs{cDirIdx,1};
            
%             datadirs = {}; datadirs{1} = cDir; 
            cd(cDir); dirInfo = dir; 
            ItemList = {}; NumItems = size(dirInfo,1); itemcount = 0;
            for cItemIdx = 1 : NumItems
                cItem = dirInfo(cItemIdx).name; 
%CG: first get all tif files together.
                if ~isempty(strfind(cItem,'.tif'))
                    itemcount = itemcount + 1; 
                    ItemList{itemcount,1} = cItem;
                end
            end
%CG: isolate tif files with some unidentified characters between keyword
%and file format string. 
            NumTifs = size(ItemList,1); specialstrs = {}; specialcount = 0;
            for cTifIdx = 1 : NumTifs
                strcount = 0; strfound = 0;
                while  strfound == 0 && strcount < size(suffixStrs,2)
                    strcount = strcount + 1;
                    if isempty(strfind(ItemList{cTifIdx},strcat(suffixStrs{strcount},'.tif'))) &&...
                            ~isempty(strfind(ItemList{cTifIdx},suffixStrs{strcount}))

                        strfound = 1; 
                    end
                end
                if strfound == 1; 
                    specialcount = specialcount + 1; 
                    specialstrs{specialcount,1} = ItemList{cTifIdx}; 
                end
            end
            NumSpecStrs = size(specialstrs,1); 
%CG: isolate the unidentified string and search from back of string to find
%first 0. Final strings should be identical once the first 0 is found. 
            IDcount = 0; unID = {}; 
            for cSpecIdx = 1 : NumSpecStrs
                strcount = 0; strfound = 0; cSpecStr = specialstrs{cSpecIdx};
                while strcount < size(suffixStrs,2) && strfound == 0
                    strcount = strcount + 1;
                    stridx = strfind(cSpecStr,suffixStrs{strcount});
                    if ~isempty(stridx)
                        tifidx = strfind(cSpecStr,'.tif'); IDcount = IDcount + 1;
                        strID = cSpecStr(stridx+numel(suffixStrs{strcount}):tifidx);
                        zeroidx = strfind(strID,'0'); 
                        if ~isempty(zeroidx)
                            strID = strID(1:zeroidx(end));
                            if IDcount == 1;
                                unID{IDcount} = strID;
                            else
                                if ~strcmp(unID{1},strID)
                                    unID{IDcount} = strID;
                                end
                            end
                        end
                    end
                end
            end
            if ~isempty(unID)
                if size(unID,2) > 1
                    errordlg('partnamesuff not found')
                elseif size(unID,2) == 1
                    partnamesuff = unID{1};
                end
            else
                partnamesuff = '';
            end
%CG: In the case of large tiff files that have been split by the
%acquisition program, the user can define "partnamesuff" such that all
%segments of the split file are converted to .mat files. 
            [~,CellName,~] = fileparts(cDir);
            if ~isempty(strfind(edit_Eval_thres.String, '-')) && isempty(strfind(edit_RawAmpthres.String, '-'))
%CG: this means the user wishes to evalutate multiple E-value thresholds
%with a single raw amplitude threshold.
                
                edit_Incr = handles.edit_Incr; ChunkLimit = str2double(handles.edit_numDataChunks.String);
                CCthresStr = edit_Eval_thres.String; rawampStr = edit_RawAmpthres.String;  
                if strcmp(ChunkLimit,'all') || strcmp(ChunkLimit,'All') || strcmp(ChunkLimit,'ALL')
                    ThresOptsID = strcat('E_val_', CCthresStr,'_Ampl',rawampStr,'_allChunks');
                    ChunkLimit = NAN; %CG: might be expecting NAN later for this case. Would make more sense to re-assign ChunkLimit='all' otherwise.
                elseif isnan(ChunkLimit)
                    ThresOptsID = strcat('E_val_', CCthresStr,'_Ampl', rawampStr,'_1Chunk');
                else
                    ThresOptsID = strcat('E_val_', CCthresStr,'_Ampl', rawampStr,'_',...
                        num2str(ChunkLimit),'Chunks');
                end
                
                ThresSet = 'E_val';
                
                hyphenIdx = strfind(CCthresStr, '-');
                start_Ethres = str2double(CCthresStr(1:hyphenIdx-1));
                end_Ethres = str2double(CCthresStr(hyphenIdx+1:end));
                Incr = edit_Incr.String; if strcmp(Incr(end),'.'); Incr(end) = []; end
                if start_Ethres > end_Ethres
                    E_thres = start_Ethres:-str2double(Incr):end_Ethres;
                else
                    E_thres = start_Ethres:str2double(Incr):end_Ethres;
                end
                Ampl_thres = str2double(edit_RawAmpthres.String);
                pause(5)
%CG: Before detecting events check that this combination of thresholds has not already been plotted
                if ~isempty(strfind(handles.osType, 'MAC')) 
                    Slashes = strfind(cDir, '/'); 
                else 
                    Slashes = strfind(cDir, '\');
                end 
                ParentDir = cDir(1:Slashes(end)-1); cd(ParentDir)
                dirInfo = dir; numItems = size(dirInfo,1); cItemIdx = 0;
                tRecs = strcat('ThresholdRecords_',CellName); plotFound = 0;
                while cItemIdx < numItems && plotFound == 0
                    cItemIdx = cItemIdx + 1; cItem = dirInfo(cItemIdx).name;
                    if strcmp(cItem, tRecs)
                        strcat(ParentDir,'/',tRecs); cd(strcat(ParentDir,'/',tRecs));
                        subDirInfo = dir; subNumItems = size(subDirInfo,1); cSubItemIdx = 0;
                        while cSubItemIdx < subNumItems && plotFound == 0
                            cSubItemIdx = cSubItemIdx + 1;
                            cSubItem = subDirInfo(cSubItemIdx).name;
                            if strcmp(cSubItem,ThresOptsID)
                                plotFound = 1;
                                text_status.String = 'Plotting detections...';
                                text_status.ForegroundColor = [0.502 0.502 0.502]; 
                                text_status.FontSize = 12;
                            end
                        end
                    end
                end
                cd(cDir)
                PTOin = struct;
                PTOin.suffixStrs = suffixStrs;
                PTOin.specifiedConditionIdx = specifiedConditionIdx;
                PTOin.ThresSet = ThresSet;
                PTOin.ThresOptsID = ThresOptsID;
                PTOin.cDir = cDir;
                
                if plotFound == 0
                    text_status.String = 'Detecting...';
                    text_status.ForegroundColor = [0.502 0.502 0.502];
                    text_status.FontSize = 12;
                    [handles,ThresOptsID] = Detection(cDir,cDirIdx,suffixStrs,...
                        CorrectShiftSwt_list{cDirIdx},E_thres,Ampl_thres,...
                        handles, partnamesuff,ChunkLimit,ThresholdDataFolderList,...
                        strlists,specifiedConditionIdx,ThresOptsID,CellName);
                    
                    PTOin.ThresOptsID = ThresOptsID;
                    PTOin.options = 'write';
                    
                    PlotThresholdOptions(PTOin)
                else
                    PTOin.options = 'read';
                    PTOin.CellName = CellName;
                    PTOin.slashType = handles.slashType;
                    PlotThresholdOptions(PTOin)
                end
            elseif isempty(strfind(edit_Eval_thres.String, '-')) && ~isempty(strfind(edit_RawAmpthres.String, '-'))
%CG: this means the user wishes to evaluate multiple raw amplitude
%thresholds with a single E-value threshold.
                edit_Incr = handles.edit_Incr; ChunkLimit = str2double(handles.edit_numDataChunks.String);
                rawampStr = edit_RawAmpthres.String; CCthresStr = edit_Eval_thres.String; 
                if strcmp(ChunkLimit,'all') || strcmp(ChunkLimit,'All') || strcmp(ChunkLimit,'ALL')
                    ThresOptsID = strcat('E_val_', CCthresStr,'_Ampl',rawampStr,'_allChunks');
                    ChunkLimit = NAN; %CG: might be expecting NAN later for this case. Would make more sense to re-assign ChunkLimit='all' otherwise.
                elseif isnan(ChunkLimit)
                    ThresOptsID = strcat('E_val_', CCthresStr,'_Ampl', rawampStr,'_1Chunk');
                else
                    ThresOptsID = strcat('E_val_', CCthresStr,'_Ampl', rawampStr,'_',...
                        num2str(ChunkLimit),'Chunks');
                end

                ThresSet = 'Ampl';

                hyphenIdx = strfind(rawampStr,'-');
                start_Ampl_thres = str2double(rawampStr(1:hyphenIdx-1));
                end_Ampl_thres = str2double(rawampStr(hyphenIdx+1:end));
                E_thres = str2double(CCthresStr);
                Incr = edit_Incr.String; if strcmp(Incr(end),'.'); Incr(end) = []; end
                if start_Ampl_thres > end_Ampl_thres
                    Ampl_thres = start_Ampl_thres:-str2double(Incr):end_Ampl_thres;
                else
                    Ampl_thres = start_Ampl_thres:str2double(Incr):end_Ampl_thres;
                end
                pause(5)
%CG: Before detecting events check that this combination of thresholds has not already been plotted
                if ~isempty(strfind(handles.osType, 'MAC')) 
                    Slashes = strfind(cDir, '/'); 
                else 
                    Slashes = strfind(cDir, '\');
                end
                ParentDir = cDir(1:Slashes(end)-1); cd(ParentDir)
                dirInfo = dir; numItems = size(dirInfo,1); cItemIdx = 0;
                tRecs = strcat('ThresholdRecords_',CellName); plotFound = 0;
                while cItemIdx < numItems && plotFound == 0
                    cItemIdx = cItemIdx + 1; cItem = dirInfo(cItemIdx).name;
                    if strcmp(cItem, tRecs)
                        strcat(ParentDir,'/',tRecs); cd(strcat(ParentDir,'/',tRecs));
                        subDirInfo = dir; subNumItems = size(subDirInfo,1); cSubItemIdx = 0;
                        while cSubItemIdx < subNumItems && plotFound == 0
                            cSubItemIdx = cSubItemIdx + 1;
                            cSubItem = subDirInfo(cSubItemIdx).name;
                            if strcmp(cSubItem,ThresOptsID)
                                plotFound = 1;
                                text_status.String = 'Plotting detections...';
                                text_status.ForegroundColor = [0.502 0.502 0.502]; 
                                text_status.FontSize = 12;
                            end
                        end
                    end
                end
                cd(cDir)
                PTOin = struct;
                PTOin.suffixStrs = suffixStrs;
                PTOin.specifiedConditionIdx = specifiedConditionIdx;
                PTOin.ThresSet = ThresSet;
                PTOin.ThresOptsID = ThresOptsID;
                PTOin.cDir = cDir;
                
                if plotFound == 0
                    text_status.String = 'Detecting...';
                    text_status.ForegroundColor = [0.502 0.502 0.502];
                    text_status.FontSize = 12;
                    [handles,ThresOptsID] = Detection(cDir,cDirIdx,suffixStrs,...
                        CorrectShiftSwt_list{cDirIdx},E_thres,Ampl_thres,...
                        handles, partnamesuff,ChunkLimit,ThresholdDataFolderList,...
                        strlists,specifiedConditionIdx,ThresOptsID,CellName);
                    
                    PTOin.ThresOptsID = ThresOptsID;
                    PTOin.options = 'write';
                    PlotThresholdOptions(PTOin)
                else
                    PTOin.options = 'read';
                    PTOin.CellName = CellName;
                    PTOin.slashType = handles.slashType;
                    PlotThresholdOptions(PTOin)
                end
            elseif isempty(strfind(edit_Eval_thres.String, '-')) && isempty(strfind(edit_RawAmpthres.String, '-'))
%CG: just a single pair of thresholds is being used to detect events.
                pause(5); ChunkLimit = NaN; ThresOptsID = '';
                E_thres = str2double(edit_Eval_thres.String); 
                Ampl_thres = str2double(edit_RawAmpthres.String);
                text_status.String = 'Detecting...';
                text_status.ForegroundColor = [0.502 0.502 0.502];
                text_status.FontSize = 12;
                %%%testing
%                 cDir = '\\C:Volumes\No2\MATLAB\AARG_Data\ReAnalysis\conotoxin\CG0706173';
                [handles, ~, ThresholdDataFolderList,strlists] = Detection(cDir,cDirIdx,suffixStrs,...
                    CorrectShiftSwt_list{cDirIdx},E_thres,Ampl_thres,...
                    handles, partnamesuff,ChunkLimit,ThresholdDataFolderList,...
                    strlists,specifiedConditionIdx,ThresOptsID,CellName);
            end
        end
        
        if CheckPDE == 1
            [CA_dataoutput] = EventDurations(ThresholdDataFolderList,strlists);
            ShowEventDurations(CA_dataoutput,strlists,suffixStrs);
        end
        if isempty(strfind(edit_Eval_thres.String, '-')) &&...
                isempty(strfind(edit_RawAmpthres.String, '-')) &&...
                    findROIs == 1
            
%CG: trigger AARG and measurement of raw fluorescence intensity as long as
%the user is not testing multiple thresholds or wanting to see what PDE
%would be appropriate.
            for cDirIdx = 1 : NumDirs
                cDir = Dirs{cDirIdx}; suffixStrs = allsuffixStrs{cDirIdx,1};
                pause(5)
                text_status.String = 'Assigning ROIs...';
                text_status.ForegroundColor = [0.502 0.502 0.502];
                text_status.FontSize = 12; datadirs{1} = cDir;
                AARG_CoreHandler('datadirs', datadirs, 'suffixStrs', suffixStrs,...
                    'partnamesuff', partnamesuff, 'ROI_sl', ROI_sl,...
                    'PDE', PDE, 'figureFate', checkbox_closeFigures.Value, 'osType', handles.osType)
                pause(10)
                
%CG: probably helps to have some long pauses when large batches of data are
%being processed
                text_status.String = 'Measuring fluorescence intensity in each ROI...';
                text_status.ForegroundColor = [0.502 0.502 0.502];
                text_status.FontSize = 12; datadirs{1} = cDir;
                pause(2)
                handles.cDir = cDir;
                GetIntensities(handles,'suffixStrs',suffixStrs,'partnamesuff', partnamesuff,...
                    'DataType', 'Raw', 'osType', handles.osType);
                pause(10)
            end
        end
    end                
    
    if latshift == 1 && detectons == 0 && CheckPDE == 0 && findROIs == 0
        xaxisshift = medianVals(:,1); yaxisshift = medianVals(:,2);
        ShiftTable = table(xaxisshift, yaxisshift, 'RowNames', CellNameList);
        ShiftTable.Properties.VariableNames = {'pixel_shift_Xaxis', 'pixel_shift_Yaxis'};
    end
    text_status.String = 'Idle'; text_status.ForegroundColor = [0.502 0.502 0.502]; text_status.FontSize = 12;
end
guidata(hObject, handles);

% --- Executes on button press in pushbutton_Abort.
function pushbutton_Abort_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Abort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox_bulkload.
function checkbox_bulkload_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_bulkload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_bulkload

% if hObject.Value == 1
%     Dirs = handles.Dirs; NumDirs = size(Dirs,2);
%     for cDirIdx = 1 : NumDirs
%         
%     end
% end

function edit_MaxLat_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxLat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxLat as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxLat as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; handles.edit_MaxLat.FontAngle = 'normal';
else
    hObject.String = 'Enter value';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_MaxLat.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_MaxLat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxLat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Eval_thres_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Eval_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Eval_thres as text
%        str2double(get(hObject,'String')) returns contents of edit_Eval_thres as a double

if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; handles.edit_Eval_thres.FontAngle = 'normal';
else
    hObject.String = 'Enter value(s)';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_Eval_thres.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_Eval_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Eval_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_MultiThresEval.
function checkbox_MultiThresEval_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MultiThresEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MultiThresEval

cb_switch = hObject.Value; %cb_otherSwitch = handles.checkbox_MultiThresRawAmp.Value;
edit_Incr = handles.edit_Incr; edit_numDataChunks = handles.edit_numDataChunks;
if cb_switch == 1
    edit_Incr.Enable = 'on'; edit_numDataChunks.Enable = 'on'; cb_otherSwitch = 0;
    handles.checkbox_MultiThresRawAmp.Value = cb_otherSwitch;
    
    handles.checkbox_CheckLateralShift.Value = 0;
    handles.checkbox_CheckLateralShift.Enable = 'off';
    handles.checkbox_detection.Value = 1;
    handles.checkbox_detection.Enable = 'on';
    handles.checkbox_checkPDE.Value = 0;
    handles.checkbox_checkPDE.Enable = 'off';
    handles.checkbox_findROIs.Value = 0;
    handles.checkbox_findROIs.Enable = 'off';
    
elseif cb_switch == 0 %&& cb_otherSwitch == 0
    edit_Incr.Enable = 'off'; edit_numDataChunks.Enable = 'off';
    
    handles.checkbox_CheckLateralShift.Enable = 'on';
    handles.checkbox_detection.Enable = 'on';
    handles.checkbox_checkPDE.Enable = 'on';
    handles.checkbox_findROIs.Enable = 'on';
    
end
handles.edit_Incr = edit_Incr; handles.edit_numDataChunks = edit_numDataChunks;
guidata(hObject, handles); 


% --- Executes on button press in checkbox_CheckLateralShift.
function checkbox_CheckLateralShift_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CheckLateralShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_CheckLateralShift


function edit_keyword_Callback(hObject, eventdata, handles)
% hObject    handle to edit_keyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_keyword as text
%        str2double(get(hObject,'String')) returns contents of edit_keyword as a double
keywordStrs = handles.keywordStrs;

kw_str = hObject.String; space_idces = strfind(kw_str, ' '); 

%CG: contents of the editable text box must be something else other than
%the default contents.
if ~strcmp(kw_str, 'Enter key word(s) here...')

%CG: if the user puts a space after the last keyword, the number of
%keywords will be: NumStrs = numel(space_idces);
    if isempty(kw_str) 
        NumStrs = numel(space_idces);
    elseif ~strcmp(kw_str(end),' ')
        NumStrs = numel(space_idces) + 1;
    elseif strcmp(kw_str(end),' ') 
        NumStrs = numel(space_idces);
    end
    
%CG: any spaces before the first keyword should be cleared
    if ~isempty(kw_str)
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
            keywordStrs{1} = kw_str;
    end
    guidata(hObject, handles);
    if size(keywordStrs,2) == 1
        disp(strcat('"',keywordStrs{1},'"'))
    elseif size(keywordStrs,2) == 2
        disp(strcat('"',keywordStrs{1}, '" and "', keywordStrs{2}, '"'))
    end
    
    hObject.ForegroundColor = [0 0 0];
    hObject.FontAngle = 'normal';
    
end

if isempty(hObject.String)
    hObject.String = 'Enter keyword(s) here';
    hObject.ForegroundColor = [0.502 0.502 0.502]; 
    hObject.FontAngle = 'italic';
end

handles.keywordStrs = keywordStrs;
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_keyword_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_keyword (see GCBO)
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

Dirs = handles.Dirs; Sz_Dirs = size(Dirs); NumDirs = Sz_Dirs(1)*Sz_Dirs(2);
checkbox_bulkload = handles.checkbox_bulkload;
% Dirs = {}; Dirs{1} = '\\C:Volumes\No2\MATLAB\AARG_Data\ReAnalysis\conotoxin\CG0706173';
if handles.pushbutton_selexps.Value == 1
    if checkbox_bulkload.Value == 0
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
    elseif checkbox_bulkload.Value == 1 
        dfcount = 0; text = ''; rDirs = {}; rDir_counter = 0;
        for cDirIdx = 1 : NumDirs
            cDir = Dirs{cDirIdx}; cd(cDir); 
            dirInfo = dir; NumItems = size(dirInfo,1); dfcount = dfcount + NumItems;

            for cItemIdx = 1 : NumItems
                cItem = dirInfo(cItemIdx).name;
                if exist(cItem) == 7 && ~strcmp(cItem,'.') && ~strcmp(cItem,'..') && isempty(strfind(cItem,'ThresholdRecords'))
                    
                    if ~isempty(strfind(handles.osType, 'MAC')) 
                        sh_cDir = strcat('/.../',cItem);
                    else
                        sh_cDir = strcat('...',cItem);
                    end
                    
                    text = strcat(text,sh_cDir,'\n');
                    rDir_counter = rDir_counter + 1; 
                    rDirs{rDir_counter} = strcat(cDir,'/',cItem);
                end
            end
        end
        hObject.Max = dfcount; handles.Dirs = rDirs;
    end
    
    hObject.ForegroundColor = [0 0 0];
    hObject.FontAngle = 'normal';
    
    hObject.String = sprintf(text);
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
function edit_selexpDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_selexpDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Incr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Incr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Incr as text
%        str2double(get(hObject,'String')) returns contents of edit_Incr as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; handles.edit_Incr.FontAngle = 'normal';
else
    hObject.String = 'Incr.';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_Incr.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_Incr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Incr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white'); 
end



function edit_numDataChunks_Callback(hObject, eventdata, handles)
% hObject    handle to edit_numDataChunks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_numDataChunks as text
%        str2double(get(hObject,'String')) returns contents of edit_numDataChunks as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; handles.edit_numDataChunks.FontAngle = 'normal';
else
    hObject.String = '# chunks';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_numDataChunks.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_numDataChunks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_numDataChunks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_SkipCheck.
function checkbox_SkipCheck_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_SkipCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_SkipCheck

function edit_RawAmpthres_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RawAmpthres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RawAmpthres as text
%        str2double(get(hObject,'String')) returns contents of edit_RawAmpthres as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; handles.edit_RawAmpthres.FontAngle = 'normal';
else
    hObject.String = 'Enter value(s)';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_RawAmpthres.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_RawAmpthres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RawAmpthres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_PDE_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PDE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PDE as text
%        str2double(get(hObject,'String')) returns contents of edit_PDE as a double

if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; handles.edit_PDE.FontAngle = 'normal';
else
    hObject.String = 'Enter value';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1);
    handles.edit_PDE.FontAngle = 'italic';
end
    
% --- Executes during object creation, after setting all properties.
function edit_PDE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PDE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ROIsl_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ROIsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ROIsl as text
%        str2double(get(hObject,'String')) returns contents of edit_ROIsl as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; pause(0.1)
    handles.edit_ROIsl.FontAngle = 'normal';
else
    hObject.String = 'Enter value';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_ROIsl.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_ROIsl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ROIsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_MultiThresRawAmp.
function checkbox_MultiThresRawAmp_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_MultiThresRawAmp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_MultiThresRawAmp

cb_switch = hObject.Value;
edit_Incr = handles.edit_Incr; edit_numDataChunks = handles.edit_numDataChunks;
if cb_switch == 1
    edit_Incr.Enable = 'on'; edit_numDataChunks.Enable = 'on'; 
    cb_otherSwitch = 0; handles.checkbox_MultiThresEval.Value = cb_otherSwitch;
    
    handles.checkbox_CheckLateralShift.Value = 0;
    handles.checkbox_CheckLateralShift.Enable = 'off';
    handles.checkbox_detection.Value = 1;
    handles.checkbox_detection.Enable = 'on';
    handles.checkbox_checkPDE.Value = 0;
    handles.checkbox_checkPDE.Enable = 'off';
    handles.checkbox_findROIs.Value = 0;
    handles.checkbox_findROIs.Enable = 'off';
    
elseif cb_switch == 0
    edit_Incr.Enable = 'off'; edit_numDataChunks.Enable = 'off'; 
    
    handles.checkbox_CheckLateralShift.Enable = 'on';
    handles.checkbox_detection.Enable = 'on';
    handles.checkbox_checkPDE.Enable = 'on';
    handles.checkbox_findROIs.Enable = 'on';
end
handles.edit_Incr = edit_Incr; handles.edit_numDataChunks = edit_numDataChunks;
guidata(hObject, handles); 


% --- Executes on button press in checkbox_checkPDE.
function checkbox_checkPDE_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_checkPDE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_checkPDE


% --- Executes on button press in checkbox_detection.
function checkbox_detection_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_detection


% --- Executes on button press in checkbox_findROIs.
function checkbox_findROIs_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_findROIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_findROIs

if hObject.Value == 1
    handles.checkbox_detection.Value = 1;
    handles.checkbox_detection.Enable = 'off';
elseif hObject.Value == 0
    handles.checkbox_detection.Enable = 'on';
end
guidata(hObject, handles); 



function edit_enterKeyword_Callback(hObject, eventdata, handles)
% hObject    handle to edit_enterKeyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_enterKeyword as text
%        str2double(get(hObject,'String')) returns contents of edit_enterKeyword as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; pause(0.1)
    handles.edit_enterKeyword.FontAngle = 'normal';
else
    hObject.String = 'Enter keyword';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_enterKeyword.FontAngle = 'italic';
end



% --- Executes during object creation, after setting all properties.
function edit_enterKeyword_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_enterKeyword (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_retrieve.
function pushbutton_retrieve_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_retrieve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function edit_pixelBinningFactor_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pixelBinningFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pixelBinningFactor as text
%        str2double(get(hObject,'String')) returns contents of edit_pixelBinningFactor as a double
if isnumeric(str2double(hObject.String)) && ~isempty(hObject.String)
    hObject.ForegroundColor = [0 0 0]; pause(0.1)
    handles.edit_pixelBinningFactor.FontAngle = 'normal';
else
    hObject.String = 'Enter value';
    hObject.ForegroundColor = [0.502 0.502 0.502]; pause(0.1)
    handles.edit_pixelBinningFactor.FontAngle = 'italic';
end

% --- Executes during object creation, after setting all properties.
function edit_pixelBinningFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pixelBinningFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox_closeFigures.
function checkbox_closeFigures_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_closeFigures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_closeFigures

