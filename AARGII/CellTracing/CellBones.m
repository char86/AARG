function CellBones(varargin)


%Outline:
%CellBones connects each AARG-generated ROI (within 5µm from the dendrite)
%with any point or points on the image designated by the user. Check the
%command line for further instructions during function call. 
%1) Initially, the user is prompted to set 'start points', which appear as
%purple 'x' markers. 
%2) with the start points in place, the user can manually insert white
%lines or "bones". In the case of neurons, these should trace out the 
%dendrites. 
%3) after a double right click from the user, CellBones automatically
%generates green connecting lines between each ROI and the point on a white
%line. To carry out this function, CellBones has the following features:
%   i) MaxDistance is a user input that specifies the minimum distance the
%   nearest point on any line must be in order for this line to be
%   processed any further as a possible connection site between the current
%   green line and a white line. 
%   ii) to build a green line, CellBones finds all neighboring pixels that
%   are closer to any point on a line than itself. Of these neighbors, the
%   one with the highest median fluorescence intensity value is chosen as
%   the next point along the green line.
%4) The user has the chance to delete ill-fitting green lines. The cell
%array containing data on all the automatically generated green lines
%(those kept and those deleted) is saved as 'unmodified_CA_ROIc'.
%5) The actual distances are measured between each ROI and its connecting
%start point. Lines connecting the selected ROI to the start point briefly flash
%blue. Note that the entire white line will flash blue, not just the
%segment between the point of intersection and the start point (or next
%point of intersection). In ImageJ: open the raw data file, set the scale
%and build a segemented line to where the ROI would be. Use the measurement
%feature in ImageJ to check that for agreement. 

% Author: Charlie J Gilbride
% version: 1.0.0

p = inputParser;
% addParameter(p,'IDString','CG')
addParameter(p,'datadir',{})
addParameter(p,'suffixStrs',{'b_TTAP2'})
addParameter(p,'MaxDistance',5) %CG: MaxDistance in um.
addParameter(p,'Resolution', 6.25) %CG: 6.25 pixels/um at 100x magnification.
addParameter(p,'Brightness', 0.75) %CG: 0 = min Brightness and 1 = max Brightness of cell
addParameter(p, 'setBrightness', 1) %CG if setBrightness = 1, then CellBones will not reduce the brightness before the user
%decides which green lines to keep. Note that keeping high brightness of
%the cell during this step can make it hard to see some green lines
addParameter(p,'Resize', 'On') %CG: should always be kept on. If there is no resizing the resolution 
%is divided by 1 and nothing changes. If resize is switched to 'off'
%without other modifications being made, errors could occur. 
addParameter(p,'Overwrite', 'Off') %CG: instead of deleting CellBonesCA file
addParameter(p,'AddSPorWL','Off') %CG: see 'Inputs' section below for description
addParameter(p,'DeleteGreenLines','Off') %CG: see 'Inputs' section below for description
addParameter(p,'SkipDistanceMeasurement','Off') %CG: see 'Inputs' section below for description

addParameter(p,'ROIColor',[1 0 0])
addParameter(p,'StartPointColor',[1 0 1])

addParameter(p,'VD','Off')
%CG: VD = verify distance; helps the user to test that CellBones is
%measuring distances correctly. 
 
parse(p,varargin{:}); param = p.Results; 

% IDString = param.IDString; 
suffixStrs = param.suffixStrs;

datadir = param.datadir; MaxDistance = param.MaxDistance; Resolution = param.Resolution;
Brightness = param.Brightness; setBrightness = param.setBrightness; 
Overwrite = param.Overwrite; Resize = param.Resize;
AddSPorWL = param.AddSPorWL; DeleteGreenLines = param.DeleteGreenLines; 
SkipDistanceMeasurement = param.SkipDistanceMeasurement; VD = param.VD;

ROIColor = param.ROIColor; StartPointColor = param.StartPointColor;

if strcmp(Overwrite, 'On'); Overwrite = 1; disp('data is being overwritten!')
elseif ~strcmp(Overwrite, 'On'); Overwrite = 0; end
    
if strcmp(AddSPorWL, 'On') || strcmp(DeleteGreenLines, 'On') || strcmp(SkipDistanceMeasurement, 'On')
    Overwrite = 1;
%CG: with Overwrite enabled, data in the CellBonesCAs file can be modified.
end

%CG: Some of the input arguments are hierarchical. It is assumed that if
%the user sets DeleteGreenLines to 'On', start points and whitelines
%already exist and the user wants these to appear as if s/he had entered
%the input argument: AddSPorWL, 'On'. 
% if strcmp(DeleteGreenLines, 'On'); AddSPorWL = 'On'; end
% if strcmp(SkipDistanceMeasurement, 'On'); AddSPorWL = 'On'; DeleteGreenLines = 'On'; end


if Brightness > 1 || Brightness < 0; errordlg('adjust brightness with input argument 1 using values between 0 and 1'); end
if isempty(Brightness); Brightness = 1; end
Brightness = 1 - Brightness;
%Contrast can be adjusted to improve visualization of dendrites. very thin
%dendrites can be hard to see if Brightness is set too high. range: 0<Brightness<=1

% if strcmp(Resize, 'On'); Resolution = Resolution/3; end

%Inputs:
%   Property                Values
%   --------                -------------------------------------------------------
%   AddSPorWL               'On' or 'Off'. Enables addition of new start points and subsequent
%                           white lines. Click 'a' to accept current start point 
%                           configuration and begin drawing white lines. Right-
%                           click to complete white line. Right-click again to
%                           accept current white line configuration. Default:
%                           'Off'.
%
%   DeleteGreenLines        'On or 'Off'. If 'On' then AddSPorWL = 'On'. After accepting start point and 
%                           white line configuration, green lines are built/rebuilt
%                           (any pre-existing green line configurations are not
%                           saved) and the user then deletes green lines as
%                           desired. Default: 'Off'.
%   
%   SkipDistanceMeasurement 'On' or 'Off'. If 'On' then AddSPorWL &
%                           DeleteGreenLines = 'On'. CellBones will use
%                           saved distance measurements from the
%                           CA_Distance variable. With SkipDistanceMeasurement
%                           set to 'On' the user will only be able to
%                           change which ROIs are excluded from further
%                           analysis.

%Outputs:
%Cell array contents (by column):   
%CA_Bones{:,1} - line handle for each white line
%CA_Bones{:,2} - bend points created by user
%CA_Bones{:,3} - all points on the line
%CA_Bones{:,4} - index values for all connected ROIs in CA_ROIc. This
%allows any ROIs connected to the current white line to be readily found in
%CA_ROIc.

%CA_ROIc contains all info about ROIs can their connectivity to white
%lines. Details below:
%CA_ROIc{:,1} - ROI linear indices taken from CA_ROIs
%CA_ROIc{:,2} - ROI patch handle
%CA_ROIc{:,3} - ROI index value in CA_ROIs
%CA_ROIc{:,4} - proximal points. Points on white lines closest to ROI in
%current row. x and y coordinates in first 2 columns. Coordinates ordered
%according to the order of white lines in CA_Bones.
%CA_ROIc{:,5} - ROIc_xy; x and y coordinates for every point along each
%green line. Will be empty if ROI is too far from any line or is too close
%to the edge of the field of view. 
%CA_ROIc{:,6} - line handle for connecting (green) line.
%CA_ROIc{:,7} - distance along the green line (in micrometers).
%CA_ROIc{:,8} - coordinates for the connecting white lines (in connecting
%order towards the connected start point). This enables the connected
%white lines in CA_Bones to be readily identified programmatically.
%CA_ROIc{:,9} - mircon distance between each line (green or white)
%intersection and the connecting start point. This information is used to
%selectively display fluorescence intensity data for each ROI that is
%closer to the connecting start point (along connecting dendrites) than the
%current ROI.
%CA_ROIc{:,10} - total distance between ROI and start point.
%CA_ROIc{:,11} - diameter of branch as estimated by the
%'diameterExtrapolation' function.

screensize = get( groot, 'Screensize' );

if isempty(datadir)
    datadir = uigetdir('Select data folder for CellBones');
end
% datadir = '/Users/cgilbri/Documents/MATLAB/Test_CellBones/CG2705161';
% datadir = '/Users/cjgilbride/Documents/DecJan2016_17/MATLAB/Test_CellBones/CG2705161';
CBfilefound = 0; abortswitch = 0;
if ischar(datadir)

    Sz_suffixStrs = size(suffixStrs); [~,CellName,~] = fileparts(datadir); TDFolders = {};tdfc = 0;
    for cCell = 1 : Sz_suffixStrs(2)
        cItem = strcat(CellName,suffixStrs{cCell},'_ThresholdData');
        tdfc = tdfc + 1; TDFolders{tdfc,1} = cItem;
        cd(strcat(datadir,'/', cItem))
        fmedian = []; cfile = 0; dirInfo = dir;
        while isempty(fmedian) 
            cfile = cfile + 1; currentfile = dirInfo(cfile).name;
            if ~isempty(strfind(currentfile, 'Chunk_1'))
                data = load(currentfile); fmedian = data.fmedian;
                try pixelBinningFactor = data.pixelBinningFactor;
                catch
                    pixelBinningFactor = '3'; %errordlg('WARNING: default pixel binning factor of applied. This is intended for debugging purposes only')
                end
            end
        end
    end
    Resolution = Resolution/str2double(pixelBinningFactor);
    cd(datadir);
    addpath(genpath(datadir))
    CellBonesFile = strcat(CellName,'_','CellBonesCAs.mat');
    
    try testing = load(CellBonesFile); CBfilefound = 1;
    catch
        CBfilefound = 0;
    end
    
    if CBfilefound == 1 && Overwrite == 1; CBfilefound = 0; end
%CG: if a CellBones file is detected, but Overwrite is enabled then
%CBfilefound should be set to 0 to allow full function call. 
end

%CG: whiteROI, anotherROI and anotherColor can be used to change the
%colour of specifc ROIs (given their CA_ROIs index) for illustrative
%purposes. Set each variable to be an empty numeric array to disable this
%feature. 
whiteROI = 246; anotherROI = 8603; anotherColor = [0.25,0.25,1];

if ischar(datadir) && CBfilefound == 0 

    TDFolders_sorted = cell(tdfc, 1);
    NameGroupSz = tdfc/Sz_suffixStrs(2);
    for cCell = 1 : Sz_suffixStrs(2)
        cItem2sort = TDFolders{cCell};
        SortIdx = (NameGroupSz*(cCell-1))+1;
        TDFolders_sorted{SortIdx, 1} = cItem2sort;
    end
    
    [~, ExpName, ~] = fileparts(datadir);
%     IDStrSta = strfind(ExpName, IDString);
            
%     if IDStrSta > 1
%         ExpName = ExpName(IDStrSta:end);
%     end 

    ThresDataList = {};
    
    for cCell = 1 : Sz_suffixStrs(2)
        tdfcounter = 0;
        for cDFIdx = 1 : tdfc
            cDF = TDFolders_sorted{cDFIdx, 1};
            if strfind(cDF, strcat(ExpName, suffixStrs{cCell}, '_ThresholdData'))
                tdfcounter = tdfcounter + 1;
            end
        end
        ThresDataList{cCell, 1} = suffixStrs{cCell};
        ThresDataList{cCell, 2} = tdfcounter;
    end

    logFluoImage=log10(fmedian); [xdim, ydim]=size(fmedian);
    logFluoImage_copy = logFluoImage;

    fh = figure('Visible', 'Off', 'NumberTitle', 'Off', 'Name', 'Connecting ROIs to dendrites GUI');
    
    ax1 = axes; logFluoImage=logFluoImage-min(logFluoImage(:));  
    
    maxval = max(logFluoImage(:));
    logFluoImage_track = logFluoImage >= Brightness*maxval;
    logFluoImage(logFluoImage_track) = Brightness*maxval;

    logFluoImage=logFluoImage/max(logFluoImage(:));                             
    logFluoImage=uint8(logFluoImage*256);                                       
    cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);
    left = screensize(3)*0.2; bot = screensize(2);
    wid = screensize(3)*0.66; hei =  screensize(4);
    colormap(cyanColorMap); fmedianRGB=ind2rgb(logFluoImage,cyanColorMap);
    Sz_fmedianRGB = size(fmedianRGB); im1 = imagesc(fmedianRGB);
    set(ax1, 'XTickLabel', '', 'TickLength', [0 0]);
    set(ax1, 'YTickLabel', '', 'TickLength', [0 0]);
    set(ax1, 'box', 'off'); set(fh, 'OuterPosition', [left, bot, wid, hei]);
    FigPos = fh.Position;
    set(fh, 'NumberTitle', 'Off')
    axis image; 
    set(fh, 'Visible', 'On');
   
    cd(strcat(datadir, '/', TDFolders_sorted{end}))
    dirInfo = dir; NumItems = size(dirInfo,1); cItem = 0;

    while cItem < NumItems
        cItem = cItem + 1; currentItem = dirInfo(cItem).name;
        Sz_currentItem = size(currentItem);
        if Sz_currentItem(2) > 7
            if strcmp(currentItem(1:8), 'CAs_Core')
                uScoreIdces = strfind(currentItem, '_');
                filekey = currentItem(1:uScoreIdces(2));
                cItem = NumItems;
            elseif strcmp(currentItem(1:8), 'CAs_CoreFour')
                uScoreIdces = strfind(currentItem, '_');
                filekey = currentItem(1:uScoreIdces(2));
                cItem = NumItems;
            end
        end
    end
        
    load(strcat(filekey, ExpName, suffixStrs{end}), 'CA_ROIs')
    Sz_CA_ROIs = size(CA_ROIs);
    if numel(Sz_CA_ROIs) == 2
        TotalCells = Sz_CA_ROIs(1)*Sz_CA_ROIs(2);
    else
        TotalCells = Sz_CA_ROIs(1)*Sz_CA_ROIs(2)*Sz_CA_ROIs(3);
    end
    cFullCell = 0; 
    mySpecialIndices = [1519,246,4487,1530,6554,8603];
    for cCell = 1 : TotalCells
        if ~isempty(CA_ROIs{cCell}) && ~ischar(CA_ROIs{cCell}) 
%CG: [DEPRECATED] subsidary ROIs will fill cells in CA_ROIs with the string 'subsid'.
%CG: A CA_ROIc might have already been saved. The user may want to redo
%CA_ROIc or add to existing data in CA_ROIc. For now, a new CA_ROIc is
%created. It might be replaced by or added to existing data later on. 
            cFullCell = cFullCell + 1;
            [Row, Col, ~] = ind2sub(size(fmedian), CA_ROIs{cCell});
            outlineVals = convhull(Row, Col);
%             ROIpatch = patch(Col(outlineVals), Row(outlineVals),...
%                     ROIColor, 'EdgeColor', ROIColor, 'Visible', 'On');
%             CA_ROIc{cFullCell,1} = CA_ROIs{cCell};
%             CA_ROIc{cFullCell,2} = ROIpatch;
%             CA_ROIc{cFullCell,3} = cCell;
            if ~isempty(anotherROI) && ~isempty(anotherColor) 
                if cCell == anotherROI
                    ROIpatch = patch(Col(outlineVals), Row(outlineVals),...
                        anotherColor, 'EdgeColor', anotherColor, 'Visible', 'On');
                elseif cCell == whiteROI
                    ROIpatch = patch(Col(outlineVals), Row(outlineVals),...
                        [1,1,1], 'EdgeColor', [1,1,1], 'Visible', 'On');
                else
                    ROIpatch = patch(Col(outlineVals), Row(outlineVals),...
                        ROIColor, 'EdgeColor', ROIColor, 'Visible', 'On');
                end
            else
                ROIpatch = patch(Col(outlineVals), Row(outlineVals),...
                        ROIColor, 'EdgeColor', ROIColor, 'Visible', 'On');
            end
            CA_ROIc{cFullCell,1} = CA_ROIs{cCell};
            CA_ROIc{cFullCell,2} = ROIpatch;
            CA_ROIc{cFullCell,3} = cCell;
        end
    end

    seltype  = ''; ExclLoop = 1; zone = 1; linehandle_error = 'noerror'; lastAction = '';
    while ExclLoop == 1 && abortswitch == 0
        if zone == 1
            k = [];
            if strcmp(AddSPorWL, 'On') || (strcmp(AddSPorWL, 'Off') &&...
                    (strcmp(DeleteGreenLines,'On') || strcmp(SkipDistanceMeasurement,'On')))
                skeleton = load(CellBonesFile);
                CA_StartPoints = skeleton.CA_StartPoints; 
                Sz_CA_StartPoints = size(CA_StartPoints);
                for cCell = 1 : Sz_CA_StartPoints(1)
                    set(CA_StartPoints{cCell}, 'Parent', ax1)
                    set(CA_StartPoints{cCell}, 'Visible', 'On')
                end
            else 
                CA_StartPoints = {};
            end

            while isempty(k) && abortswitch == 0
%CG: while-loop controlling start point placement. Break out of this loop
%by pressing 'a' key, before returning to first loop and proceeding with
%line building. 
                title(sprintf(strcat('Instructions (Step 1): Left-click to set start points. Right-click to delete previous.\n',...
                'Press "a" key to accept and move on to tracing dendrites')), 'FontSize', 16)
                try
                    if (strcmp(DeleteGreenLines, 'On') || strcmp(SkipDistanceMeasurement,'On')) &&...
                            strcmp(AddSPorWL,'Off')
                        k = 1;
                    elseif strcmp(DeleteGreenLines, 'Off') && strcmp(SkipDistanceMeasurement,'Off')
                        k = waitforbuttonpress;
                    end
                    if k == 1
                        if (strcmp(DeleteGreenLines, 'On') || strcmp(SkipDistanceMeasurement,'On')) &&...
                                strcmp(AddSPorWL,'Off')
                            KeyPressRes = 'a';
%CG: jump to next step if the user wants to DeleteGreenLines or wants to
%repeat distance measurement while explicitly commanding that no new start
%points or white lines should be added. 
                        elseif strcmp(DeleteGreenLines, 'Off') && strcmp(SkipDistanceMeasurement,'Off')
                            KeyPressRes = get(fh,'CurrentCharacter');
                        end
                        if ~strcmp(KeyPressRes, 'a')
                            disp('press "a" key to continue or "z" to remove previous start point...')
                            disp('...otherwise left-click to set start points')
                            k = [];
                        elseif strcmp(KeyPressRes, 'a') && isempty(CA_StartPoints)
                            disp('set start points to begin tracing the dendritic tree')
                        end
                    elseif k == 0
                        mouseclick = get(fh, 'SelectionType');
                        if strcmp(mouseclick, 'normal')
                            cObj = gco; ClassDef = class(cObj);
                            if strcmp(ClassDef, 'matlab.graphics.primitive.Image')
                                cp = get(gca, 'CurrentPoint'); x_co = round(cp(1,2));
                                y_co = round(cp(2,1));

                                if isempty(CA_StartPoints)
                                    cPatch = 1;
                                else
                                    Sz_CA_StartPoints = size(CA_StartPoints);
                                    cPatch = Sz_CA_StartPoints(1)+1;
                                end
                                CA_StartPoints{cPatch,1} = patch(y_co, x_co, StartPointColor,...
                                    'Marker', 'x', 'MarkerSize', 18, 'MarkerFaceColor', StartPointColor,...
                                    'MarkerEdgeColor', StartPointColor, 'LineWidth', 2, 'Visible', 'On');

                            elseif strcmp(ClassDef, 'matlab.graphics.primitive.Patch')
                                x_co = get(cObj, 'XData'); y_co = get(cObj, 'YData');

                                Sz_CA_StartPoints = size(CA_StartPoints); cPatch = 0;

                                while cPatch < Sz_CA_StartPoints(1)
                                    cPatch = cPatch + 1;
                                    CheckX = get(CA_StartPoints{cPatch}, 'XData') - x_co;
                                    CheckY = get(CA_StartPoints{cPatch}, 'YData') - y_co;
                                    if sum(CheckX) == 0 && sum(CheckY) == 0
                                        set(CA_StartPoints{cPatch,1}, 'Visible', 'Off')
                                        if Sz_CA_StartPoints(1) == 1
                                            CA_StartPoints = {};
                                        else
                                            CA_StartPoints(cPatch,:) = [];
                                            cPatch = Sz_CA_StartPoints(1);
                                        end
                                    end
                                end
                            end
                        elseif strcmp(mouseclick, 'alt')

                            set(CA_StartPoints{end,1}, 'Visible', 'Off')
                            Sz_CA_StartPoints = size(CA_StartPoints);
                            if Sz_CA_StartPoints(1) == 1
                                CA_StartPoints = {};
                            else
                                CA_StartPoints(end,:) = [];
                            end

                        end
                        k = [];
                    end
                catch 
%CG: in case the figure window is closed, then abortswitch = 1, which will
%stop CellBones without a crash. 
                    abortswitch = 1;
                end
            end
            if abortswitch == 0
                cd(datadir)

                try testing = load(CellBonesFile, 'CA_Bones'); save(CellBonesFile, 'CA_StartPoints', '-append');
%CG: If there is any saved 'CA_Bones' variable, CellBones append
%CA_StartPoints rather than potentially make an undesired overwrite.
                catch
                    save(CellBonesFile, 'CA_StartPoints'); 
                end
                
                if (strcmp(AddSPorWL, 'Off') && strcmp(DeleteGreenLines, 'Off') &&...
                        strcmp(SkipDistanceMeasurement,'Off'))
                    CA_Bones = {};
                elseif (strcmp(AddSPorWL, 'On') && strcmp(DeleteGreenLines, 'Off') &&...
                        strcmp(SkipDistanceMeasurement,'Off'))
                    try CA_Bones = skeleton.CA_Bones; 
                        CA_Bones = skeleton.CA_Bones; 
                        Sz_CA_Bones = size(CA_Bones);
                        for cBoneIdx = 1 : Sz_CA_Bones(1)
                            set(CA_Bones{cBoneIdx, 1}, 'Parent', ax1)
                            set(CA_Bones{cBoneIdx, 1}, 'Visible', 'on')
                            set(CA_Bones{cBoneIdx, 1}, 'LineWidth', 2)
                            set(CA_Bones{cBoneIdx, 1}, 'Color', 'w')
                        end
                    catch
                        CA_Bones = {};
                    end
                    
                elseif strcmp(AddSPorWL, 'On') || ((strcmp(DeleteGreenLines, 'On') || strcmp(SkipDistanceMeasurement,'On')) &&...
                        strcmp(AddSPorWL,'Off'))
                    CA_Bones = skeleton.CA_Bones; 
                    Sz_CA_Bones = size(CA_Bones);
                    for cBoneIdx = 1 : Sz_CA_Bones(1)
                        set(CA_Bones{cBoneIdx, 1}, 'Parent', ax1)
                        set(CA_Bones{cBoneIdx, 1}, 'Visible', 'on')
                        set(CA_Bones{cBoneIdx, 1}, 'LineWidth', 2)
                        set(CA_Bones{cBoneIdx, 1}, 'Color', 'w')
                    end
                end
                if strcmp(AddSPorWL, 'On') || (strcmp(AddSPorWL, 'Off') &&  strcmp(DeleteGreenLines, 'Off') &&...
                        strcmp(SkipDistanceMeasurement,'Off'))
                    title(sprintf(strcat('Instructions (Step 2): Left-click on start point or white line to begin tracing a dendrite.\n',...
                    'Build line with further left-clicks. Right-click to stop building current white line.\n',...
                    'Left-click on start point or white line again to build new white line.\n',...
                    'Alternatively, press "a" to accept and move to next step')), 'FontSize', 16)

                    axishandle = gca; 
                    [xPoly, yPoly, linehandle, linehandle_error, seltype, abortswitch, lastAction] = getline_CG_CellBones(axishandle, 0,...
                        CA_Bones, seltype, StartPointColor,linehandle_error,abortswitch, lastAction, fh);

                    if ~strcmp(linehandle_error, 'error1') &&...
                        ~strcmp(linehandle_error, 'error2') &&...
                            ~isempty(xPoly) && (~strcmp(seltype, 'finish') ||...
                            strcmp(seltype, 'finish') && strcmp(lastAction, 'a-press-but-no-alt-click')) &&...
                                ~isempty(linehandle) && abortswitch == 0
%CG: if strcmp(seltype, 'finish') && strcmp(lastAction,
%'a-press-but-no-alt-click'), then the user has pressed a without closing
%the final line first. 
%
                        if isempty(CA_Bones)
                            push = 0;
                        else
                            Sz_CA_Bones = size(CA_Bones);
                            push = Sz_CA_Bones(1);
                        end

                        CA_Bones{push+1,1} = linehandle; CA_Bones{end,2} = [xPoly yPoly];
                        yRange = max(yPoly) - min(yPoly); xRange = max(xPoly) - min(xPoly);
                        if xRange > yRange
                            if xPoly(1) < xPoly(end)
                                AllXorY = xPoly(1):xPoly(end);
                            else
                                AllXorY = xPoly(end):xPoly(1);
                            end
                        else
                            if yPoly(1) < yPoly(end)
                                AllXorY = yPoly(1):yPoly(end);
                            else
                                AllXorY = yPoly(end):yPoly(1);
                            end
                        end                        

                        t = 0:(1/(numel(AllXorY))):1; 
                        if numel(xPoly) > 1 && numel(yPoly) > 1 
                            Allxy = interparc(t, xPoly, yPoly);
                            Allxy = round(Allxy);
                        else 
                            Allxy = [];
                        end
                        
                        if isempty(Allxy)
                            CA_Bones(end, :) = [];
                        elseif ~isempty(Allxy)   
                            CA_Bones{end,3} = Allxy;
                        end
                    elseif strcmp(linehandle_error, 'error1')
                        title(sprintf(strcat('line must start at the end of an existing line or on a start point.\n',...
                            'Or the user must press the "a" key to proceed to the next step')), 'FontSize', 16)
                        
                    elseif strcmp(linehandle_error, 'error2')
                        title(sprintf(strcat('finish building last white line (by right-clicking) before pressing "a"')), 'FontSize', 16)
                    elseif strcmp(seltype, 'finish')
                        ExclLoop = 0;
                    end
                end
            end
        elseif zone > 1 && abortswitch == 0
            axishandle = gca;

            [xPoly, yPoly, linehandle, linehandle_error, seltype, abortswitch, lastAction] = getline_CG_CellBones(axishandle, 0,...
                CA_Bones, seltype, StartPointColor,linehandle_error,abortswitch, lastAction, fh);

            if ~strcmp(linehandle_error, 'error1') &&...
                ~strcmp(linehandle_error, 'error2') &&...
                    ~isempty(xPoly) && (~strcmp(seltype, 'finish') ||...
                        strcmp(seltype, 'finish') && strcmp(lastAction, 'a-press-but-no-alt-click')) &&...
                            ~isempty(linehandle)

                CA_Bones{end+1,1} = linehandle; CA_Bones{end,2} = [xPoly yPoly];
                yRange = max(yPoly) - min(yPoly); xRange = max(xPoly) - min(xPoly);

                if xRange > yRange
                    if xPoly(1) < xPoly(end)
                        AllXorY = xPoly(1):xPoly(end);
                    else
                        AllXorY = xPoly(end):xPoly(1);
                    end
                else
                    if yPoly(1) < yPoly(end)
                        AllXorY = yPoly(1):yPoly(end);
                    else
                        AllXorY = yPoly(end):yPoly(1);
                    end
                end

                t = 0:(1/(numel(AllXorY))):1; 
                Allxy = interparc(t, xPoly, yPoly);
                Allxy = round(Allxy);
                if isempty(Allxy)
                    CA_Bones(end, :) = [];
                elseif ~isempty(Allxy)   
                    CA_Bones{end,3} = Allxy;
                end
            elseif strcmp(linehandle_error, 'error1')
                title(sprintf(strcat('line must start at the end of an existing line or on a start point.\n',...
                    'Or the user must press the "a" key to proceed to the next step')), 'FontSize', 16)
                
            elseif strcmp(linehandle_error, 'error2')
                title(sprintf(strcat('finish building last white line (by right-clicking) before pressing "a"')), 'FontSize', 16)
            elseif strcmp(seltype, 'finish')
                ExclLoop = 0;
            end                
        end

        if zone >= 2 && ~strcmp(linehandle, 'error1')
            zone = zone + 1; 
        else
            zone = zone + 1;
        end
        
        if strcmp(AddSPorWL, 'Off') && (strcmp(DeleteGreenLines, 'On') ||...
                strcmp(SkipDistanceMeasurement,'On'))
            ExclLoop = 0;
        end
        if strcmp(lastAction,'a')
            ExclLoop = 0;
        end

    end
    if abortswitch == 0
        cd(datadir)
        save(CellBonesFile, 'CA_Bones', 'Resize', '-append')

%CG: User has set a line along all the dendrites s/he wants to include in
%the analysis. CellBones will now join the centre points of each ROI to the
%most appropriate nearest point on one of these lines. 

%CG: Initially, CellBones identifies the top NumProximalPnts for each ROI,
%with each proximal point being located on a separate line. Lines that are
%more than MaxDistance are excluded from this list. So NumProximalPnts is a
%maximal value as well rather than the number of proximal points the user
%can expect to have in the list - it depends on the number of independent
%lines the user has drawn. The chance that a ROI is truly connected to a
%dendrite with a line more than 20um away is extremely low. MaxDistance
%helps to reduce computations, by removing computations that are, beyond
%reasonable doubt, completely pointless. 

%CG: for each ROI, calculate the distance (in ?m) between the centre point
%of the ROI and each point along each line within range. Pick the closest
%point along each line. Exclude any points that are out of range before
%carrying out further computations.

        if (strcmp(SkipDistanceMeasurement, 'On') || strcmp(DeleteGreenLines, 'On')) && strcmp(AddSPorWL,'Off')
            NumPatches = size(CA_ROIc,1);
            for cROIidx = 1 : NumPatches
                set(CA_ROIc{cROIidx,2},'Visible', 'Off')
            end

%             CA_ROIc = skeleton.CA_ROIc;
            CA_ROIc = cat(2, CA_ROIc, skeleton.CA_ROIc(:,4:end));
            Sz_CA_ROIc = size(CA_ROIc);
            for cROIidx = 1 : Sz_CA_ROIc(1)
                set(CA_ROIc{cROIidx, 2}, 'Parent', ax1); 
                set(CA_ROIc{cROIidx, 2}, 'Visible', 'On')
                set(CA_ROIc{cROIidx, 6}, 'Parent', ax1);
            end
            

        elseif strcmp(DeleteGreenLines, 'Off') 

            Sz_CA_Bones = size(CA_Bones); NumBones = Sz_CA_Bones(1);
            Sz_CA_ROIc = size(CA_ROIc); TotalNumROIs = Sz_CA_ROIc(1);

            for cROIidx = 1 : TotalNumROIs

                cROI = CA_ROIc{cROIidx, 1}; cROI_cp = cROI(round(numel(cROI)/2));
%CG: extract centre point of current ROI. 
                [y_ROI_cp, x_ROI_cp] = ind2sub([xdim, ydim], cROI_cp);

                x_proximals = []; y_proximals = [];
                for cBoneidx = 1 : NumBones
                    cBone = CA_Bones{cBoneidx, 3}; Sz_cBone = size(cBone);
                    cBoneDistances = [];

                    for cpnt = 1 : Sz_cBone(1)
                        x_bone = cBone(cpnt,1); y_bone = cBone(cpnt,2);
                        x_ROI_cp_diff = x_bone - x_ROI_cp;
                        y_ROI_cp_diff = y_bone - y_ROI_cp;
                        hyp_d = hypot(x_ROI_cp_diff, y_ROI_cp_diff);
                        real_d = hyp_d/Resolution;

                        if real_d <= MaxDistance
                            if isempty(cBoneDistances)
                                cBoneDistances(1,1) = real_d; cBoneDistances(1,2) = cpnt;
                            else
                                cBoneDistances(end+1,1) = real_d; cBoneDistances(end,2) = cpnt;
                            end
                        end
                    end
                    if ~isempty(cBoneDistances)
                        [~, MinIdx] = min(cBoneDistances(:,1));
                        if isempty(x_proximals)
                            y_proximals = cBone(cBoneDistances(MinIdx, 2),1);
                            x_proximals = cBone(cBoneDistances(MinIdx, 2),2);
                        else
                            y_proximals(end+1,1) = cBone(cBoneDistances(MinIdx, 2),1);
                            x_proximals(end+1,1) = cBone(cBoneDistances(MinIdx, 2),2);
                        end
                    end

%CG: add the closest points along the current line to the appropriate row
%in column 4 of the CA_Bones cell array.

%CG: note that the order in which "bones" are entered into the
%subCA_proximals cell array are the same order they appear in the CA_Bones
%cell array. 

                end
                CA_ROIc{cROIidx, 4} = [x_proximals y_proximals];
            end

%CG: now it is necessary to connect the centre pixel of each ROI with the 
%nearest points on the bones. Find the neighbouring pixels of the current
%pixel (this is the centre pixel at start) and measure only the ones closer
%to the target pixel (target pixel = current proximal point) for
%fluorescence intensity (this is the intensity value seen in the figure
%image - fmedianRGB). Of the closer pixels, select the one with the highest
%fluorescence intensity to be the current pixel. 

            neighbors_sl = 3; FirstROI = 0; 
            
            wb = waitbar(0, ''); wbPos = wb.Position;
            wb_left = FigPos(1)+FigPos(3)*0.18; wb_bot = FigPos(2); 
            wb_width = wbPos(3); wb_height = wbPos(4);
            set(wb, 'Position', [wb_left wb_bot wb_width wb_height]);
            
            for cROIidx = 1 : TotalNumROIs
                xyProximals = CA_ROIc{cROIidx, 4}; LastDistance = []; TotalDistance = 0;
                if ~isempty(xyProximals)
                    Sz_CA_Bones = size(CA_Bones); NumBonesInc = Sz_CA_Bones(1);
                    cROI = CA_ROIc{cROIidx, 1}; 
                    currpix = cROI(round(numel(cROI)/2));
                    [x_currpix, y_currpix] = ind2sub([xdim, ydim], currpix);
                    FirstROI = FirstROI + 1; LineCoos = [];
                    yLineCoos = []; xLineCoos = [];
                    for cBoneIdx = 1 : NumBonesInc

                        if isempty(LineCoos)
                            LineCoos = CA_Bones{cBoneIdx, 3};
                            yLineCoos = LineCoos(:,1);
                            xLineCoos = LineCoos(:,2);
                        else
                            LineCoos_temp = CA_Bones{cBoneIdx, 3}; 
                            LineCoos = cat(1,LineCoos, LineCoos_temp);
                            yLineCoos_temp = LineCoos(:,1);
                            yLineCoos = cat(1,yLineCoos, yLineCoos_temp);
                            xLineCoos_temp = LineCoos(:,2);
                            xLineCoos = cat(1,xLineCoos, xLineCoos_temp);
                        end

                    end
                    x_tgts = xLineCoos; y_tgts = yLineCoos; connectedIdces = [];
                    connX = []; connY = []; lengthcount = 0; ROIc_xy = [];
                    cIdx = 0; twol = 0;
                    while isempty(find(twol == 2,1))
                        x_CNhoodDiag = x_currpix-floor(neighbors_sl/2):1:x_currpix+floor(neighbors_sl/2);
                        y_CNhoodDiag = y_currpix-floor(neighbors_sl/2):1:y_currpix+floor(neighbors_sl/2);
%CG: find the neighbors of the central point (just like in AARG), except
%only 3x3. 
                        for cEl = 1 : neighbors_sl
                             x_CNhood(cEl:neighbors_sl:neighbors_sl^2) = x_CNhoodDiag(cEl);
                        end

                        for cEl = 1 : neighbors_sl
                            if cEl == 1
                                y_CNhood(cEl:1:neighbors_sl) = y_CNhoodDiag(cEl);
                            elseif cEl > 1
                                y_CNhood(((cEl-1)*neighbors_sl)+1:1:neighbors_sl*cEl) = y_CNhoodDiag(cEl);
                            end
                        end
                        if ~isempty(find(x_CNhood < 1, 1)) || ~isempty(find(y_CNhood < 1, 1)) ||...
                                ~isempty(find(x_CNhood > xdim, 1)) || ~isempty(find(y_CNhood > ydim, 1))
%CG: if the line hits an edge before reaching the white line, give no green
%line for this ROI. 
                            ROIc_xy = []; twol = 2;
                            break
                        end
                        CA_Minis = {}; currpix_ds_opts = [];

                        for cBoneIdx = 1 : NumBonesInc

                            LineCoos = CA_Bones{cBoneIdx, 3};
                            yLineCoos = LineCoos(:,1); xLineCoos = LineCoos(:,2);

                            iMinDists = []; iMinDistIdces = [];
                            for cEl = 1 : numel(x_CNhood)
                                cX = x_CNhood(cEl); cY = y_CNhood(cEl);

                                cXs = ones(numel(xLineCoos), 1);
                                cXs = cXs.*cX;
                                cYs = ones(numel(yLineCoos), 1);
                                cYs = cYs.*cY;

                                Xs = xLineCoos - cXs; Ys = yLineCoos - cYs;
                                hs = hypot(Xs, Ys);

                                [sMinDists, sMinDistIdces] = min(hs);
                                if isempty(iMinDists)
                                    iMinDists = sMinDists;
                                    iMinDistIdces = sMinDistIdces;
                                else
                                    iMinDists(end+1) = sMinDists;
                                    iMinDistIdces(end+1) = sMinDistIdces;
                                end

                            end

                            cXs = x_currpix - xLineCoos(iMinDistIdces(round(numel(x_CNhood)/2)));
                            cYs = y_currpix - yLineCoos(iMinDistIdces(round(numel(y_CNhood)/2)));

                            if isempty(currpix_ds_opts)
                                currpix_ds_opts = hypot(cXs, cYs);
                            else
                                currpix_ds_opts(end+1) = hypot(cXs, cYs);
                            end
                            if isempty(CA_Minis)
                                CA_Minis{1,1} = iMinDists;
                            else
                                CA_Minis{end+1,1} = iMinDists;
                            end
                        end

                        [currpix_ds, currpix_mindx] = min(currpix_ds_opts);

                        if isempty(LastDistance)
                            LastDistance = currpix_ds;
                        else
                            ShiftDistance = LastDistance - currpix_ds;
                            TotalDistance = TotalDistance + ShiftDistance;
                            LastDistance = currpix_ds;
                        end

                        iMinDists = CA_Minis{currpix_mindx,1};

                        if currpix_ds ~= 0 

                            x_CNhood(round(numel(x_CNhood)/2)) = [];
                            y_CNhood(round(numel(x_CNhood)/2)) = [];
                            iMinDists(round(numel(x_CNhood)/2)) = [];
                            iMinDistIdces(round(numel(x_CNhood)/2)) = [];
%CG: remove the current point from the arrays. Need only the neighbors now.
                            All_CNhood = sub2ind([xdim ydim], x_CNhood, y_CNhood);

                            if ~isempty(connectedIdces)

                                preidx = sub2ind([xdim ydim], connX, connY);
                                preidxlog = All_CNhood == preidx;

                                x_CNhood(preidxlog) = []; y_CNhood(preidxlog) = [];
                                iMinDists(preidxlog) = []; iMinDistIdces(preidxlog) = [];
                            end
                            pIdces = find(iMinDists < currpix_ds);
                            if currpix_ds == 1 && min(iMinDists) == 1
%CG: if the current pixel is very close to the line, it can happen that
%pIdces will fall empty because both min(MinDists) = 1 and currpix_ds = 1.
%It is better not to change pIdces = find(MinDists < currpix_ds) to pIdces
%= find(MinDists <= currpix_ds) because this will cause the connecting line
%to meander towards the dendrite line. When pIdces is empty in this
%instance, it should only be when the current pixel is almost touching the
%dendritic line and all neighbors with 1 should be equal distance to the
%line. 
                                pIdces = find(iMinDists == 1);
                            end
                            valid_CNhood = sub2ind([xdim ydim], x_CNhood(pIdces),y_CNhood(pIdces));
%CG: get the maximum fluorescence intensity value of all eligible pixels. 
%The one with the highest intensity will be the next current pixel. 

                            connX = x_currpix; connY = y_currpix;

                            [~, MaxIdx] = max(fmedian(valid_CNhood));

                            x_currpix = x_CNhood(pIdces(MaxIdx));
                            y_currpix = y_CNhood(pIdces(MaxIdx));

                            if isempty(ROIc_xy)
                                ROIc_xy = [x_currpix y_currpix];
                            else
                                ROIc_xy(end+1,1) = x_currpix;
                                ROIc_xy(end,2) = y_currpix;
                            end
                            lengthcount = lengthcount + 1;

                            connectedIdces = cIdx;
                            cIdx = sub2ind([xdim, ydim], x_currpix, y_currpix);

                            xl = x_tgts == x_currpix; yl = y_tgts == y_currpix;
                            twol = xl + yl;

                        else
                            cIdx = sub2ind([xdim, ydim], x_currpix, y_currpix);
                            xl = x_tgts == x_currpix; yl = y_tgts == y_currpix;
                            twol = xl + yl;
                        end

                        if isempty(cIdx); cIdx = 0; end
                    end

                    if lengthcount == 0 
%CG: in this case, the ROI has most likely fallen on the dendritic line.
                       ROIc_xy = [x_currpix y_currpix];
                    end

                    if ~isempty(ROIc_xy)
%CG: ROIc_xy can still be empty at this point if the connecting line
%touches the edge of the image before reaching the target dendritic line. 
                        CA_ROIc{cROIidx, 5} = ROIc_xy;
                        xpts = ROIc_xy(:,1); 
                        ypts = ROIc_xy(:,2);
                        lh = line(ypts,xpts,'Color',[0 1 0],'LineWidth',2);
                        CA_ROIc{cROIidx, 6} = lh; 
                        CA_ROIc{cROIidx, 7} = (TotalDistance/Resolution);
                        set(CA_ROIc{cROIidx, 6}, 'Parent', ax1);
                    end

                end
                waitbar(cROIidx/TotalNumROIs, wb)  
            end
            close(wb); cd(datadir); save(CellBonesFile, 'CA_ROIc', '-append')
        end
    end
%CG: refine the lines. the user can exclude green lines that are not
%appropriately formed. This can happen especially when the ROIs are
%clustering.
    if abortswitch == 0
        unmodified_CA_ROIc = CA_ROIc; Sz_CA_ROIc = size(CA_ROIc);
        title(sprintf(strcat('Instructions (Step 3): Connections between ROIs and white lines have been automatically created.\n',...
            'These connections are represented by green lines. The user can select which green lines to discard by left-\n',...
            'clicking on the target line. Left-clicking once turns the line grey. Left-click on grey lines to delete the line.\n',...
            'Note that grey lines do not represent exclusions. Right-click to undo. Press "k" to delete all grey lines.\n',...
            'Press "a" to accept and move on to next stage.')), 'FontSize', 16)
        k = []; CA_LastAction = {};

        if setBrightness == 1
%CG: use an image with lower brightess to make the green lines easier to
%see.
            im2 = imagesc('CData',fmedianRGB);
            axis image;
        else
            logFluoImage = logFluoImage_copy;
            logFluoImage=logFluoImage-min(logFluoImage(:)); 
            logFluoImage=logFluoImage/max(logFluoImage(:));                             
            logFluoImage=uint8(logFluoImage*256);                                       
            cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);
            colormap(cyanColorMap); fmedianRGB=ind2rgb(logFluoImage,cyanColorMap);
            im2 = imagesc('CData',fmedianRGB);
            axis image;
        end
%CG: make sure all startpoints, lines etc... are on top.
        uistack(im2, 'bottom'); uistack(im1, 'bottom')

        while isempty(k) && abortswitch == 0
        
            try
                if strcmp(SkipDistanceMeasurement, 'On')
                    k = 1;
                elseif strcmp(SkipDistanceMeasurement, 'Off')
                    k = waitforbuttonpress;
                end
                if k == 1 
                    if strcmp(SkipDistanceMeasurement, 'On')
                        KeyPressRes = 'a';
                    elseif strcmp(SkipDistanceMeasurement, 'Off')    
                        KeyPressRes = get(fh,'CurrentCharacter'); pause(1)
                    end
                    if strcmp(KeyPressRes, 'a')
                        CA_LastAction = {};
                    elseif strcmp(KeyPressRes, 'k')
%CG: this action deletes all lines coloured grey. This gives the user the option
%to grey all the lines s/he might want to delete and then deleting them all
%in one go. 
                        for row = 1 : Sz_CA_ROIc(1)
                            cLine = CA_ROIc{row, 6};
                            if ~isempty(cLine)
                                cLineColor = get(cLine, 'Color');
                                GreyStatus = cLineColor - [0.5 0.5 0.5];
                                if sum(GreyStatus) == 0
                                    set(cLine, 'Visible', 'Off')
                                    CA_ROIc(row, 5) = {[]}; CA_ROIc(row, 6) = {[]};
                                    disp(strcat('line for ROI @ row "', num2str(row), '" deleted'))
                                end
                            end
                        end
                        k = [];
%CG: press 'a' or 'z' (or possibly 'y', which is more convenient than 'z'
%on a German keyboard). 'a' accepts everything as is and any modifications
%are saved. An unmodified version of the green lines is also saved. 'z'
%will undo the previous action. 
                    end

                elseif k == 0
                    mouseclick = get(fh, 'SelectionType');
                    if strcmp(mouseclick, 'normal')
                        cObj = gco; ClassDef = class(cObj); PointerStatus = get(fh, 'Pointer');
                        if strcmp(ClassDef, 'matlab.graphics.primitive.Line') && strcmp(PointerStatus, 'arrow')
                            cLineColor = get(gco, 'Color'); GreenStatus = cLineColor - [0 1 0];
                            GreyStatus = cLineColor - [0.5 0.5 0.5];

                            if sum(GreenStatus) == 0
%CG: current line has not been previously clicked and should be made grey.   
                                set(cObj, 'Color', [0.5 0.5 0.5])
                                xdata = get(cObj, 'XData'); ydata = get(cObj, 'YData');
                                Sz_xdata = size(xdata); co_numpoints = Sz_xdata(2);

                                for cROI = 1 : Sz_CA_ROIc(1)
                                    cLineVals = CA_ROIc{cROI, 5};
                                    Sz_cLineVals = size(cLineVals);
                                    cROI_numpoints = Sz_cLineVals(1);
                                    if ~isempty(cLineVals) && cROI_numpoints == co_numpoints
                                        xAns = cLineVals(:,1) - ydata';
                                        yAns = cLineVals(:,2) - xdata';
                                        if sum(xAns) == 0 && sum(yAns) == 0
                                            CA_LastAction{1,1} = cROI;
                                        end
                                    end
                                end
                            elseif sum(GreyStatus) == 0
%CG: current line has already been selected and should be removed. 
                                set(cObj, 'Visible', 'Off'); xdata = get(cObj, 'XData');
                                ydata = get(cObj, 'YData'); Sz_xdata = size(xdata);
                                co_numpoints = Sz_xdata(2);

                                for cROI = 1 : Sz_CA_ROIc(1)
                                    cLineVals = CA_ROIc{cROI, 5};
                                    Sz_cLineVals = size(cLineVals);
                                    cROI_numpoints = Sz_cLineVals(1);
                                    if ~isempty(cLineVals) && cROI_numpoints == co_numpoints
                                        xAns = cLineVals(:,1) - ydata';
                                        yAns = cLineVals(:,2) - xdata';

                                        if sum(xAns) == 0 && sum(yAns) == 0
                                            CA_LastAction{1,1} = cROI;
                                            CA_LastAction{2,1} = CA_ROIc{cROI, 5};
                                            CA_LastAction{3,1} = CA_ROIc{cROI, 6};

                                            CA_ROIc(cROI, 5) = {[]}; CA_ROIc(cROI, 6) = {[]};
                                            disp(strcat('line for ROI @ row "', num2str(cROI), '" deleted'))
                                        end                               
                                    end
                                end
                            end
                            k = [];
                        else
                            k = [];
                        end

                    elseif strcmp(mouseclick, 'alt') && strcmp(PointerStatus, 'arrow')
                        k = [];
%CG: undo previous action. If the last action was to delete a grey line, 
%then the line will be put back as a green line. If the last action was to
%make the line grey, it will be coloured green again.
                        if ~isempty(CA_LastAction)
                            Sz_CA_LastAction = size(CA_LastAction);
                            NumCells = Sz_CA_LastAction(1)*Sz_CA_LastAction(2);
                            if NumCells > 1
                                cROI = CA_LastAction{1,1};

                                CA_ROIc{cROI, 5} = CA_LastAction{2,1}; 
                                CA_ROIc{cROI, 6} = CA_LastAction{3,1};

                                set(CA_ROIc{cROI, 6}, 'Parent', ax1);
                                set(CA_ROIc{cROI, 6}, 'Color', [0 1 0]);
                                set(CA_ROIc{cROI, 6}, 'Visible', 'On');
                            elseif NumCells == 1
                                cROI = CA_LastAction{1,1};
                                cLineColor = get(CA_ROIc{cROI, 6}, 'Color');
                                GreenStatus = cLineColor - [0 1 0];
                                GreyStatus = cLineColor - [0.5 0.5 0.5];

                                if sum(GreenStatus) == 0
                                    set(CA_ROIc{CA_LastAction{1,1},6}, 'Color', [0.5 0.5 0.5])
                                elseif sum(GreyStatus) == 0
                                    set(CA_ROIc{CA_LastAction{1,1},6}, 'Color', [0 1 0])
                                end
                            end
                        end
                    else
                        k = [];
                    end
                else
                    k = [];
                end
            catch
                abortswitch = 1;
            end
        end
    
        cd(datadir)
        save(CellBonesFile, 'CA_ROIc', 'unmodified_CA_ROIc', '-append')
        Sz_CA_ROIc = size(CA_ROIc);
    end
    
    if abortswitch == 0
        if strcmp(SkipDistanceMeasurement, 'On')
            CA_Distances = skeleton.CA_Distances;
            for cROIidx = 1 : Sz_CA_ROIc(1)
                if isempty(CA_ROIc{cROIidx,6})
                    set(CA_ROIc{cROIidx,2}, 'FaceColor', [0 0 0])
                end
            end
        elseif strcmp(SkipDistanceMeasurement, 'Off')
%CG: connecting ROI lines have been refined. Now CellBones calculates the
%distance for each connected ROI to its connecting starting point. Once
%calculated, the distances are inserted into CA_Distances, which has the
%same format as CA_ROIs. Do this for each connected ROI in turn. 
            CA_Distances = cell(Sz_CA_ROIs);
%CG: Connection points need to be identified between the green lines and
%white lines, and then between white lines and another white line if the
%current white line is not connected directly to a start point. 
            NumROIs = Sz_CA_ROIc(1);
            title(sprintf(strcat('Connectivity accepted.',...
                '\n Measuring distances between ROIs and connected start points.\n',...
                'ROIs given a black spot will be excluded from analysis'),'FontSize', 16))
%CG: ROI on row 123 of CA_ROIc should be connected to Bone 8 (8th row
%CA_Bones). 
            wb = waitbar(0, ''); wbPos = wb.Position;
            wb_left = FigPos(1)+FigPos(3)*0.18; wb_bot = FigPos(2); 
            wb_width = wbPos(3); wb_height = wbPos(4);
            set(wb, 'Position', [wb_left wb_bot wb_width wb_height]);

            Sz_CA_Bones = size(CA_Bones);
            if Sz_CA_Bones(2) >= 4; CA_Bones(:,4) = {[]}; end
            for cROIidx = 1 : NumROIs

                set(CA_ROIc{cROIidx,6}, 'Color', [0 0 1])
                set(CA_ROIc{cROIidx,2}, 'FaceColor', [0 0 0])
                GreenLine = CA_ROIc{cROIidx,5};
                if ~isempty(GreenLine)
                    GreenLineDistance = CA_ROIc{cROIidx,7};
                    ConnCoo = zeros(1,2);
                    ConnCoo(1) = GreenLine(end,1);
                    ConnCoo(2) = GreenLine(end,2);
%CG: coordinates of GreenLine need to be reversed
%CG: (10Feb21) why do the coordinates of the greenlines need to be
%reversed? This causes problems with the sample data where the field of
%view is not square. 
%CG: gwConnCoo is the coordinate location for the connection point between 
%the current green line and the connecting white line. Next, find the white
%line that is the connecting line in CA_Bones. 
                    Top = ConnCoo;
%                     disp(strcat('ROI name = "', num2str(cROIidx), '"; ConnCoo = "[', num2str(ConnCoo(1)),',',num2str(ConnCoo(2)),']"',...
%                         '; Sz_fmedianRGB = "[', num2str(Sz_fmedianRGB(1)), ',', num2str(Sz_fmedianRGB(2)),']"'))
                 
                    ConnIdx = sub2ind([Sz_fmedianRGB(1),Sz_fmedianRGB(2)], ConnCoo(1), ConnCoo(2));
                    GrandTotalDistance = GreenLineDistance;
                    patchidx = []; BoneCollector = []; convidx = [];
                    while isempty(patchidx)

                        if ~isempty(convidx)
                            set(CA_Bones{convidx,1}, 'Color', [1 1 1])
                            pause(0.25)
                        end
                        cBoneidx = 0;
                        while cBoneidx < Sz_CA_Bones(1)
%CG: test each line in the CA_Bones list to find the contact point between 
%green line, or previous white line.
                            cBoneidx = cBoneidx + 1;
                            if ~isempty(find(BoneCollector == cBoneidx,1))
                                cBoneidx = cBoneidx + 1;
                            end

                            cBoneCoos = CA_Bones{cBoneidx, 3};
%CG: (10Feb21) after 'unreversing' assignment for the coordinates of the
%greenlines above, "xidx = cBoneCoos(:,1) == ConnCoo(1,1);" was changed to:
%"xidx = cBoneCoos(:,1) == ConnCoo(1,2);". Same applies for the next 3
%lines of code. 
                            xidx = cBoneCoos(:,1) == ConnCoo(1,2);
                            yidx = cBoneCoos(:,2) == ConnCoo(1,1);

                            xidx_alt = cBoneCoos(:,2) == ConnCoo(1,2);
                            yidx_alt = cBoneCoos(:,1) == ConnCoo(1,1);

                            sumidx = xidx + yidx;

                            if ~isempty(find(sumidx == 2,1)) 
                                convidx = cBoneidx;
                                cBoneidx = Sz_CA_Bones(1);
                                bingoindex = find(sumidx == 2);
                            end

                        end
                        if ~isempty(convidx)

                            set(CA_Bones{convidx,1}, 'Color', [0 0 1])
                            pause(0.1)

                            dx = diff(cBoneCoos(1:bingoindex,2));
                            dy = diff(cBoneCoos(1:bingoindex,1));

                            hyp = hypot(dx, dy); WhiteLineDistance = sum(hyp);
                            GrandTotalDistance = GrandTotalDistance + WhiteLineDistance;
%CG: test if the current white line hits a patch.
                            bumend = cBoneCoos(1,:); topend = cBoneCoos(bingoindex,:);

                            Sz_CA_StartPoints = size(CA_StartPoints);
                            cPatch = 0;

                            while cPatch < Sz_CA_StartPoints(1)
                                cPatch = cPatch + 1; 
                                PatchX = get(CA_StartPoints{cPatch}, 'XData');
                                PatchY = get(CA_StartPoints{cPatch}, 'YData');

                                PatchX_range = PatchX-10:1:PatchX+10;
                                PatchY_range = PatchY-10:1:PatchY+10;
%CG: PatchX and PatchY will only ever be single coordinates. The patch
%itself covers more than one point, so a line could be at the patch but be 
%just missing the centre point. The amount of leeway a line gets is
%hard-coded here into PatchX_range and PatchY_range. 
                                if (~isempty(find(PatchX_range == topend(1,1),1)) && ~isempty(find(PatchY_range == topend(1,2),1))) ||...
                                     (~isempty(find(PatchX_range == bumend(1,1),1)) && ~isempty(find(PatchY_range == bumend(1,2),1)))  
                                    patchidx = cPatch;
                                    set(CA_Bones{convidx,1}, 'Color', [1 1 1])
                                    cPatch = Sz_CA_StartPoints(1);
                                end
                            end

%CG: if patchidx remains empty, the current white line does
%not connect with a start point and must connect with another white line,
%which is closer to a start point. The 'bumend' would be the starting point
%for any white line. This will be the point where the current white line
%connects to the start point or the next white line in the chain closer to
%the start point. 
                            if isempty(patchidx)
                                ConnCoo = bumend;
                            end

%CG: it would be most convenient to display ROI distances based on which
%connecting white lines are selected by the user. It is not convenient to
%select ROIs individually. (it is the first pass through the isempty(patchidx)
%while loop when isempty(BoneCollector). If isempty(BoneCollector) is not
%part of the condition for entering a value into the 4th column of
%CA_Bones, then each ROI index will be assigned to all rows in CA_Bones for
%which the corresponding white line is directly or indirectly connected to
%the current ROIs green line. The cROIidx should be entered into the row
%only for the white line it is directly connected to (throw the ROI's green
%line).
                            if Sz_CA_Bones(2) > 3 && isempty(BoneCollector)
                                if isempty(CA_Bones{convidx,4})
                                    CA_Bones{convidx,4} = cROIidx;
                                else
                                    ROIindices = CA_Bones{convidx,4};
                                    ROIindices = cat(2,ROIindices,cROIidx);
                                    CA_Bones{convidx,4} = ROIindices;
                                end
                            elseif Sz_CA_Bones(2) <= 3 && isempty(BoneCollector)
                                CA_Bones{convidx,4} = cROIidx;
                            end
                            Sz_CA_Bones = size(CA_Bones);
%CG: BoneCollector stores location info for the white lines that connect 
%the current ROI to a start point. 
                            if isempty(BoneCollector)
                                BoneCollector = convidx;
                            else
                                BoneCollector(end+1) = convidx;
                            end
                        elseif isempty(convidx)
                            patchidx = NaN;
                        end

                    end
                    if ~isempty(convidx)
                        CA_Distances{CA_ROIc{cROIidx,3}} = GrandTotalDistance/Resolution;
                        CA_ROIc{cROIidx, 8} = BoneCollector; cIntBone = 0;
                        iDistances = zeros(1,numel(BoneCollector));
                        while cIntBone < numel(BoneCollector)
%CG: this loop collects the distances for each point of intersection
%between white lines. 
                            cIntBone = cIntBone + 1;
                            cBone = CA_Bones{BoneCollector(cIntBone),3};

                            TopX_logic = cBone(:,1) == Top(1);
                            TopY_logic = cBone(:,2) == Top(2);

                            TopX_logic_alt = cBone(:,2) == Top(1);
                            TopY_logic_alt = cBone(:,1) == Top(2);

                            sumtop = TopX_logic + TopY_logic;

                            sumtop_alt = TopX_logic_alt + TopY_logic_alt;

                            if ~isempty(find(sumtop == 2,1)) && isempty(find(sumtop_alt == 2,1)) 
                               bingoidx = find(sumtop == 2);
                            elseif isempty(find(sumtop == 2,1)) && ~isempty(find(sumtop_alt == 2,1))
                               bingoidx = find(sumtop_alt == 2);
                            elseif ~isempty(find(sumtop == 2,1)) && ~isempty(find(sumtop_alt == 2,1))
                               bingoidx = find(sumtop == 2);
                            end

                            dx = diff(cBone(1:bingoidx,1));
                            dy = diff(cBone(1:bingoidx,2));

                            Top = cBone(1,:); hyp = hypot(dx, dy);
                            iDistances(cIntBone) = sum(hyp)/Resolution;
                        end
                        iD_update = iDistances;
                        for cidx = 1 : numel(iDistances)
                            iD_update(cidx) = sum(iDistances(cidx:end));
                        end
                        CA_ROIc{cROIidx, 9} = iD_update;
                        CA_ROIc{cROIidx, 10} = GrandTotalDistance/Resolution;
                    end
                end
                if ~isempty(GreenLine)
                    if ~isempty(convidx)
                        currentROI = CA_ROIc{cROIidx,2};
                        if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                            set(CA_ROIc{cROIidx,2}, 'FaceColor', [1 1 1])
                        elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                            set(CA_ROIc{cROIidx,2}, 'FaceColor', anotherColor)
                        else
                            set(CA_ROIc{cROIidx,2}, 'FaceColor', [1,0,0])
                        end
                        set(CA_ROIc{cROIidx,6}, 'Color', [0 1 0])
                    end
                end
                waitbar(cROIidx/NumROIs, wb) 
            end
            pause(1)
            save(CellBonesFile, 'CA_Bones', 'CA_ROIc', 'CA_Distances', '-append')

            close(wb)
    
%CG: verify that distances are correct. 
        end
    end
    if abortswitch == 0
        if strcmp(VD, 'On')
            disp('verify distances'); k = [];
            while isempty(k)
                k = waitforbuttonpress;
                if k == 1 
                    KeyPressRes = get(fh,'CurrentCharacter');
                    if strcmp(KeyPressRes, 'a')
                    else
                        k = [];
                    end
                elseif k == 0
                    k = []; mouseclick = get(fh, 'SelectionType');
                    if strcmp(mouseclick, 'normal')
                        cObj = gco; ClassDef = class(cObj);
                        PointerStatus = get(fh, 'Pointer');
                        if strcmp(ClassDef, 'matlab.graphics.primitive.Patch') && strcmp(PointerStatus, 'arrow')

                            xdata = get(cObj, 'XData'); ydata = get(cObj, 'YData');
                            tgtidx = sub2ind([xdim ydim], ydata, xdata);
                            Sz_CA_ROIc = size(CA_ROIc); cROI = 0;
                            while cROI < Sz_CA_ROIc(1)
                                cROI = cROI + 1;
                                if ~isempty(CA_ROIc{cROI, 6})
                                    iROI = CA_PatchH_ROIs{cROI, 1};
                                    xdata_iROI = get(iROI, 'XData');
                                    ydata_iROI = get(iROI, 'YData');
                                    idx_iROI = sub2ind([xdim ydim], ydata_iROI, xdata_iROI);
                                    uCat = unique(cat(2, idx_iROI, tgtidx));

                                    if numel(uCat) <= numel(tgtidx)
                                        set(CA_PatchH_ROIs{cROI, 1}, 'FaceColor', [0 0 1])
                                        pause(0.25)
                                        set(CA_ROIc{cROI, 6}, 'Color', [0 0 1])

                                        BoneCollection = CA_ROIc{cROI, 8};
                                        Sz_BoneCollection = size(BoneCollection);
                                        NumBones = Sz_BoneCollection(1)*Sz_BoneCollection(2);

                                        for cBone = 1 : NumBones
                                            pause(0.25); set(CA_Bones{BoneCollection(cBone),1}, 'Color', [0 0 1])
                                        end

                                        distance = CA_Distances{CA_ROIc{cROI, 3}};
                                        disp(strcat('ROI index: "', num2str(cROI), '"'))
                                        disp(strcat('total distance to connected start point: "', num2str(distance), '"um'))

                                        cROI = Sz_CA_ROIc(1);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        datadir_char = datadir; datadir = {}; datadir{1}=datadir_char;
        SaveDend = 0; 
        disbran = struct;
        disbran.datadir = datadir; disbran.SaveDend = SaveDend; disbran.Resize = Resize;
        disbran.fmedianRGB = fmedianRGB;
        set(fh, 'Visible', 'off')
        skeleton.CA_StartPoints = CA_StartPoints; skeleton.CA_Bones = CA_Bones; 
        skeleton.CA_ROIc = CA_ROIc; skeleton.CA_Distances = CA_Distances;

        [CA_id,CA_Distances,subsidROIList] = PickROIs(disbran, skeleton, whiteROI, anotherROI, anotherColor);
        skeleton.CA_id = CA_id; skeleton.CA_Distances = CA_Distances; skeleton.subsidROIList = subsidROIList;
        save(CellBonesFile, 'CA_id', 'CA_Distances', 'subsidROIList','-append')
        
%CG: inverting list of children for current figure puts the ROIs on top of
%the green lines.
        figchi = get(gca, 'Children');
        inv_figchi = flipud(figchi); inv_figchi(1) = [];
        inv_figchi = cat(1, inv_figchi, figchi(end));
        set(gca, 'Children', inv_figchi)

        if ~isempty(CA_id)
            [CA_Bones] = DiscardDeselectedROIs(skeleton, whiteROI, anotherROI, anotherColor);
            skeleton.CA_Bones = CA_Bones;
            save(CellBonesFile, 'CA_Bones','-append')

            title('Done!', 'FontSize',16)
        else 
            title('Why not select some ROIs? Currently, you have selected no data for further analysis', 'FontSize', 16)
        end
    end
elseif CBfilefound == 1 && Overwrite == 0
    errordlg('CellBones call aborted to avoid unintended data loss at this stage. Delete CellBonesCAs file or use input parameters to prevent abort')
end
if abortswitch == 1
    errordlg('Closing the gui window aborts analysis. Doing this can result in data loss!')
end    
