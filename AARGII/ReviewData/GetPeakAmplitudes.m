function GetPeakAmplitudes(varargin)

p = inputParser;

addParameter(p,'suffixStrs',{'a_baseline', 'b_TTAP2'})
addParameter(p,'IntensityDataType','raw') %CG: use 'raw' or 'normalized'
addParameter(p,'ShowTrace',[])
%CG: controls which traces will be displayed. Values correspond to the
%level of activity. if ShowTrace = [1 2 3], the traces for the 1st, 2nd and
%3rd most active ROIs will be displayed.
addParameter(p,'SelectTrace',[])
%CG: similar to ShowTrace, except SelectTrace will show the trace for which
%the CA_ROIs index value is given. 
%CG: for ROIs in supplementary figure 5: [1519 1530 4487 6554 8603]
addParameter(p,'ConditionSpec',[1])
%CG: Specify condition for which ShowTrace or SelectTrace should apply.
addParameter(p,'TraceYLim',[600 1800])
%CG: set the lower and upper limits for trace plot(s). 
addParameter(p,'OrderedDimLimit',[])
%CG: useful for code development/debugging
addParameter(p,'OrderOfROIs','MostActive-LeastActive')
%CG: control how ROIs appear in heatmap: (top of graph)-(bottom of graph). 
%Available options: 'random', 'MostActive-LeastActive' or
%'LeastActive-MostActive'.
addParameter(p,'haf','On')
%CG: if haf (HighlightAmplitudeFrames) is 'On' then a second heatmap is
%created where only frames between baseline and peak amplitude times will
%be shown - all other values will be set to zero. This is useful for
%checking that the automated baseline and peak detection algorithm is
%working.
addParameter(p,'minval',0)
%CG: If there are ROIs established after baseline, these will have invalid
%normalized fluorescence intensity values of 0. Once these ROIs have been
%established they will have valid normalized intensity values (i.e. values
%greater than zero), so the minimum value will never reach zero. In such
%cases, the simplest solution is to force the minimum normalized intensity
%value to be zero acorss all conditions.
addParameter(p,'maxtick',30)
%CG: set maxtick to empty in order to find the maximum intensity value. The
%script will report the highest value in the command line. Note: if the
%highest intensity value for the current cell does not exceed the value for
%a preceeding cell, the highest value for the preceeding cell will be
%reported.
addParameter(p,'DisplayHeatmaps','Off')
%CG: if 'On' then GetPeakAmplitudes calls GetHeatMaps.m and heatmap.m to
%plot the first heatmap showing changes in fluorescence intensity
%(normalized or raw data according to the definition of IntensityDataType)
%for each ROI included in the data analysis. 
addParameter(p,'xScale','minutes') %CG: 'minutes' or default will be in seconds;
addParameter(p,'xSpaceCtrl',0.1) %CG: control spacing along the x-axis
addParameter(p,'vZoom','Off')
%CG: turn vertical zoom on the heatmaps 'On' or 'Off'
addParameter(p,'vZoomLim',[1 15])
%CG: set limits of the y-axis. The range (vZoomLim) indicates the number of ROIs to be
%displayed.
addParameter(p,'hZoom','Off')
%CG: turn horizontal zoom on the heatmaps 'On' or 'Off'
addParameter(p,'hZoomLim',[])
%CG: Set limits of the x-axis. go to the xScale variable to check units.
addParameter(p,'FrameRate',20)
%CG: only applicable for the heatmap.
addParameter(p,'BinWidth',[]) 
%CG: BinWidth empty means BinWidth = 1/FrameRate; i.e. there will be no
%binning.
addParameter(p,'ColorbarLabel','arbitary units (a.u.)')
%CG: color bar label for first heatmap.
addParameter(p,'cmap','jet')
%CG: colormap for heatmap. 'jet' recommended to get full visible spectrum.
addParameter(p,'OnlyInts','On')
%CG: Set to 'On' to remove non-integers from colorbar display. Any other
%string will include non-integer values. 0 also included. Must be string.
addParameter(p,'Dirs',{})
parse(p,varargin{:}); param = p.Results; 

suffixStrs = param.suffixStrs; IntensityDataType = param.IntensityDataType;
ShowTrace = param.ShowTrace;  SelectTrace = param.SelectTrace; ConditionSpec = param.ConditionSpec;
TraceYLim = param.TraceYLim; OrderedDimLimit = param.OrderedDimLimit; 
OrderOfROIs = param.OrderOfROIs; haf = param.haf;minval = param.minval; 
maxtick = param.maxtick; DisplayHeatmaps = param.DisplayHeatmaps;
xScale = param.xScale; xSpaceCtrl = param.xSpaceCtrl; vZoom = param.vZoom;
vZoomLim = param.vZoomLim; hZoom = param.hZoom; hZoomLim = param.hZoomLim;
FrameRate = param.FrameRate; BinWidth = param.BinWidth; ColorbarLabel = param.ColorbarLabel;
cmap = param.cmap; OnlyInts = param.OnlyInts; Dirs = param.Dirs;

markerSize = 34; %CG: standard marker size = 18. 
lineWidth = 0.5; %CG: default line width = 0.5;

if isempty(BinWidth); BinWidth = 1/FrameRate; end
Sz_suffixStrs = size(suffixStrs); ConditionNames = {'baseline', '+conotoxin'};

if isempty(maxtick); DispMaxTick = 'On';
else maxtick = maxtick + 1; DispMaxTick = 'Off'; end

if isempty(maxtick); DisplayHeatmaps = 'Off'; end

FpB = BinWidth*FrameRate;
BinIntStore = {}; xAxLabel_CA = {}; xlab_CA = {}; CA_ROIs_store = {};

%Outputs
%CA_Measure contains measurements for baseline and peak amplitude of
%fluorescence intensity for every ROI that was selected for further
%analysis in CellBones.
%Details for CA_Measure contents:
%CA_Measure{:,1}: row number.
%CA_Measure{:,2}: Index value locating the current ROI in CA_ROIs etc. 
%CA_Measure{:,3}: baseline values for each trace segment* within current roi. 
%CA_Measure{:,4}: peak values for each trace segment within current roi.
%CA_Measure{:,5}: exact baseline values of each detection within current trace.
%CA_Measure{:,6}: exact peak values of each detection within current trace.
%CA_Measure{:,7}: average baseline values of each detection within current trace.
%CA_Measure{:,8}: average peak values of each detection within current trace.
%CA_Measure{:,9}: peak amplitude values (peak-baseline) of each detection within current trace.
%CA_Measure{:,10}: number of elements contained within each trace segment.
%CA_Measure{:,11}: y-values for each trace segment.
%CA_Measure{:,12}: x-values for each trace segment. 


%* a 'trace' is the complete fluorescence signal (cfs) for a single roi collected
%in one condition. A 'trace segment' is a small section of the cfs where
%one event was detected. If there were 10 detections within the current roi
%or cfs, then there will be 10 trace segments for the current roi or cfs.
scrsz=get(0,'ScreenSize');
%Dirs = uipickfiles('Prompt', 'Select data folders');
% Dirs = {}; Dirs{1} ='/Users/cjgilbride/Documents/DecJan2016_17/MATLAB/AARG_pub_datatest/CG0706173';

if iscell(Dirs)
    
    if strcmp(xScale, 'minutes'); xlabel_str = 'Time (minutes)';
    elseif ~strcmp(xScale, 'minutes'); xlabel_str = 'Time (s)'; end
    
    Slashes = strfind(Dirs{1}, '/'); LastSlash = max(Slashes);
    FirstDir = Dirs{1}; %ParentDir = FirstDir(1 : LastSlash-1);
    
    NumOfDirs = size(Dirs,2); NumOfConditions = Sz_suffixStrs(1)*Sz_suffixStrs(2);
    for cDir = 1 : NumOfDirs
        
        datadir = Dirs{cDir}; OrderedDimLimit = []; 
        ca = 0; dataStore = {}; cd(datadir); [~, CellName, ~] = fileparts(datadir);
        try 
            if strcmp(IntensityDataType,'raw')
                load(strcat(CellName, '_hmvarcol_raw.mat'), 'dataStore', 'intensityStore', 'BinIntStore', 'ROIidxstore', 'es')
            elseif strcmp(IntensityDataType,'normalized')
                load(strcat(CellName, '_hmvarcol_norm.mat'), 'dataStore', 'intensityStore', 'BinIntStore', 'ROIidxstore', 'es')
            end
            CellName = es.CellName; NumROIsRIS = es.NumROIsRIS;
            r_co = es.r_co; c_co = es.c_co; z_co = es.z_co;
            savedIntensityDataType = es.IntensityDataType; edges = es.edges;
            if ~strcmp(IntensityDataType, savedIntensityDataType); BreakTry; end
            bones = load(strcat(CellName, '_CellBonesCAs.mat'));
            breaktry
            CA_ROIc = bones.CA_ROIc; ROIcIdces = cell2mat(CA_ROIc(:,3));
        catch
            bones = load(strcat(CellName, '_CellBonesCAs.mat'));
            lists = load(strcat(CellName,'_FileNameLists.mat'));
            ExpNameList = lists.ExpNameList; NumStrFilesList = lists.NumStrFilesList;
            
            NumEventsPerROI = []; ROIidxstore = []; tor = 0;
            CA_ROIc = bones.CA_ROIc; ROIcIdces = cell2mat(CA_ROIc(:,3));
            for cConditionIdx = 1 : NumOfConditions
                cdtidx = find(NumStrFilesList == cConditionIdx);
                partdfname = ExpNameList{cdtidx(end)}; tifidx = strfind(partdfname, '.tif');
                tgtfolder = strcat(partdfname(1:tifidx-1), '_ThresholdData');
                TDDir = strcat(datadir, '/', tgtfolder); 
                if iscell(TDDir)
                    cd(TDDir{1}); fileorder = load(strcat(tgtfolder{1}, '_FileOrder.mat'));
                else cd(TDDir); fileorder = load(strcat(tgtfolder, '_FileOrder.mat'));
                end 

                TotalCCSList = fileorder.TotalCCSList; TotalFrames = TotalCCSList(end);
%CG: get the total number of frames in each condition.                
                if cConditionIdx == 1
                    
                    try
                        cd(strcat(datadir, '/', CellName, '_fIntensity'))
                        intensities = load(strcat(CellName, suffixStrs{cConditionIdx}, '_IntensityFiles.mat'));
                        lastIntensities = load(strcat(CellName, suffixStrs{end}, '_IntensityFiles.mat'));
                        if strcmp(IntensityDataType, 'raw')
                            CA_Intensity = intensities.CA_IntensityRaw;
                            lastCA_Intensity = lastIntensities.CA_IntensityRaw;
                        elseif strcmp(IntensityDataType, 'normalized')
                            CA_Intensity = intensities.CA_Intensity170;
                            lastCA_Intensity = intensities.CA_Intensity170;
                        end
                    catch
                        errordlg('ERROR: call GetIntensities before calling GetPeakAmplitudes (or perhaps a Keyword is mis-spelled?)')
                    end
                    cd(strcat(datadir, '/', CellName, suffixStrs{end}, '_ThresholdData'))
                    lastdata = load(strcat('CAs_Core_', CellName, suffixStrs{end}, '.mat'));
                    cd(strcat(datadir, '/', CellName, suffixStrs{1, cConditionIdx}, '_ThresholdData'))
                elseif cConditionIdx > 1
                    cd(strcat(datadir, '/', CellName,  '_fIntensity'))
                    intensities = load(strcat(CellName, suffixStrs{cConditionIdx}, '_IntensityFiles.mat'));
                    if strcmp(IntensityDataType, 'raw')
                        CA_Intensity = intensities.CA_IntensityRaw;
                    elseif strcmp(IntensityDataType, 'normalized')
                        CA_Intensity = intensities.CA_Intensity170;
                    end
                end
                cd(strcat(datadir, '/', CellName, suffixStrs{cConditionIdx}, '_ThresholdData'))
                data = load(strcat('CAs_Core_', CellName, suffixStrs{cConditionIdx}, '.mat'));
                    
                CA_ROIs = data.CA_ROIs; CA_NumEvents = data.CA_NumEvents; CA_ROIs_store{cConditionIdx} = CA_ROIs;
                Sz_CA_NumEvents = size(CA_NumEvents); Sz_CA_Intensity = size(CA_Intensity); CA_ROIsEnd = lastdata.CA_ROIs;
                dataStore{cConditionIdx} = data; intensityStore{cConditionIdx} = intensities;
                if tor == 0
                    if cConditionIdx == 1 && isempty(OrderedDimLimit) 
                        
                        Sz_lastCA_Intensity = size(lastCA_Intensity);
                        if numel(Sz_lastCA_Intensity) == 2
                            DimLimit = Sz_lastCA_Intensity(1)*Sz_lastCA_Intensity(2);
%                             DimLimit = Sz_CA_Intensity(1)*Sz_CA_Intensity(2);
                        elseif numel(Sz_lastCA_Intensity) == 3
                            DimLimit = Sz_lastCA_Intensity(1)*Sz_lastCA_Intensity(2)*Sz_lastCA_Intensity(3);
%                             DimLimit = Sz_CA_Intensity(1)*Sz_CA_Intensity(2)*Sz_CA_Intensity(3);
                        end
                        OrderedDimLimit = randperm(DimLimit);
%CG: eventhough the ROIs might eventually be ordered with the most active
%ROIs placed higher up in the list, an underlying pattern will emerge that
%has no biological meaning. It just reflects single ROIs having only one
%event and how established ROIs are organized in the cell array. 
                    end
%                     DimLimit = Sz_CA_Intensity(1)*Sz_CA_Intensity(2)*Sz_CA_Intensity(3);
                    for gem = 1 : DimLimit

                        [rO, cO, zO] = ind2sub(Sz_lastCA_Intensity, OrderedDimLimit(gem));
%                         [rO, cO, zO] = ind2sub(Sz_CA_Intensity, OrderedDimLimit(gem));
                        if numel(Sz_CA_Intensity) == 2; ThirdSzLim = 1; else ThirdSzLim = Sz_CA_Intensity(3); end
                        if rO <= Sz_CA_Intensity(1) && cO <= Sz_CA_Intensity(2) && zO <= ThirdSzLim
                            cInt = CA_Intensity{rO, cO, zO};
                            cROIlast = CA_ROIsEnd{rO, cO, zO};

                            if rO <= Sz_CA_NumEvents(1) && cO <= Sz_CA_NumEvents(2) && zO <= ThirdSzLim  
                                cNum = CA_NumEvents{rO, cO, zO};
                            else cNum = [];
                            end
                            ROIcrow = ROIcIdces == OrderedDimLimit(gem);
                            if ~isempty(cInt) && ~isempty(cNum) && ~isempty(CA_ROIc{ROIcrow,6})
                                if cConditionIdx == 1
                                    if isempty(NumEventsPerROI)
                                        NumEventsPerROI = CA_NumEvents{rO, cO, zO};
                                        ROIidxstore = OrderedDimLimit(gem);
                                    elseif ~isempty(NumEventsPerROI)
                                        NumEventsPerROI(end+1, 1) =  CA_NumEvents{rO, cO, zO};
                                        ROIidxstore(end+1, 1) = OrderedDimLimit(gem);
                                    end
                                else
                                    MatchDex = find(ROIidxstore(:,1) == OrderedDimLimit(gem));
                                    NumEventsPerROI(MatchDex, cConditionIdx) =  CA_NumEvents{rO, cO, zO};
                                    ROIidxstore(MatchDex, cConditionIdx) = OrderedDimLimit(gem);
                                end
                            elseif ~isempty(cInt) && isempty(cNum) && ~isempty(CA_ROIc{ROIcrow,6})
%CG: cNum can be empty if the ROI was established by activity in another
%condition and there was no overlapping activity in this same ROI in the
%current condition. 
                                if cConditionIdx == 1
                                    if isempty(NumEventsPerROI)
                                        NumEventsPerROI = 0;
                                        ROIidxstore = OrderedDimLimit(gem);
                                    elseif ~isempty(NumEventsPerROI)
                                        NumEventsPerROI(end+1, 1) =  0;
                                        ROIidxstore(end+1, 1) = OrderedDimLimit(gem);
                                    end
                                else
                                    MatchDex = find(ROIidxstore(:,1) == OrderedDimLimit(gem));
                                    NumEventsPerROI(MatchDex, cConditionIdx) =  0;
                                    ROIidxstore(MatchDex, cConditionIdx) = OrderedDimLimit(gem);
                                end
                            elseif isempty(cInt) && ~isempty(cNum)% && ~strcmp(cNum, 'subsid') && ~strcmp(cROIlast, 'subsid')
                                disp(num2str([rO, cO, zO]))
                                errordlg('CA_Intensity cell should not be empty while CA_NumEvents is non-empty')
                            
%CG: if strcmp(cROIlast, 'subsid') then an event in a non-baseline
%condition has found the largest event of a pre-established ROI to be
%overlapped entirely with an event occurring in the current condition. 
                            elseif isempty(cInt) && ~isempty(cNum) && ~strcmp(cNum, 'subsid') && strcmp(cROIlast, 'subsid')
                                ca = ca + 1;
                            end
                        end
                    end
                    if cConditionIdx == 1 && ~strcmp(OrderOfROIs, 'random')
                        if strcmp(OrderOfROIs, 'MostActive-LeastActive')
                            [~, sortingIdces] = sort(NumEventsPerROI, 'descend');
                        elseif strcmp(OrderOfROIs, 'LeastActive-MostActive')
                            [~, sortingIdces] = sort(NumEventsPerROI, 'ascend');
                        end
                        ROIidxstore = ROIidxstore(sortingIdces); NumROIsRIS = size(ROIidxstore,1);
                        [r_co, c_co, z_co] = ind2sub(Sz_CA_Intensity, ROIidxstore);
                    end

                    TotalBrun = numel(ROIidxstore(:,1));
                    edges = FpB : FpB : TotalFrames; edges = cat(2, 1, edges); edges = unique(edges);
                    if ~strcmp(hZoom, 'On')
                        hZoomEndPnt = numel(edges)-1; BinIntData = zeros(TotalBrun, numel(edges)-1);
                        startedge = 1; endedge = hZoomEndPnt;
                    elseif strcmp(hZoom, 'On')
                        dhZoomLim = hZoomLim(2) - hZoomLim(1);
                        if strcmp(xScale, 'minutes')
                            dhZoomLim = dhZoomLim*60; hZoomEndPnt = round((dhZoomLim*FrameRate));
                        else
                            hZoomEndPnt = round((dhZoomLim*FrameRate));
                        end
                        BinIntData = zeros(TotalBrun, hZoomEndPnt);
                        startedge = hZoomLim(1)*FrameRate; endedge = hZoomLim(2)*FrameRate;
                    end

                    for brunette = 1 : TotalBrun

                        [r1, c1, z1] = ind2sub(Sz_CA_Intensity, ROIidxstore(brunette));

                        if numel(Sz_CA_Intensity) == 2; ThirdSzLim = 1;
                        elseif numel(Sz_CA_Intensity) ~= 2; ThirdSzLim = Sz_CA_Intensity(3);
                        end

                        if r1 <= Sz_CA_Intensity(1) && c1 <= Sz_CA_Intensity(2) && z1 <= ThirdSzLim
                            CurrentInts = CA_Intensity{r1, c1, z1};
                        else CurrentInts = [];
                        end

                        if ~isempty(CurrentInts)
                            cter = 0;
                            for cElement = startedge : endedge
                                cpoint = edges(cElement);
                                if cpoint+1 > numel(CurrentInts); cpoint = cpoint-1; end
                                cter = cter + 1;
                                if cpoint >= numel(CurrentInts); BinIntData(brunette, cter) = 0;
%CG: BinIntData is larger by two frames than CurrentInts. Need to fix this later. Won't
%cause serious issues a this stage. 
                                else BinIntData(brunette, cter) = mean(CurrentInts(cpoint:cpoint+1));
                                end
                            end
                        end
                    end
                    if strcmp(vZoom, 'On'); BinIntData = BinIntData(vZoomLim(1):vZoomLim(2), :);
                        BinIntStore{cConditionIdx} = BinIntData;
                    else BinIntStore{cConditionIdx} = BinIntData;
                    end
                end
            end
            cd(datadir)
            
            es = struct;
            es.CellName = CellName; es.NumROIsRIS = NumROIsRIS; 
            es.r_co = r_co; es.c_co = c_co; es.z_co = z_co; es.IntensityDataType = IntensityDataType;
            es.edges = edges;
            if strcmp(IntensityDataType,'raw')
                save(strcat(CellName, '_hmvarcol_raw.mat'), 'dataStore', 'intensityStore', 'BinIntStore', 'ROIidxstore', 'es')
            elseif strcmp(IntensityDataType,'normalized')
                save(strcat(CellName, '_hmvarcol_norm.mat'), 'dataStore', 'intensityStore', 'BinIntStore', 'ROIidxstore', 'es')
            end
        end

        ROIidxstore = ROIidxstore(:,1); BinIntStore_amps = {};
%CG: both columns of ROIidxstore were always identical. ROIidxstore needs
%only to be a vector array.  
      
        h1 = []; CA_Measure = {};
%CG: CA_Measure stores data collected from traces (e.g. peak amplitude). 
        for cConditionIdx = 1 : NumOfConditions
            data = dataStore{cConditionIdx}; intensities = intensityStore{cConditionIdx};
            CA_Intensity = intensities.CA_IntensityRaw; CA_Indices = data.CA_Indices;
            CA_ROIs = data.CA_ROIs;
            BinIntData = BinIntStore{cConditionIdx};
            bidlog = logical(BinIntData);
            bidlog(:) = 0;

            ss = size(CA_Indices); 
            if numel(ss) == 2; CAmax = ss(1)*ss(2); elseif numel(ss) == 3; CAmax = ss(1)*ss(2)*ss(3); end
            
            CA_Mcount = 0;
            for BID_row = 1 : NumROIsRIS

                CA_ROIc_row_idx = find(ROIcIdces == ROIidxstore(BID_row));
%CG: ROIindexstore provides a useful list of the ROIs displayed in the 
%heatmap, but ROIs in CA_ROIc are not arranged in the same order. Further, if no green
%line handle exists for the current ROI in CA_ROIc{CA_ROIc_row_idx,6} then
%the ROI was not selected for further analysis. 
                if CAmax >= ROIidxstore(BID_row)
                    peak = CA_Indices{r_co(BID_row),c_co(BID_row),z_co(BID_row)};
                end
                ROI_active = 0; 
                if ~isempty(CA_ROIc_row_idx)
                    if ~isempty(CA_ROIc{CA_ROIc_row_idx,6}); ROI_active = 1; end
                end
                if ROI_active == 1 && CAmax >= ROIidxstore(BID_row) && ~isempty(peak)
%CG: if isempty(peak) && cCondition > 1 this should mean that the current
%ROI was established in condition 1 and hosted no activity during the
%current condition.
                    CA_Mcount = CA_Mcount + 1; 
                    CA_Measure{CA_Mcount,1} = BID_row; CA_Measure{CA_Mcount,2} = ROIidxstore(BID_row);
                    inttrace = CA_Intensity{r_co(BID_row),c_co(BID_row),z_co(BID_row)};
                    
                    AmpVals=[];

                    if ~isempty(CA_ROIs{r_co(BID_row),c_co(BID_row),z_co(BID_row)})

                        freqcount = 0; base = []; avBaseVals = []; avPeakVals = [];
                        BaseIdces = []; PeakIdces = []; peakframes = peak(1,:);
                        CA_ew = {}; CA_rt = {}; valrt = []; ewnum = []; pre_pp=[]; 
                        SegBaseIdces = []; SegPeakIdces = [];
                        for cpp = 1 : numel(peakframes)
                            if strcmp(IntensityDataType, 'raw')
                                intderiv = diff(inttrace(1:peakframes(cpp)));
                                possbases = find(intderiv <= 0);  d_possbases = diff(possbases);
                                d_possbases_idx = find(d_possbases == 1); 
                                if ~isempty(d_possbases_idx)
                                    possbases = possbases(d_possbases_idx(end)+1);
                                elseif isempty(d_possbases_idx)
                                    possbases = find(intderiv <= 0);
                                end
%CG: looking for the point where the x axis is intersected works when the
%data is normalized, but not with the raw data. Find the point closest to
%the peak where the rate of change is 0 or opposite in sign to the
%derivative during the rising phase for a least 2 consecutive frames. 
                            elseif strcmp(IntensityDataType, 'normalized')
                                possbases = find(inttrace(1:peakframes(cpp)) <= 0); 
                            end
                            if ~isempty(possbases); base = [base, possbases(end)]; 
                            elseif isempty(possbases); base = [base,1]; end
                            ew_exp = 1*FrameRate;
%CG: making ew_exp too large can be problematic for curve fitting in ReviewFit  
                            Sz_BinIntData = size(BinIntData);
                            if peakframes(cpp)+ew_exp > Sz_BinIntData(2) && base(cpp)-ew_exp >= 1
                                ew  = BinIntData(BID_row, base(cpp)-ew_exp:end);
                                ewidx = base(cpp)-ew_exp:1:Sz_BinIntData(2);
                                bidlog(BID_row, base(cpp)-ew_exp:end) = 1;
                            elseif peakframes(cpp)+ew_exp <= Sz_BinIntData(2) && base(cpp)-ew_exp < 1
                                ew = BinIntData(BID_row, 1:peakframes(cpp)+ew_exp);
                                ewidx = 1:1:peakframes(cpp)+ew_exp;
                                bidlog(BID_row, 1:peakframes(cpp)+ew_exp) = 1;
                            else
                                ew = BinIntData(BID_row, base(cpp)-ew_exp:peakframes(cpp)+ew_exp);
                                ewidx = base(cpp)-ew_exp:1:peakframes(cpp)+ew_exp;
                                bidlog(BID_row, base(cpp)-ew_exp:peakframes(cpp)+ew_exp) = 1;
%CG: bidlog will be used to build a second activity map to show where
%events have been detected. 
                            end
                            peakvalidx = find(ewidx==peakframes(cpp));
                            [~, basevalidx] = sort(ew(1:peakvalidx), 'ascend'); %[~, peakvalidx] = sort(ew, 'descend');
%                             peakvalidx = find(ewidx==peakframes(cpp));                            

                            if peakvalidx > basevalidx(1) && peakvalidx ~= numel(ew)
%CG: Peak should be ahead of the baseline value and the peak should not be
%the last value in the trace (in this last case it is ambiguous whether or
%not a peak has been detected).
                                if basevalidx(1) > 1
                                    mbb = mean(ew(basevalidx(1)-1:basevalidx(1):basevalidx(1)+1)); 
                                elseif basevalidx(1) == 1
                                    mbb = mean(ew(basevalidx(1))); 
                                end
                                mpp = ew(peakvalidx);
                                seg_mbb = basevalidx(1); 
                            else
                                mbb = NaN; mpp = NaN; seg_mbb = NaN; ewmore = []; 
                            end
%CG: avoid measuring same event more than once by check if the previous top
%peak values were the same as the current values. Very unlikely that they
%would be identical if it is not the same event.
                            if ~isempty(pre_pp) && ~isnan(mpp)
                                if sum(pre_pp - mpp) == 0; mpp = NaN; end
                            end
                            pre_pp = mpp;
                            
                            ii1=find(ew(1:peakvalidx(1))<mbb);
%CG: find the index value of ew closest to the maximum peak value which
%falls below the average baseline value (bb). 
                            if ~isempty(ii1)
                                rt = (peakvalidx(1) - ii1(end))/FrameRate;
                            else
                                rt = (peakvalidx(1) - basevalidx(1))/FrameRate;
                            end

                            valrt = [valrt, [rt]]; AmpVals=[AmpVals, [mpp-mbb]];
                            avBaseVals = [avBaseVals, mbb]; avPeakVals = [avPeakVals, mpp];
                            BaseIdces = [BaseIdces, ewidx(basevalidx(1))]; PeakIdces = [PeakIdces, ewidx(peakvalidx(1))];
                            SegBaseIdces = [SegBaseIdces, seg_mbb]; SegPeakIdces = [SegPeakIdces, peakvalidx(1)];
                            if ~isnan(mpp); CA_ew{1, cpp} = ew; else CA_ew{1, cpp} = []; end; CA_rt{1,cpp} = rt;
                            ewnum = [ewnum, [numel(ew)]];
                            pre_basevalidx = ewidx(basevalidx(1));
                            pre_peakvalidx = ewidx(peakvalidx(1));
                            CA_ewidx{1, cpp} = ewidx;
                            freqcount = freqcount + 1;
                                    
                        end

                        CA_Measure{CA_Mcount, 3} = SegBaseIdces;
                        CA_Measure{CA_Mcount, 4} = SegPeakIdces;
                        CA_Measure{CA_Mcount, 5} = BaseIdces;
                        CA_Measure{CA_Mcount, 6} = PeakIdces;
                        CA_Measure{CA_Mcount, 7} = avBaseVals;
                        CA_Measure{CA_Mcount, 8} = avPeakVals;
                        CA_Measure{CA_Mcount, 9} = AmpVals;
                        CA_Measure{CA_Mcount, 10} = ewnum;
                        CA_Measure{CA_Mcount, 11} = CA_ew;
                        CA_Measure{CA_Mcount, 12} = CA_ewidx;

                    elseif isempty(CA_ROIs{r_co(BID_row),c_co(BID_row),z_co(BID_row)})
%CG: if isempty(CA_ROIs{r_co(BID_row),c_co(BID_row),z_co(BID_row)}) then
%and cCondition = 1, it probably means that the current ROI was established
%only in the second condition.
                        CA_Measure{CA_Mcount, 3} = [];
                        CA_Measure{CA_Mcount, 4} = [];
                        CA_Measure{CA_Mcount, 5} = [];
                        CA_Measure{CA_Mcount, 6} = [];
                        CA_Measure{CA_Mcount, 7} = [];
                        CA_Measure{CA_Mcount, 8} = [];
                        CA_Measure{CA_Mcount, 9} = [];
                        CA_Measure{CA_Mcount, 10} = [];
                        CA_Measure{CA_Mcount, 11} = [];
                        CA_Measure{CA_Mcount, 12} = [];
                    end
                elseif ROI_active == 0 && CAmax < ROIidxstore(BID_row)
                    CA_Mcount = CA_Mcount + 1; 
                    CA_Measure{CA_Mcount,1} = BID_row;
                    CA_Measure{CA_Mcount,2} = ROIidxstore(BID_row);
                    CA_Measure{CA_Mcount, 3} = [];
                    CA_Measure{CA_Mcount, 4} = [];
                    CA_Measure{CA_Mcount, 5} = [];
                    CA_Measure{CA_Mcount, 6} = [];
                    CA_Measure{CA_Mcount, 7} = [];
                    CA_Measure{CA_Mcount, 8} = [];
                    CA_Measure{CA_Mcount, 9} = [];
                    CA_Measure{CA_Mcount, 10} = [];
                    CA_Measure{CA_Mcount, 11} = [];
                    CA_Measure{CA_Mcount, 12} = [];
                elseif isempty(peak) && cConditionIdx > 1
                    CA_Mcount = CA_Mcount + 1; 
                    CA_Measure{CA_Mcount,1} = BID_row;
                    CA_Measure{CA_Mcount,2} = ROIidxstore(BID_row);
                    CA_Measure{CA_Mcount, 3} = [];
                    CA_Measure{CA_Mcount, 4} = [];
                    CA_Measure{CA_Mcount, 5} = [];
                    CA_Measure{CA_Mcount, 6} = [];
                    CA_Measure{CA_Mcount, 7} = [];
                    CA_Measure{CA_Mcount, 8} = [];
                    CA_Measure{CA_Mcount, 9} = [];
                    CA_Measure{CA_Mcount, 10} = [];
                    CA_Measure{CA_Mcount, 11} = [];
                    CA_Measure{CA_Mcount, 12} = [];
                end
            end

            if ~isempty(find(ConditionSpec==cConditionIdx,1))
                if ~isempty(ShowTrace)
                    for ctraceidx = 1 : numel(ShowTrace)
                        ctrace = ShowTrace(ctraceidx);
                        plotCompleteTrace(ctrace, CA_Measure, BinIntData, ROIidxstore,... 
                            lineWidth, markerSize, ShowTrace, [], suffixStrs, ctraceidx, CellName,...
                            cConditionIdx, 'ShowTrace', TraceYLim)
                    end
                    plotHistogram(h1, ctrace, ctraceidx, CA_Measure, cConditionIdx, CA_Mcount, ShowTrace)
                end
                if ~isempty(SelectTrace)
                    for cSelectTraceIdx = 1 : numel(SelectTrace)
                        ctrace = find(cell2mat(CA_Measure(:,2)) == SelectTrace(cSelectTraceIdx));
%                         ctrace = SelectTrace(ctraceidx);
                        plotCompleteTrace(ctrace, CA_Measure, BinIntData, ROIidxstore,... 
                            lineWidth, markerSize, [], SelectTrace,suffixStrs, cSelectTraceIdx, CellName,...
                            cConditionIdx, 'SelectTrace', TraceYLim)
                    end
                end
                    
            end
            
            BinIntData_binary = BinIntData;
            invbid = bidlog == 0;
            BinIntData_binary(invbid) = 0;
            BinIntData_binary(bidlog) = 1;
            BinIntStore_binary{cConditionIdx} = BinIntData_binary;
            if strcmp(IntensityDataType, 'raw')
                try
                    measures = load(strcat(CellName, suffixStrs{cConditionIdx},'_MeasuresRAW.mat'));
                    CA_Measure = measures.CA_Measure; IntensityDataType = measures.IntensityDataType;

                catch
                    save(strcat(CellName, suffixStrs{cConditionIdx},'_MeasuresRAW.mat'), 'CA_Measure', 'IntensityDataType');
                end
            elseif strcmp(IntensityDataType, 'normalized')
                
                try
                    measures = load(strcat(CellName, suffixStrs{cConditionIdx},'_MeasuresNorm.mat'));
                    CA_Measure = measures.CA_Measure; IntensityDataType = measures.IntensityDataType;
                catch
                    save(strcat(CellName, suffixStrs{cConditionIdx},'_MeasuresNorm.mat'), 'CA_Measure', 'IntensityDataType');
                end
            end
        end
%             MCV = [];
%             if strcmp(DisplayHeatmaps, 'On')
%                 hmvars = struct;
%                 hmvars.datadir = datadir;
%                 hmvars.suffixStrs = suffixStrs;
%                 hmvars.ConditionNames = ConditionNames;
%                 hmvars.OrderOfROIs = OrderOfROIs;
%                 hmvars.FpB = FpB;
%                 hmvars.cmap = cmap;
%                 hmvars.OnlyInts = OnlyInts;
%                 hmvars.BinIntStore = BinIntStore;
%                 hmvars.edges = edges;
%                 hmvars.CellName = CellName;
%                 hmvars.scrsz = scrsz;
%                 hmvars.xSpaceCtrl = xSpaceCtrl;
%                 hmvars.hZoom = hZoom;
%                 hmvars.xScale = xScale;
%                 hmvars.hZoomLim = hZoomLim;
%                 hmvars.FrameRate = FrameRate;
%                 hmvars.maxtick = maxtick;
%                 hmvars.DispMaxTick = DispMaxTick;
%                 hmvars.minval = minval;
%                 hmvars.MCV = MCV;
%                 hmvars.xlabel_str = xlabel_str;
%                 hmvars.ColorbarLabel = ColorbarLabel;
% %CG: show heatmap for CA_Intensity data. 
%                 GetHeatMaps(hmvars)
%                 
%                 grayCtrl = linspace(1,0,maxtick);
%                 grayCtrl = grayCtrl';
%                 cmap = cat(2, grayCtrl, grayCtrl, grayCtrl);
%                 hmvars.cmap = cmap;
%                 hmvars.maxtick = 1;
%                 hmvars.BinIntStore = BinIntStore_binary;
% %CG: heatmaps shows where activity was detected.                 
%                 GetHeatMaps(hmvars)
%                 
%             end
    end
end
end

function plotCompleteTrace(ctrace, CA_Measure, BinIntData, ROIidxstore,... 
    lineWidth, markerSize, ShowTrace, SelectTrace, suffixStrs, ctraceidx, CellName,...
    cConditionIdx, TraceType, TraceYLim)

    ybases = CA_Measure{ctrace,7}; ypeaks = CA_Measure{ctrace,8};
    xbases = CA_Measure{ctrace,5}; xpeaks = CA_Measure{ctrace,6};
    figure('Name', strcat('CA_ROIs index: "', num2str(ROIidxstore(CA_Measure{ctrace,1})),'"'));
    plot(BinIntData(CA_Measure{ctrace,1},:),'k', 'LineWidth', lineWidth)
    ax = gca; ax.TickDir = 'out'; 
    if ~isempty(ybases)
        hold on
        plot(xbases,ybases,'.','Color', [0 0.75 0], 'MarkerSize', markerSize)
        plot(xpeaks,ypeaks,'r.','MarkerSize', markerSize)
        hold off
        if strcmp(TraceType, 'ShowTrace')
            if ShowTrace(ctraceidx) == 1
                titlestr = '1st most active ROI for "';
            elseif ShowTrace(ctraceidx) == 2
                titlestr = '2nd most active ROI for "';
            elseif ShowTrace(ctraceidx) == 3
                titlestr = '3rd most active ROI for "';
            else
                titlestr = strcat(num2str(ShowTrace(ctraceidx)), 'th most active ROI for "');
            end
            title(strcat(titlestr, CellName, '"; Condition "', suffixStrs{cConditionIdx}, '"'))
        elseif strcmp(TraceType, 'SelectTrace')
            titlestr = strcat('Trace from ROI CA_ROIs{',num2str(SelectTrace), '}. "'); 
            title(strcat(titlestr, CellName, '"; Condition "', suffixStrs{cConditionIdx}, '"')) 
        end
        ax1 = gca;
        if ~isempty(TraceYLim); ax1.YLim = TraceYLim; end
    end
end

function plotHistogram(h1, ctrace, ctraceidx, CA_Measure, cConditionIdx, CA_Mcount, ShowTrace)

    
    figure;
    if isempty(h1)
        h1 = histogram(CA_Measure{ctrace,9});
    else
        histogram(CA_Measure{ctrace,9},h1.BinEdges)
    end
    ax1i = gca;
    ax1i.YLim = [0 25];
    ylabel('Number of events')
    xlabel('Amplitude (r.f.u.)')
    if ShowTrace(ctraceidx) == 1 
        ornind = 'st'; 
    elseif ShowTrace(ctraceidx) == 2; 
        ornind = 'nd';
    elseif ShowTrace(ctraceidx) == 3; 
        ornind = 'rd'; 
    else
        ornind = 'th';
    end

    title(strcat('"',num2str(ShowTrace(ctraceidx)),ornind,'" most active ROI in condition: "', num2str(cConditionIdx), '"')) 
    CA_ewexpan = CA_Measure{CA_Mcount,11};
    numevents = size(CA_ewexpan, 2);
end




