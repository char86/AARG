function AARG_CoreHandler(varargin)

% Outline
% AARG_CoreHandler should be called by AARGI once all necessary outputs
% from the Detection function have been saved in the matlab path.
% AARG_CoreHandler manages the flow of data to AARG_Core which carries out
% the automatic ROI assignment. 

% Known issue: AARG_CoreHandler calls AARG_Core in two independent
% for-loops. The first for-loop handles the first experiment condition
% only, while the remaining conditions are handled in the second for-loop.
% This separation is not necessary. This will not significantly impact
% performance but affects the maintainability of the code if there is no
% logical reason for having significantly more lines of code than
% necessary.

% Author: Charlie J Gilbride
% version: 1.0.0

p = inputParser;

addParameter(p,'datadirs',{})
addParameter(p,'suffixStrs',{'a_baseline', 'b_mib'})
addParameter(p,'partnamesuff','')
%CG: partitioned file string
addParameter(p,'ROI_sl', 3)   
%CG: 'ROI_sl' must be an odd number starting from 3. This sets the size of
%the ROIs, which will have the dimensions: 'ROI_sl x ROI_sl' (units: pixel/elements) 
addParameter(p,'PDE',70)
%CG: IEI is the inter-event interval. Frequency distributions are displayed by
%'DetectEvents' call, which should help the user to decide an appropriate
%value.
addParameter(p,'CoreVid','Off')
addParameter(p,'MovieFileType','Uncompressed AVI')
addParameter(p,'ScaleBarSwitch','Off')

addParameter(p,'RefineROIs','On')
%CG: Switching RefineROIs to 'Off' will cause the first event detected at
%any given location to establish a ROI. Subsequent larger events will not
%reset the established ROI when RefineROIs is 'Off'. This feature was added
%while developing the AARG algorithm and is not intended to be applied by
%the general user.
addParameter(p,'figureFate',1)
%CG: Leave or close figures after analysis. Leaving figures behind after
%analysis can be useful for making publication figures.
addParameter(p, 'osType', '')
parse(p,varargin{:}); param = p.Results; 

% IDString = param.IDString; 
datadirs = param.datadirs;
partnamesuff = param.partnamesuff; ROI_sl = param.ROI_sl; 
% CoreOneImage = param.CoreOneImage; CoreTwoImage = param.CoreTwoImage;
CoreVid = param.CoreVid; MovieFileType = param.MovieFileType;
ScaleBarSwitch = param.ScaleBarSwitch; RefineROIs = param.RefineROIs; PDE = param.PDE;
figFate = param.figureFate; osType = param.osType;

if isempty(PDE)
    disp('Post-Detection Exclusion needs to be set during AARG_CoreHandler call')
    datadirs = 1;
end
% if nargin == 0
suffixStrs = param.suffixStrs;
% end
warning('off', 'all')
%CG: AARG_CoreHandler will attempt to add a path which may or may not exist under a
%'try' within a try-catch loop. To avoid unnecessary user concern or
%confusion, warnings have been turned off. 
fmedianRGB_blank = []; fmedianRGB = [];

if strcmp(RefineROIs, 'Off')
    disp(strcat('WARNING: AARG will not reposition ROIs according to CoM of largest event. R',...
        'edefine RefineROIs variable if appropriate'))
end

if ~strcmp(RefineROIs, 'Off')
    strinsert = '_';
else
    strinsert = '_nrr_';
end

%CG: control size and position of the figure containing image of cell.
left_expfactor = 0.25; bot_expfactor = -0.5;
width_expfactor = 0.4; height_expfactor = 0.4;

%CG: Parts2Take is for reducing the amount of data processed/analysed. 
Parts2Take = [];

screenSize = get(0,'ScreenSize'); screenWidth = screenSize(3); screenHeight = screenSize(4);

if iscell(datadirs)
    Sz_datadirs = size(datadirs); NumDirs = Sz_datadirs(1)*Sz_datadirs(2);
    for cDir = 1 : NumDirs
        datadir = datadirs{cDir}; [~, FileName, ~] = fileparts(datadir);
        
        if ~isempty(strfind(osType, 'MAC')) 
            Slashes = strfind(datadir, '/');
        else 
            Slashes = strfind(datadir, '\');
        end 
        
        tgtdir = datadir(1:Slashes(end)-1); cd(tgtdir); addpath(tgtdir)
    
        cd(datadir); dirInfo = dir; NumItems = size(dirInfo,1); TDFolders = {};
        tdfc = 0; Sz_suffixStrs = size(suffixStrs);
        NumExpr = strcat(partnamesuff, '\d*'); suffixStrs_copy = suffixStrs;
        for cCell = 1 : Sz_suffixStrs(2)
            for cItemIdx = 1 : NumItems
                cItem = dirInfo(cItemIdx).name;
                
                if numel(cItem) > numel(datadir(Slashes(end)+1:end))
                    partnameExpr = regexp(cItem, NumExpr, 'match');
                    suffixID = suffixStrs_copy{cCell};
                    TDFolderName = strcat(suffixID, '_ThresholdData');
                    TDFolderName_ext = strcat(suffixID, partnameExpr, '_ThresholdData');
                    if iscell(TDFolderName_ext) && ~isempty(TDFolderName_ext)
                        TDFolderName_ext = TDFolderName_ext{1};
                    end
                    if iscell(TDFolderName) && ~isempty(TDFolderName)
                        TDFolderName = TDFolderName{1};
                    end

                    if ~isempty(strfind(cItem, TDFolderName)) && ~isempty(strfind(cItem, '_ThresholdData'))
                        tdfc = tdfc + 1; TDFolders{tdfc,1} = cItem;
                    elseif ~isempty(partnameExpr)
                        if ~isempty(strfind(cItem, TDFolderName_ext))
                            tdfc = tdfc + 1; TDFolders{tdfc,1} = cItem;
                        end
                    end
                    if ~isempty(strfind(cItem, 'fmedian')); ValidfmedianFile = cItem; end
                end
            end
            suffixStrs_copy{cCell} = 'x';
        end
%CG: the TDFolders contents must be sorted. Each name should be sorted
%according to the suffixStr and those not containing "partnamesuff"
%probably need to go at the top.

        TDFolders_sorted = cell(tdfc, 1);
        dExpr1 = strcat('(?<=', partnamesuff, ')\d');
        dExpr2 = '\d*(?=_ThresholdData)'; NameGroupSz = tdfc/Sz_suffixStrs(2);
        for cCell = 1 : Sz_suffixStrs(2)
            for cItemIdx = 1 : tdfc
%CG: an embedded for-loop is necessary here, in case the files have been
%partitioned during acquisition. In which case, the size of TDFolders will
%not be equal to the size of suffixStrs. 
                cItem2sort = TDFolders{cItemIdx};
                partnameExpr1 = regexp(cItem2sort, dExpr1, 'match');
                partnameExpr2 = regexp(cItem2sort, dExpr2, 'match');

                if isempty(partnameExpr1) && isempty(partnameExpr2)
                    partnameNum2 = [];
                elseif ~isempty(partnameExpr1) && ~isempty(partnameExpr2)
                    partnameNum2 = partnameExpr2{1};
                elseif isempty(partnameExpr1) && ~isempty(partnameExpr2)
                    partnameNum2 = [];
                end

                if isempty(partnameNum2) && ~isempty(strfind(cItem2sort, strcat(suffixStrs{cCell}, '_ThresholdData')))
                    SortIdx = (NameGroupSz*(cCell-1))+1;
                    TDFolders_sorted{SortIdx, 1} = cItem2sort;
                elseif ~isempty(partnameNum2) && ~isempty(strfind(cItem2sort, strcat(suffixStrs{cCell}, partnamesuff)))

                    Nulls = strfind(partnameNum2, '0');
                    if ~isempty(Nulls); partnameNum2(1:Nulls(end)) = []; end
                    
                    pnOrder = str2double(partnameNum2);
                    if cCell == 1; SortIdx = pnOrder;
                    else SortIdx = (NameGroupSz*(cCell-1))+pnOrder;
                    end
                    TDFolders_sorted{SortIdx, 1} = cItem2sort;
                end
            end
        end

        [~, ExpName, ~] = fileparts(datadir); 

%CG: If there are multiple ThresholdData folders, then it will be necessary
%to loop through the folders with the appropriate AARG_Core function. 

        ThresDataList = {};

        for cCell = 1 : Sz_suffixStrs(2)
            tdfcounter = 0;
            for cDFIdx = 1 : tdfc
                cDF = TDFolders_sorted{cDFIdx, 1}; dExpr = strcat('(?<=', partnamesuff, ')\d');
                Targetstr = strcat(ExpName, suffixStrs{cCell}, partnamesuff, regexp(cDF, dExpr, 'match'), '_ThresholdData');
                if ~isempty(Targetstr) && iscell(Targetstr)
                    Targetstr = Targetstr{1};
                end
                if strfind(cDF, strcat(ExpName, suffixStrs{cCell}, '_ThresholdData'))
                    tdfcounter = tdfcounter + 1;
                elseif ~isempty(Targetstr)
                    if strfind(cDF,Targetstr); tdfcounter = tdfcounter + 1; end
                end
            end
            ThresDataList{cCell, 1} = suffixStrs{cCell}; ThresDataList{cCell, 2} = tdfcounter;
        end
        Cellc = 1; Sz_ThresDataList = size(ThresDataList);
        chunkcount = 0;
        for cFolderIdx = 1 : tdfc 
            ZeroList = 0; cFolder = TDFolders_sorted{cFolderIdx};
            TDLoc = strfind(cFolder, '_ThresholdData');
            TDDir = strcat(datadir, '/', cFolder);
            addpath(TDDir); cd(TDDir)

            if cFolderIdx == 1
                ZeroList = 1; Cellc = Cellc + 1;
            elseif ~mod((cFolderIdx-1),ThresDataList{Cellc, 2}) && cFolderIdx ~= tdfc 
                if ~isempty(strfind(cFolder, strcat(ExpName, suffixStrs{Cellc}, '_ThresholdData')))
                    Cellc = Cellc + 1; ZeroList = 1;
                end
            elseif ~mod((cFolderIdx-1),ThresDataList{Cellc, 2}) && cFolderIdx == tdfc 
                ZeroList = 1;
            end
            if Sz_ThresDataList(1) < Cellc; Cellc = Sz_ThresDataList(1); end

            try
                BreakTry
                load(strcat(cFolder, ''), 'TDFileList', 'TotalChunkSize', 'TotalCCSList');
            catch
                dirInfo = dir; NumFiles = size(dirInfo,1); TDFileList = {};
                TDFileC = 0;

                for cFileIdx = 1 : NumFiles

                    cFile = dirInfo(cFileIdx).name; Sz_cFile = size(cFile);

                    if Sz_cFile(2) > 2
                        if strfind(cFile, 'Thresholds.mat')
                            TDFileC = TDFileC + 1;
                            for cCellIdx = 1 : Sz_suffixStrs(2)
                                if ~isempty(strfind(cFile, strcat(suffixStrs{cCellIdx}, '_Chunk')))
                                    TDFileList{TDFileC, 1} = cFile;
                                elseif ~isempty(strfind(cFile, strcat(cFolder(1:TDLoc-1), '_Chunk')))
                                    TDFileList{TDFileC, 1} = cFile;
                                end
                            end
                        end
                    end
                end
%CG: (10Feb21) if the chunk label hits 10 then the file order gets messed up. To
%ensure correct order of the files, the cells need to be indexed and
%sorted. 

%CG: (10Feb21) grab the chunk number from the file name to create the
%index. 
                numberOfDataFiles = size(TDFileList,1); 
                indexDataFileList = zeros(numberOfDataFiles,1);
                indexDataFileList(:) = NaN;
                for cFileIdx = 1 : numberOfDataFiles
                    numstr = regexp(TDFileList{cFileIdx}, '(?<=Chunk_)\d*', 'match');
                    indexDataFileList(cFileIdx,1) = str2double(numstr{1});
                end
%                 indexDataFileList = cell2mat(indexDataFileCellArray); 
                [~, idices] = sort(indexDataFileList, 'ascend');
                TDFileList = TDFileList(idices);
                
                Sz_TDFileList = size(TDFileList); NumIDs = Sz_TDFileList(1)*Sz_TDFileList(2);
                TotalChunkSize = 0;

                if ZeroList == 1
                    LastcChunkSize = 0; TotalCCSList = [];
                else
                    LastcChunkSize = TotalCCSList(end); TotalCCSList = [];
                end
                ChunkSizeList = [];

                for cIDIdx = 1 : NumIDs

                    mm = matfile(TDFileList{cIDIdx}); ChunkSize = mm.ChunkSize;
                    chunkcount = chunkcount + ChunkSize;
                    TotalChunkSize = TotalChunkSize + ChunkSize;
                    if cIDIdx == 1 
                        ChunkSizeList(1, 1) = ChunkSize;
                        TotalCCSList(1, 1) = LastcChunkSize + 1;
                        TotalCCSList(2, 1) = ChunkSize + (LastcChunkSize + 1);
                    elseif cIDIdx == NumIDs
                        ChunkSizeList(end+1, 1) = ChunkSize;
                        TotalCCSList(end+1, 1) = ChunkSizeList(cIDIdx, 1) + TotalCCSList(cIDIdx, 1) + 1;
                        TotalCCSList(end) = TotalCCSList(end)-1;
                    else
                        ChunkSizeList(end+1, 1) = ChunkSize;
                        TotalCCSList(end+1, 1) = ChunkSizeList(cIDIdx, 1) + TotalCCSList(cIDIdx, 1);
                    end
                end
                save(strcat(cFolder, '_FileOrder.mat'), 'TDFileList', 'TotalChunkSize', 'TotalCCSList');
            end
            
        end
        
        NumTDFs_baseline = ThresDataList{1, 2};
        
        Sz_ThresDataList = size(ThresDataList);
        nExps = Sz_ThresDataList(1);
%CG: nBaseExps is the number of experiments except the baseline.
%AARG_Core needs to loop through all the non-baseline experiments.

        fileorder = load(strcat(TDFolders_sorted{NumTDFs_baseline}, '_FileOrder.mat'));
        TotalCCSList = fileorder.TotalCCSList; gTotalChunkSize = TotalCCSList(end);
        cd(datadir)
        LockOutc1 = 0; %largestMajAxLen = []; largestEvent_andIdx = [];
        %testMtx = zeros(170,170);
        %%%%%%testing
        %maxMaximaProjImg = zeros(170,170);
        %%%%%%testing  
        suppMat = []; suppMat_ROIs = []; suppCA = {}; suppCA_zSpace = {};
  
%%%%%%testing 
%         NumTDFs_baseline = 1;
        StartVal = 1;
%CG: use StartVal to select which split file to work on. 
%%%%%%testing   
        if (NumTDFs_baseline == 1 || StartVal ~= 1) && ~isempty(Parts2Take); 
            disp('WARNING: Not all data being used!!'); 
        end
        
        for cFolderIdx = StartVal : NumTDFs_baseline 
            cFolder = TDFolders_sorted{cFolderIdx}; TDDir = strcat(datadir, '/', cFolder);
            cd(TDDir)
            if cFolderIdx == StartVal
                if strcmp(RefineROIs, 'Off')
                    RefineOffDir = strcat(TDDir, '/RefineROIsOff_',FileName, suffixStrs{1});
                    try
                        cd(RefineOffDir)
                    catch
                        mkdir(strcat('RefineROIsOff_',FileName, suffixStrs{1}));
                        addpath(RefineOffDir)
                    end
                else
                    RefineOffDir = '';
                end
                OriginTD = TDDir;
            end
            load(strcat(cFolder, '_FileOrder.mat'), 'TDFileList');
            DefString = strcat(ExpName, suffixStrs{1});
            if cFolderIdx == StartVal
                CA_ROIs = cell (1);
                CA_NumEvents = cell(1);
                CA_EventSize = cell(1);
                CA_FrameIEI = cell(1);
                CA_Indices = cell(1);
                CA_MaximaSpread = cell(1);
                CA_PatchHanROIs = cell(1);
                CA_PatchHanEventSize = cell(1);
                
                CA_PreEstROIs = {}; TotalPreCASize = []; Sz_CA_PreEstROIs = [];
                CA_PreEstES = {};
                
                subCA_EventSize = cell(1);
                subCA_Indices = cell(1);
                subCA_MaximaSpread = cell(1);

                ROI_Template_Mtx = []; timestepTemplate_Mtx = [];
                compTemplate_Mtx = []; zzTemplate_Mtx = [];
            end
            
            try 
                if strcmp(RefineROIs, 'Off')
                    cd(RefineOffDir)
                end
                if cFolderIdx == StartVal && LockOutc1 == 0
                    load(strcat(TDDir, '/', 'CAs', strinsert, 'Core_', DefString, '.mat'), 'CA_ROIs', 'CA_EventSize',...
                        'CA_NumEvents', 'CA_FrameIEI','CA_Indices','ProperCASize');
                    if figFate == 0
                        figH_ROIs4Show = figure('Name', strcat(DefString, ': Established ROIs'), 'NumberTitle', 'Off');
                        data = load(TDFileList{1}); fmedian = data.fmedian;
                        figure(figH_ROIs4Show); logFluoImage=log10(fmedian);                                                
                        logFluoImage=logFluoImage-min(logFluoImage(:));                             
                        logFluoImage=logFluoImage/max(logFluoImage(:));                             
                        logFluoImage=uint8(logFluoImage*256);                                       
                        cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);        
                        colormap(cyanColorMap);                                                     
                        fmedianRGB_blank=ind2rgb(logFluoImage,cyanColorMap); 
                        Sz_CA_ROIs = size(CA_ROIs); TCs = Sz_CA_ROIs(1)*Sz_CA_ROIs(2);
                        for cc = 1 : TCs
                            if ~isempty(CA_ROIs{cc})
                                fmedianRGB_blank(CA_ROIs{cc}) = 1;
                            end
                        end
                        imshow(fmedianRGB_blank)
                    end
                    LockOutc1 = 1;
                else
                    BreakTry
                end
                
            catch

                if LockOutc1 == 0
                    if cFolderIdx == StartVal
                        EstROIs_figH = figure('Name', strcat(DefString, ': Established ROIs (baseline)'), 'NumberTitle', 'Off');
                        h = waitbar(0, '');
%                         FigPos = EstROIs_figH.Position;
%                         left_exp = screenWidth*left_expfactor; bot_exp = screenHeight*bot_expfactor;
                        
%                         left_exp = FigPos(1)*-left_expfactor; bot_exp = FigPos(2)*-bot_expfactor;
%                         width_exp = FigPos(1)*width_expfactor; height_exp = FigPos(2)*height_expfactor;
%                         NewFigPos = FigPos+[left_exp, bot_exp, width_exp, height_exp];
                        [NewFigPos] = roiDisplayPosition(screenWidth,screenHeight);
                    end
                    load(strcat(cFolder, '_FileOrder.mat'), 'TDFileList', 'TotalChunkSize', 'TotalCCSList');
                    
                    if ~isempty(Parts2Take)
                        TDFileList_2 = TDFileList;TDFileList={};
                        TotalCCSList_2 = TotalCCSList;
                        TotalChunkSize = 0;
                        TotalCCSList = TotalCCSList_2(Parts2Take(1));
                        for cPart = 1 : numel(Parts2Take)
                            TDFileList{cPart,1} = TDFileList_2{Parts2Take(cPart)};
                            TotalCCSList = [TotalCCSList;TotalCCSList_2(Parts2Take(cPart)+1)];
                            TotalChunkSize = TotalChunkSize + ...
                                TotalCCSList_2(Parts2Take(cPart)+1)-TotalCCSList_2(Parts2Take(cPart));
                        end
                    else
                        TotalCCSList_2 = [];
                    end
                    timestep = TotalCCSList(1);
                    coi = struct;
%CG: coi = Core One Input. Structure used to help manage CoreOne input
%arguments.
                    coi.ROI_sl = ROI_sl;
                    coi.gTotalChunkSize = gTotalChunkSize;
                    coi.TotalCCSList = TotalCCSList;
                    coi.TotalCCSList_2 = TotalCCSList_2;
                    coi.DefString = DefString;
                    coi.TDFileList = TDFileList;
                    coi.CoreVid = CoreVid;
                    coi.MovieFileType = MovieFileType;
                    coi.timestep = timestep;
                    coi.cFolderIdx = cFolderIdx;
                    coi.StartVal = StartVal;
                    coi.EstROIs_figH = EstROIs_figH;
                    coi.NewFigPos = NewFigPos;
                    coi.ROI_Template_Mtx = ROI_Template_Mtx;
                    coi.timestepTemplate_Mtx = timestepTemplate_Mtx;
                    coi.compTemplate_Mtx = compTemplate_Mtx;
                    coi.zzTemplate_Mtx = zzTemplate_Mtx;
                    coi.CA_ROIs = CA_ROIs;
                    coi.CA_EventSize = CA_EventSize;
                    coi.CA_NumEvents = CA_NumEvents;
                    coi.CA_FrameIEI = CA_FrameIEI;
                    coi.CA_Indices = CA_Indices;
                    coi.CA_MaximaSpread = CA_MaximaSpread;
                    coi.CA_PatchHanROIs = CA_PatchHanROIs;
                    coi.CA_PatchHanEventSize = CA_PatchHanEventSize;
                    coi.fmedianRGB_blank = fmedianRGB_blank;
                    coi.fmedianRGB = fmedianRGB;
                    coi.ScaleBarSwitch = ScaleBarSwitch;
                    coi.h = h;
                    coi.RefineROIs = RefineROIs;
                    coi.datadir = datadir;
                    coi.TDDir = TDDir;
%                     coi.largestEvent_andIdx = largestEvent_andIdx;
%                     coi.largestMajAxLen = largestMajAxLen;
%                     coi.testMtx = testMtx;
                    coi.CA_PreEstROIs = CA_PreEstROIs;
                    coi.TotalPreCASize = TotalPreCASize;
                    coi.Sz_CA_PreEstROIs = Sz_CA_PreEstROIs;
                    coi.CA_PreEstES = CA_PreEstES;
                    coi.suppMat = suppMat;
                    coi.suppMat_ROIs = suppMat_ROIs;
                    coi.suppCA = suppCA;
                    coi.suppCA_zSpace = suppCA_zSpace;
                    coi.subCA_EventSize = subCA_EventSize;
                    coi.subCA_Indices = subCA_Indices;
                    coi.subCA_MaximaSpread = subCA_MaximaSpread;
                    coi.IEI = PDE;
                    coi.startpnt = NaN;
                    
                    %%%%%testing
                    %coi.maxMaximaProjImg = maxMaximaProjImg;
                    %coi.subsidROIList = subsidROIList;
                    %%%%%testing
                    [coo] = AARG_Core(coi);
%CG: coo = Core One Output.                                      
                    h = coo.h;
                    CA_ROIs = coo.CA_ROIs;
                    CA_NumEvents = coo.CA_NumEvents;
                    CA_EventSize = coo.CA_EventSize;
                    CA_FrameIEI = coo.CA_FrameIEI;
                    CA_Indices = coo.CA_Indices;
                    CA_MaximaSpread = coo.CA_MaximaSpread;
                    ProperCASize = coo.ProperCASize;
                    EstROIs_ax = coo.EstROIs_ax;
                    CA_PatchHanROIs = coo.CA_PatchHanROIs;
                    CA_PatchHanEventSize = coo.CA_PatchHanEventSize;
                    EstROIs_figH = coo.EstROIs_figH;
                    ROI_Template_Mtx = coo.ROI_Template_Mtx;
                    timestepTemplate_Mtx = coo.timestepTemplate_Mtx;
                    compTemplate_Mtx = coo.compTemplate_Mtx;
                    zzTemplate_Mtx = coo.zzTemplate_Mtx;
                    CA_PreEstROIs = coo.CA_PreEstROIs;
                    TotalPreCASize = coo.TotalPreCASize;
                    Sz_CA_PreEstROIs = coo.Sz_CA_PreEstROIs;
                    fmedianRGB_blank = coo.fmedianRGB_blank;
                    fmedianRGB = coo.fmedianRGB;
%                     largestEvent_andIdx = coo.largestEvent_andIdx;
%                     largestMajAxLen = coo.largestMajAxLen;
                    suppMat = coo.suppMat;
                    suppCA_zSpace = coo.suppCA_zSpace;
                    suppMat_ROIs = coo.suppMat_ROIs;
                    suppCA = coo.suppCA;
                    subCA_EventSize = coo.subCA_EventSize;
                    subCA_Indices = coo.subCA_Indices;
                    subCA_MaximaSpread = coo.subCA_MaximaSpread;
                    %%%%%testing
%                     maxMaximaProjImg = coo.maxMaximaProjImg;
                    %%%%%testing
                end
            end

        end
       
        
        if strcmp(RefineROIs, 'Off'); cd(RefineOffDir); else cd(OriginTD); end

        if LockOutc1 == 0
            
            if figFate == 0
                Sz_CA_ROIs = size(CA_ROIs); TCs = Sz_CA_ROIs(1)*Sz_CA_ROIs(2);
                for cc = 1 : TCs
                    if ~isempty(CA_ROIs{cc})
                        fmedianRGB_blank(CA_ROIs{cc}) = 1;
                    end
                end
                figH_ROIs4Show = figure('Name', strcat(DefString, ': Established ROIs'), 'NumberTitle', 'Off');
                imshow(fmedianRGB_blank)
            end
            
            close(h); close(EstROIs_figH)
            
            save(strcat('CAs', strinsert, 'Core_', DefString, '.mat'), 'CA_ROIs', 'CA_EventSize',...
                   'CA_NumEvents', 'CA_FrameIEI', 'CA_Indices', 'CA_MaximaSpread', 'ProperCASize',...
                   'subCA_EventSize','subCA_Indices','subCA_MaximaSpread');

        end 

        if iscell(ProperCASize); ProperCASize = size(ProperCASize); end
        addpath(OriginTD);    
        strinsert = '_';
%CG: strinsert is only needed for first AARG_Core call. From this point it
%should be defined as underscore.
        NumTDFs_baseline = ThresDataList{1, 2}; NumTDFs_nbaseline = NumTDFs_baseline; 
        cumulTD = [];
        
        for cnbExp = 2 : nExps

            CA_PreEstROIs = CA_ROIs; CA_PreEstES = CA_EventSize;
%CG: 'CA_PreEstROIs' contains all the pre-established ROIs and gets
%up-dated with every condition (i.e. after each cycle of the cnbExp = 2 :
%nExps for-loop
            NumTDFs_nbaseline = sum([ThresDataList{cnbExp, 2},NumTDFs_nbaseline]);
            cFolder = TDFolders_sorted{NumTDFs_nbaseline};
            cd(strcat(datadir, '/', cFolder)) 
            load(strcat(cFolder, '_FileOrder.mat'), 'TotalCCSList');
            gTotalChunkSize = TotalCCSList(end);

            LockOutc2 = 0;
%CG: Second call for AARG_Core takes the experimental data. Red ROIs depict ROIs that were
%established in by the first call of AARG_Core or a preceding condition. Yellow ROIs are ones
%appearing at locations that were silent in the baseline data or preceding condition. 
            if isempty(cumulTD)
                cumulTD = ThresDataList{cnbExp, 2} + NumTDFs_baseline;
                startpnt = NumTDFs_baseline;
            else 
                cumulTD = ThresDataList{cnbExp, 2} + cumulTD;
                startpnt = cFolderIdx;
            end
%             largestEvent_andIdx = [];largestMajAxLen = [];
%%%%%%%testing
%             cumulTD = 4;
            StartVal = startpnt+1;
%%%%%%%testing  
            if StartVal ~= startpnt+1; 
                disp('WARNING: Not all data being used!!'); 
            end
            
            for cFolderIdx = StartVal : cumulTD
%CG: the second column of the ThresDataList will give the number total
%number of raw data files for each condition (with the condition name being
%stored in the 1st column). When the acquisition was large and the
%acquisition software automatically partitions the data file, the numbers
%stored in the cells of the second column will be >1.
                cFolder = TDFolders_sorted{cFolderIdx}; cd(strcat(datadir, '/', cFolder))
                if cFolderIdx == startpnt+1
                    c4_start = cFolderIdx; OriginTD = strcat(datadir, '/', cFolder);
                end
                load(strcat(cFolder, '_FileOrder.mat'), 'TDFileList');
                TDDir = strcat(datadir, '/', cFolder); cd(TDDir);

                if cFolderIdx == startpnt+1

                    DefString = strcat(ExpName, suffixStrs{cnbExp});
                    CA_ROIs = CA_PreEstROIs;
                    CA_NumEvents = cell(1);
                    CA_EventSize = CA_PreEstES;
                    CA_FrameIEI = cell(1);
                    CA_Indices = cell(1);
                    CA_MaximaSpread = cell(1);
                    CA_PatchHanROIs = cell(1);
                    CA_PatchHanEventSize = cell(1);
                    
                    subCA_EventSize = cell(1);
                    subCA_Indices = cell(1);
                    subCA_MaximaSpread = cell(1);
                    
                    ROI_Template_Mtx = []; timestepTemplate_Mtx = [];
                    compTemplate_Mtx = []; zzTemplate_Mtx = [];
                    Sz_CA_PreEstROIs = size(CA_PreEstROIs);
                    TotalPreCASize = [];

                end
                try 
%                     BreakTry
                    if cFolderIdx == StartVal && LockOutc2 == 0;
                        load(strcat('CAs_', 'Core_', DefString, '.mat'), 'CA_ROIs',...
                            'CA_EventSize', 'CA_NumEvents', 'CA_FrameIEI', 'ProperCASize');
                        
                        if figFate == 0
%CG: to replicate the image of the cell with established ROIs as a 2D
%image, exactly as it appears after AARG_Core completes, the ROI indices
%need to be changed to red for all the ROIs established before 2nd and all
%subsequent AARG_Core calls (i.e. pre-Established ROIs). All the indices
%covered by ROIs established in the most recent call of AARG_Core should be
%changed to yellow.
                            if numel(Sz_CA_PreEstROIs) == 2
                                TCs_preEst = Sz_CA_PreEstROIs(1)*Sz_CA_PreEstROIs(2);
                            elseif numel(Sz_CA_PreEstROIs) == 3
                                TCs_preEst = Sz_CA_PreEstROIs(1)*Sz_CA_PreEstROIs(2)*Sz_CA_PreEstROIs(3);
                            end

                            for cc = 1 : TCs_preEst
                                if ~isempty(CA_PreEstROIs{cc})
                                    fmedianRGB_blank(CA_PreEstROIs{cc}) = 1;
                                end
                            end

                            [MaxPreEst_row, MaxPreEst_col, MaxPreEst_z] = ind2sub(Sz_CA_PreEstROIs, TCs_preEst);

                            Sz_CA_ROIs = size(CA_ROIs);
                            if numel(Sz_CA_ROIs) == 2
                                TCs = Sz_CA_ROIs(1)*Sz_CA_ROIs(2);
                            elseif numel(Sz_CA_ROIs) == 3
                                TCs = Sz_CA_ROIs(1)*Sz_CA_ROIs(2)*Sz_CA_ROIs(3);
                            end
                            Sz_fmedianRGB_blank = size(fmedianRGB_blank);
                            for cc = 1 : TCs
%CG: When there are pre-established ROIs to take into consideration, it is
%essential the cells are interrogated using subscripts. the CA_PreEstROI
%cell array will have a different size compared to the CA_ROIs cell array
%and indices mean different cells, when arrays have different sizes. 
                                [c_row, c_col, c_z] = ind2sub(Sz_CA_ROIs, cc);
                                if c_row <= MaxPreEst_row && c_col <= MaxPreEst_col && c_z <= MaxPreEst_z
                                    if ~isempty(CA_ROIs{c_row, c_col, c_z})
                                        [rsubs, csubs, ~] = ind2sub(Sz_fmedianRGB_blank, CA_ROIs{c_row, c_col, c_z});

                                        if isempty(CA_PreEstROIs{c_row, c_col, c_z})
                                            zsubs = ones(1, numel(rsubs));
                                            zsubs = zsubs.*2;

                                            TwoDidces = sub2ind(Sz_fmedianRGB_blank, rsubs, csubs, zsubs);
                                            fmedianRGB_blank(CA_ROIs{c_row, c_col, c_z}) = 1;
                                            fmedianRGB_blank(TwoDidces) = 1;
                                        end
                                    end
                                else
                                    if ~isempty(CA_ROIs{c_row, c_col, c_z})

                                        fmedianRGB_blank(CA_ROIs{c_row, c_col, c_z}) = 1;
        %CG: if the same cell is empty in CA_PreEstROIs, this means that the current ROI,
        %was established in Core Three and should appear yellow. 
                                        [rsubs, csubs, ~] = ind2sub(Sz_fmedianRGB_blank, CA_ROIs{c_row, c_col, c_z});
                                        zsubs = ones(1, numel(rsubs)); zsubs = zsubs.*2;
                                        TwoDidces = sub2ind(Sz_fmedianRGB_blank, rsubs, csubs, zsubs);
                                        fmedianRGB_blank(TwoDidces) = 1;

                                    end
                                end
                            end


                            figH_ROIs4Show = figure('Name', strcat(DefString, ': Established ROIs'),...
                                'NumberTitle','off');
                            imshow(fmedianRGB_blank)
        %                     saveas(figH_ROIs4Show, strcat(DefString, '_Core'), 'pdf')
        %                     saveas(figH_ROIs4Show, strcat(DefString, '_Core'), 'fig')

        %                     saveas(EstROIs_figH, strcat(DefString, '_Core'), 'pdf')
        %                     saveas(EstROIs_figH, strcat(DefString, '_Core'), 'fig')
                        end
                        LockOutc2 = 1;
                    else
                        BreakTry
                    end
                    
                catch

                    if LockOutc2 == 0
                        if cFolderIdx == c4_start
                            EstROIs_figH = figure('Name', strcat(DefString, ': Established ROIs'), 'NumberTitle', 'Off');
                            h = waitbar(0, '');
                            
                            FigPos = EstROIs_figH.Position;
                            left_exp = FigPos(1)*-left_expfactor; bot_exp = FigPos(2)*-bot_expfactor;
                            width_exp = FigPos(1)*width_expfactor; height_exp = FigPos(2)*height_expfactor;
                            NewFigPos = FigPos+[left_exp, bot_exp, width_exp, height_exp];
                        end
                        load(strcat(cFolder, '_FileOrder.mat'), 'TDFileList_sorted', 'TotalChunkSize', 'TotalCCSList');

                        if ~isempty(Parts2Take)
                            TDFileList_2 = TDFileList;TDFileList={};
                            TotalCCSList_2 = TotalCCSList;
                            TotalChunkSize = 0;
                            TotalCCSList = TotalCCSList_2(Parts2Take(1));
                            for cPart = 1 : numel(Parts2Take)
                                TDFileList{cPart,1} = TDFileList_2{Parts2Take(cPart)};
                                TotalCCSList = [TotalCCSList;TotalCCSList_2(Parts2Take(cPart)+1)];
                                TotalChunkSize = TotalChunkSize + ...
                                    TotalCCSList_2(Parts2Take(cPart)+1)-TotalCCSList_2(Parts2Take(cPart));
                            end
                        else
                            TotalCCSList_2 = [];
                        end
                        timestep = TotalCCSList(1);
                        
                        c2i = struct;
                        c2i.ROI_sl = ROI_sl;
                        c2i.gTotalChunkSize = gTotalChunkSize;
                        c2i.TotalCCSList = TotalCCSList;
                        c2i.TotalCCSList_2 = TotalCCSList_2;
                        c2i.TDFileList = TDFileList;
                        c2i.DefString = DefString;
                        c2i.CoreVid = CoreVid;
                        c2i.MovieFileType = MovieFileType;
                        c2i.timestep = timestep;
                        c2i.cFolderIdx = cFolderIdx;
                        c2i.StartVal = StartVal;
                        c2i.EstROIs_figH = EstROIs_figH;
                        c2i.NewFigPos = NewFigPos;
                        c2i.h = h;
                        c2i.ROI_Template_Mtx = ROI_Template_Mtx;
                        c2i.timestepTemplate_Mtx = timestepTemplate_Mtx;
                        c2i.compTemplate_Mtx = compTemplate_Mtx;
                        c2i.zzTemplate_Mtx = zzTemplate_Mtx;
                        c2i.CA_ROIs = CA_ROIs;
                        c2i.CA_EventSize = CA_EventSize;
                        c2i.CA_NumEvents = CA_NumEvents;
                        c2i.CA_FrameIEI = CA_FrameIEI;
                        c2i.CA_Indices = CA_Indices;
                        c2i.CA_MaximaSpread = CA_MaximaSpread;
                        c2i.CA_PatchHanROIs = CA_PatchHanROIs;
                        c2i.CA_PatchHanEventSize = CA_PatchHanEventSize;
                        c2i.CA_PreEstROIs = CA_PreEstROIs;
                        c2i.TotalPreCASize = TotalPreCASize;
                        c2i.Sz_CA_PreEstROIs = Sz_CA_PreEstROIs;
                        c2i.CA_PreEstES = CA_PreEstES;
                        c2i.fmedianRGB_blank = fmedianRGB_blank;
                        c2i.fmedianRGB = fmedianRGB;
                        c2i.ScaleBarSwitch = ScaleBarSwitch;
                        c2i.TDDir = TDDir;
%                         c2i.largestEvent_andIdx = largestEvent_andIdx;
%                         c2i.largestMajAxLen = largestMajAxLen;
                        c2i.suppMat = suppMat;
                        c2i.suppMat_ROIs = suppMat_ROIs;
                        c2i.suppCA = suppCA;
                        c2i.suppCA_zSpace = suppCA_zSpace;
                        c2i.subCA_Indices = subCA_Indices;
                        c2i.subCA_EventSize = subCA_EventSize;
                        c2i.subCA_MaximaSpread = subCA_MaximaSpread;
                        c2i.IEI = PDE;
                        c2i.startpnt = startpnt;

                        [c2po] = AARG_Core(c2i);

                        CA_ROIs = c2po.CA_ROIs;
                        CA_EventSize = c2po.CA_EventSize;
                        CA_NumEvents = c2po.CA_NumEvents;
                        CA_FrameIEI = c2po.CA_FrameIEI;
                        CA_Indices = c2po.CA_Indices;
                        CA_MaximaSpread = c2po.CA_MaximaSpread;
                        ProperCASize = c2po.ProperCASize;
                        EstROIs_ax = c2po.EstROIs_ax;
                        CA_PatchHanROIs = c2po.CA_PatchHanROIs;
                        CA_PatchHanEventSize = c2po.CA_PatchHanEventSize;
                        EstROIs_figH = c2po.EstROIs_figH;
                        h = c2po.h;
                        ROI_Template_Mtx = c2po.ROI_Template_Mtx;
                        timestepTemplate_Mtx = c2po.timestepTemplate_Mtx;
                        compTemplate_Mtx = c2po.compTemplate_Mtx;
                        zzTemplate_Mtx = c2po.zzTemplate_Mtx;
                        CA_PreEstROIs = c2po.CA_PreEstROIs;
                        TotalPreCASize = c2po.TotalPreCASize;
                        Sz_CA_PreEstROIs = c2po.Sz_CA_PreEstROIs;
                        fmedianRGB_blank = c2po.fmedianRGB_blank;
                        fmedianRGB = c2po.fmedianRGB;
%                         largestEvent_andIdx = c2po.largestEvent_andIdx;
%                         largestMajAxLen = c2po.largestMajAxLen;
                        suppMat = c2po.suppMat;
                        suppMat_ROIs = c2po.suppMat_ROIs;
                        suppCA = c2po.suppCA;
                        suppCA_zSpace = c2po.suppCA_zSpace;
                        subCA_EventSize = c2po.subCA_EventSize;
                        subCA_Indices = c2po.subCA_Indices;
                        subCA_MaximaSpread = c2po.subCA_MaximaSpread;
                    end
                end

            end 

            cd(OriginTD)
            if LockOutc2 == 0
                if figFate == 0
                    
                    data = load(TDFileList{1}); fmedian = data.fmedian;
                    logFluoImage=log10(fmedian);                                                
                    logFluoImage=logFluoImage-min(logFluoImage(:));                             
                    logFluoImage=logFluoImage/max(logFluoImage(:));                             
                    logFluoImage=uint8(logFluoImage*256);                                       
                    cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);        
                    colormap(cyanColorMap);                                                     
                    fmedianRGB_blank=ind2rgb(logFluoImage,cyanColorMap);  
                    if numel(Sz_CA_PreEstROIs) == 2
                        TCs_preEst = Sz_CA_PreEstROIs(1)*Sz_CA_PreEstROIs(2);
                    elseif numel(Sz_CA_PreEstROIs) == 3
                        TCs_preEst = Sz_CA_PreEstROIs(1)*Sz_CA_PreEstROIs(2)*Sz_CA_PreEstROIs(3);
                    end

                    for cc = 1 : TCs_preEst
                        if ~isempty(CA_PreEstROIs{cc})
                            fmedianRGB_blank(CA_PreEstROIs{cc}) = 1;
                        end
                    end
                    
                    [MaxPreEst_row, MaxPreEst_col, MaxPreEst_z] = ind2sub(Sz_CA_PreEstROIs, TCs_preEst);

                    Sz_CA_ROIs = size(CA_ROIs);
                    if numel(Sz_CA_ROIs) == 2
                        TCs = Sz_CA_ROIs(1)*Sz_CA_ROIs(2);
                    elseif numel(Sz_CA_ROIs) == 3
                        TCs = Sz_CA_ROIs(1)*Sz_CA_ROIs(2)*Sz_CA_ROIs(3);
                    end
                    Sz_fmedianRGB_blank = size(fmedianRGB_blank);
                    for cc = 1 : TCs
%CG: When there are pre-established ROIs to take into consideration, it is
%essential the cells are interrogated using subscripts rather than linear
%indices. the CA_PreEstROI cell array will have a different size compared
%to the CA_ROIs cell array and indices mean different cells, when arrays
%have different sizes.
                        [c_row, c_col, c_z] = ind2sub(Sz_CA_ROIs, cc);
                        if c_row <= MaxPreEst_row && c_col <= MaxPreEst_col && c_z <= MaxPreEst_z
                            if ~isempty(CA_ROIs{c_row, c_col, c_z})
                                [rsubs, csubs, ~] = ind2sub(Sz_fmedianRGB_blank, CA_ROIs{c_row, c_col, c_z});

                                if isempty(CA_PreEstROIs{c_row, c_col, c_z})
                                    zsubs = ones(1, numel(rsubs));
                                    zsubs = zsubs.*2;

                                    TwoDidces = sub2ind(Sz_fmedianRGB_blank, rsubs, csubs, zsubs);
                                    fmedianRGB_blank(CA_ROIs{c_row, c_col, c_z}) = 1;
                                    fmedianRGB_blank(TwoDidces) = 1;
                                end
                            end
                        else
                            if ~isempty(CA_ROIs{c_row, c_col, c_z})

                                fmedianRGB_blank(CA_ROIs{c_row, c_col, c_z}) = 1;
%CG: if the same cell is empty in CA_PreEstROIs, this means that the current ROI,
%was established in Core Three and should appear yellow. 
                                [rsubs, csubs, ~] = ind2sub(Sz_fmedianRGB_blank, CA_ROIs{c_row, c_col, c_z});
                                zsubs = ones(1, numel(rsubs)); zsubs = zsubs.*2;
                                TwoDidces = sub2ind(Sz_fmedianRGB_blank, rsubs, csubs, zsubs);
                                fmedianRGB_blank(TwoDidces) = 1;

                            end
                        end
                    end
                    figH_ROIs4Show = figure('Name', strcat(DefString, ': Established ROIs'), 'NumberTitle', 'Off');
                    figure(figH_ROIs4Show); imshow(fmedianRGB_blank)
%                     saveas(figH_ROIs4Show, strcat(DefString, '_Core'), 'pdf')
%                     saveas(figH_ROIs4Show, strcat(DefString, '_Core'), 'fig')

%                     saveas(EstROIs_figH, strcat(DefString, '_Core'), 'pdf')
%                     saveas(EstROIs_figH, strcat(DefString, '_Core'), 'fig')
                end
                close(h)
                close(EstROIs_figH)
                cd(OriginTD)
                save(strcat('CAs_', 'Core_', DefString, '.mat'), 'CA_ROIs', 'CA_EventSize',...
                   'CA_NumEvents', 'CA_FrameIEI', 'CA_Indices', 'CA_MaximaSpread', 'CA_PreEstROIs',...
                   'ProperCASize', 'subCA_EventSize','subCA_Indices','subCA_MaximaSpread');
               
            end
            cd(tgtdir)
        end
    end
%     disp(strcat('analysis of "', FileName, '" is complete'))
end