function [handles,ThresOptsID,ThresholdDataFolderList,strlists] = ...
    Detection(cDir,cDirIdx,suffixStrs,CorrectShiftSwt,E_thres,Ampl_thres,...
    handles,partnamesuff,ChunkLimit,ThresholdDataFolderList,strlists,...
    specifiedConditionIdx,ThresOptsID,CellName)
% Outline
% Detection organizes data files for activity detection by subfunction 
% ApplyCrossCorr. ApplyCrossCorr reads raw data directly
% from tif files and applies cross-correlation thresholding algorithm.
% Suprathreshold events are stored in chunks of approximately 500 frames.

% Testing different thresholds: The user can specify multiple thresholds.
% The results will be stored in the ThresholdOptions folder rather than the
% ThresholdData folder.

% Author: Charlie J Gilbride
% version: 1.0.0

if nargin == 0
    AllDirs = uipickfiles('Prompt', 'Select data folders');
    CorrectShiftSwt = 'Off';
    suffixStrs = {'a_baseline', 'b_TTAP2'};
    E_thres = 5;
    E_Incr = 0;
    Ampl_thres = 3.5;
    AmplIncr = 0.5;
    NumCombos = 1;
    partnamesuff = '-file00'; specifiedConditionIdx = [1]; cDirIdx = 1;
    ChunkLimit = 8; handles = struct;
elseif nargin > 0    
    
    AllDirs = {}; AllDirs{1} = cDir;
    if numel(E_thres) > numel(Ampl_thres)
       NumCombos = numel(E_thres); E_Incr = unique(diff(E_thres));
    else
        NumCombos = numel(Ampl_thres); AmplIncr = unique(diff(Ampl_thres));
    end
end

if NumCombos <= 1
    E_valThreshold = E_thres;
    AmplThreshold = Ampl_thres;
    E_valThresholdList = E_valThreshold;
    Ampl_ThresholdList = AmplThreshold;
    disp(strcat('E value threshold =  "', num2str(E_valThreshold), '"'));
    disp(strcat('Amplitude threshold =  "', num2str(AmplThreshold), '"'));
    ThresOptsID = '';
else
    E_valThreshold = E_thres;
    AmplThreshold = Ampl_thres;
    E_valThresholdList = E_valThreshold;
    Ampl_ThresholdList = AmplThreshold;
    if numel(E_thres) == 1
        DStr_e = 'Fix';
        DStr_ampl = 'Flex';
    elseif numel(Ampl_thres) == 1
        DStr_e = 'Flex';
        DStr_ampl = 'Fix';
    end
%CG: use 'Flex' to incrementally increase Activity and Size Thresholds from
%their '_start' values by 'ActIncr' or 'SizeIncr'. Use 'Fix' to keep one
%threshold all the same size across all comparisons.

    if nargin == 0
        if strcmp(DStr_e, 'Fix') && strcmp(DStr_ampl, 'Flex')
            ThresOptsID = strcat('E_val_', num2str(E_thres),...
                '_Ampl', num2str(Ampl_thres(1)),'-', num2str(Ampl_thres(end)));
        elseif strcmp(DStr_e, 'Flex') && strcmp(DStr_ampl, 'Fix')
            ThresOptsID = strcat('E_val_', num2str(E_thres(1)),'-',num2str(E_thres(end)),...
                '_Ampl', num2str(Ampl_thres));
        elseif strcmp(DStr_e, 'Fix') && strcmp(DStr_ampl, 'Fix')
            ThresOptsID = strcat('E_val_', num2str(E_thres), '_Ampl', num2str(Ampl_thres));
        end
    
        if strcmp(DStr_e, 'Flex')
            E_valThresholdList = (E_thres : E_Incr : (E_Incr*(NumCombos-1)) + E_thres);
        elseif strcmp(DStr_e, 'Fix')
            E_valThresholdList = ones(1, NumCombos).*E_thres;
        end

        if strcmp(DStr_ampl, 'Flex')
            Ampl_ThresholdList = (Ampl_thres : AmplIncr : (AmplIncr*(NumCombos-1)) + Ampl_thres);
        elseif strcmp(DStr_ampl, 'Fix')
            Ampl_ThresholdList = ones(1, NumCombos).*Ampl_thres;
        end
    else
        if strcmp(DStr_e, 'Flex') && strcmp(DStr_ampl, 'Fix')
            Ampl_ThresholdList = ones(1, NumCombos).*Ampl_thres;
        elseif strcmp(DStr_e, 'Fix') && strcmp(DStr_ampl, 'Flex')
            E_valThresholdList = ones(1, NumCombos).*E_thres;
        end

    end
end
if numel(unique(E_valThresholdList)) > 1 && numel(unique(Ampl_ThresholdList)) == 1
    NumThresholds = numel(E_valThresholdList); 
elseif numel(unique(E_valThresholdList)) == 1 && numel(unique(Ampl_ThresholdList)) > 1
    NumThresholds = numel(Ampl_ThresholdList); 
else
    NumThresholds = 1;
end
ThresholdDataFolderStr = '_ThresholdData'; 

if ~isempty(partnamesuff)
    NumExpr = strcat(partnamesuff, '\d*');
else 
    NumExpr = '';
end

NumOfWantedChunks = [];
%CG: make NumOfWantedChunks empty, unless doing code dev.

if isempty(specifiedConditionIdx)
    Sz_suffixStrs = size(suffixStrs);
    NumConditions = Sz_suffixStrs(1)*Sz_suffixStrs(2);
elseif ~isempty(specifiedConditionIdx)
    NumConditions = 1;
end
    
dispstat('', 'init'); statusStr = '';

if nargin > 0 
    if numel(E_valThresholdList) == 1
        if numel(AmplThreshold) > 1
            statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
            '"; Amplitude threshold =  "', num2str(AmplThreshold(1)),'-', num2str(AmplThreshold(end)), '"');
        else
            statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                '"; Amplitude threshold =  "', num2str(AmplThreshold), '"');
        end
    else
        if numel(E_valThreshold) > 1
            statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold(1)),...
                '-',num2str(E_valThreshold(end)),'"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "',...
                CellName,'"');
        else
            statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                '"; Amplitude threshold =  "', num2str(AmplThreshold(1)), '-', num2str(AmplThreshold(end)), '" for "',...
                CellName,'"');
        end
    end
    pause(5)
    handles.text_status.String = sprintf(statusStr);
end

if iscell(AllDirs)
    
    FirstDir = AllDirs{1}; 
    if ~isempty(strfind(handles.osType, 'MAC')) 
        Slashes = strfind(cDir, '/');
    else 
        Slashes = strfind(cDir, '\');
    end 
    LastSlash = max(Slashes); DestDir = FirstDir(1 : LastSlash-1);
    cd(DestDir); [~,ParentFolder,~] = fileparts(DestDir);
    disp(strcat('Processing data for "', ParentFolder, '" group'))
    
    [~, CellName, ~] = fileparts(AllDirs{1}); 
    try
%CG: In theory it is good to save FileNameLists.mat because it is a
%light-weight file but it can take a bit of time for matlab to put it
%together. However, if FileNameLists is wrongly configured and the user
%runs AARGI again, this will throw an error.
        DoNotTry
        
%         lists = load(strcat(CellName, '_FileNameLists.mat'));
%         ExpNameList = lists.ExpNameList;
%         partnamenumList = lists.partnamenumList;
%         NumStrFilesList = lists.NumStrFilesList;
%         numExpNameListFiles = lists.numExpNameListFiles;
% 
%         if ~isempty(specifiedConditionIdx)
%             newExpNameList = {}; newpartnamenumList = {}; newNumStrFilesList = []; 
%             Sz_ExpNameList = size(ExpNameList);
%             for cCellIdx = 1 : Sz_ExpNameList(1)
%                 if ~isempty(strfind(ExpNameList{cCellIdx}, suffixStrs{specifiedConditionIdx}))
%                     newExpNameList{cCellIdx,1} = ExpNameList{cCellIdx};
%                     newpartnamenumList{cCellIdx,1} = partnamenumList{cCellIdx};
%                     newNumStrFilesList = [newNumStrFilesList;NumStrFilesList(cCellIdx)];
%                 end
%             end
%             newnumExpNameListFiles = numel(newNumStrFilesList);
%         end
        
    catch

        ExpNameList = {}; 

        CurrentDir = AllDirs{1}; cd(CurrentDir);
        numExpNameListFiles = 0; dirInfo=dir; numfiles=size(dirInfo,1);                                                   

        dExprA = strcat('(?<=', partnamesuff, ')\d','.tif');
        
        for cFile = 1 : numfiles
            
            currentFile = dirInfo(cFile).name;
            
            cStr = 0; quitWhile = 0; 
            if ~isempty(specifiedConditionIdx) 
                quitWhile = 1; dExprB = strcat(suffixStrs{specifiedConditionIdx},'.tif');
            else
                while cStr < NumConditions && quitWhile == 0
                    cStr = cStr + 1;
                    if numel(currentFile) >= numel(suffixStrs{cStr})
                        if ~isempty(strfind(currentFile,suffixStrs{cStr}))
                            quitWhile = 1; dExprB = strcat(suffixStrs{cStr},'.tif');
                        end
                    end
                end
            end
            if quitWhile == 1
                statusExprA = regexp(currentFile, dExprA, 'match');
                statusExprB = strfind(currentFile,dExprB);

                if ~isempty(statusExprA) || ~isempty(statusExprB)

                    if ~isempty(NumExpr)
                        partnameExpr = regexp(currentFile, NumExpr, 'match');
                    else
                        partnameExpr = '';
                    end
                    if isempty(partnameExpr)

                        ExpName = currentFile;
                        numExpNameListFiles = numExpNameListFiles + 1;
                        ExpNameList{numExpNameListFiles, 1} = ExpName;

                        partnamenumList{numExpNameListFiles, 1} = 1;

                    elseif ~isempty(partnameExpr)  

                        ExpName = currentFile;
                        numExpNameListFiles = numExpNameListFiles + 1;
                        ExpNameList{numExpNameListFiles, 1} = ExpName;

                        partnameStr = partnameExpr{1};
                        Nulls = strfind(partnameStr, '0');
                        if ~isempty(Nulls)
                            partnameStr(1:Nulls(end)) = [];
                        end
                        partnamenumList{numExpNameListFiles, 1} = str2double(partnameStr);

                    end
                end
            end
        end
        ExpNameList_sorted = cell(numExpNameListFiles, 1);
        partnamenumList_sorted = cell(numExpNameListFiles, 1);
        dExpr1 = strcat('(?<=', partnamesuff, ')\d');

        NumStrFilesList = zeros(numExpNameListFiles,1);
        if ~isempty(specifiedConditionIdx)
            ExpNameList_sorted = {};
%             ExpNameList_sorted{1} = ExpNameList{specifiedConditionIdx,1};
            for cItemIdx = 1 : numExpNameListFiles
                fn = ExpNameList{cItemIdx,1};
                if ~isempty(strfind(fn, suffixStrs{specifiedConditionIdx}))
                    NumStrFilesList(cItemIdx,1) = specifiedConditionIdx;
                end
            end
            CondIdx = find(NumStrFilesList(:,1) == specifiedConditionIdx);
            LastCondIdx = CondIdx(end); NumCondIdx = numel(CondIdx);
            for cItemIdx = 1 : numExpNameListFiles

                if ~isempty(strfind(ExpNameList{cItemIdx, 1}, suffixStrs{specifiedConditionIdx}))
                    cItem2sort_ExpName = ExpNameList{cItemIdx, 1};
                    partnameExpr1 = regexp(cItem2sort_ExpName, dExpr1, 'match');
                    if isempty(partnameExpr1); partnameExpr1 = '1'; end

                    if ~isempty(partnameExpr1)
                        cPartNameNum = partnamenumList{cItemIdx, 1};
                        SortIdx = (LastCondIdx-NumCondIdx)+cPartNameNum;
                        ExpNameList_sorted{SortIdx, 1} = ExpNameList{cItemIdx, 1};
                        partnamenumList_sorted{SortIdx, 1} = partnamenumList{cItemIdx, 1};
                    end
                end
            end
            
        else
            for cCell = 1 : NumConditions
                for cItemIdx = 1 : numExpNameListFiles
                    fn = ExpNameList{cItemIdx,1};
                    if ~isempty(strfind(fn, suffixStrs{cCell}))
                        NumStrFilesList(cItemIdx,1) = cCell;
                    end
                end
            end
            for cCell = 1 : NumConditions
                CondIdx = find(NumStrFilesList(:,1) == cCell);
                LastCondIdx = CondIdx(end); NumCondIdx = numel(CondIdx);
                for cItemIdx = 1 : numExpNameListFiles

                    if ~isempty(strfind(ExpNameList{cItemIdx, 1}, suffixStrs{cCell}))
                        cItem2sort_ExpName = ExpNameList{cItemIdx, 1};
                        partnameExpr1 = regexp(cItem2sort_ExpName, dExpr1, 'match');
                        if isempty(partnameExpr1); partnameExpr1 = '1'; end

                        if ~isempty(partnameExpr1)
                            cPartNameNum = partnamenumList{cItemIdx, 1};
                            SortIdx = (LastCondIdx-NumCondIdx)+cPartNameNum;
                            ExpNameList_sorted{SortIdx, 1} = ExpNameList{cItemIdx, 1};
                            partnamenumList_sorted{SortIdx, 1} = partnamenumList{cItemIdx, 1};
                        end
                    end
                end
            end
        end
        ExpNameList = ExpNameList_sorted;
        partnamenumList = partnamenumList_sorted;
        save(strcat(CellName, '_FileNameLists.mat'),'ExpNameList','partnamenumList',...
            'numExpNameListFiles','NumStrFilesList')

    end
    strlists(cDirIdx).ExpNameList = ExpNameList;
    strlists(cDirIdx).partnamenumList = partnamenumList;
    strlists(cDirIdx).numExpNameListFiles = numExpNameListFiles;
    strlists(cDirIdx).NumStrFilesList = NumStrFilesList;
    
%CG: count number of data chunks to analyze across all directories 
    NumFilesToThres = size(ExpNameList,1);
    if isnan(ChunkLimit)
        TotalChunkNum = 0;
        ExpNameList = strlists(cDirIdx).ExpNameList;
        for cFileIdx = 1 : NumFilesToThres
            FileTif = ExpNameList{cFileIdx,1}; 
            InfoImage = imfinfo(FileTif);
            NumberImages = length(InfoImage); NumChunks = round(NumberImages/500);
            TotalChunkNum = TotalChunkNum + NumChunks;
            dispstat(strcat('Counting chunks: "', num2str(TotalChunkNum*NumThresholds),'"'))
            if ~isempty(statusStr)
                statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                    '"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "',CellName,'"');
                statusStr = strcat(statusStr, '\nCounting chunks (x no. of thresholds): "', num2str(TotalChunkNum*NumThresholds),'"');
                pause(1); handles.text_status.String = sprintf(statusStr);

            end
        end
        ChunkLimit = TotalChunkNum;
    end
%%%%%%%%testing only    
%     TotalChunkNum = 48; ChunkLimit = 48;
%%%%%%%%testing only    

    if ~isnan(ChunkLimit); TotalChunkNum = ChunkLimit*NumCombos; end
    
    dispstat(strcat('Total chunks: "', num2str(TotalChunkNum),'"'))
    if ~isempty(statusStr)
        statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                '"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "', CellName,'"');
       
        statusStr = strcat(statusStr, '. Counting chunks: "', num2str(TotalChunkNum*NumThresholds),'"');
        pause(1); handles.text_status.String = sprintf(statusStr);

    end
    ChunkTracker = 0;
%CG: collect names of data folder directories is ThresholdDataFolderList
%and pass to EventDurations function.

    for cthreshold = 1 : NumThresholds
        
        E_valThreshold = E_valThresholdList(cthreshold);
        AmplThreshold = Ampl_ThresholdList(cthreshold);
        ChunkLimitTracker = 1;
        ExpNameList = strlists(cDirIdx).ExpNameList;
        CurrentDir = AllDirs{1}; cd(CurrentDir); shiftkey = []; Ref = [];
        DE_in = struct;
        DE_in.ExpNameList = ExpNameList;
        DE_in.E_valThresholdList = E_valThresholdList;
        DE_in.E_valThreshold = E_valThreshold;
        DE_in.AmplThreshold = AmplThreshold;
        DE_in.TotalChunkNum = TotalChunkNum;
        DE_in.NumThresholds = NumThresholds;
        DE_in.ChunkLimit = ChunkLimit;
        DE_in.CurrentDir = CurrentDir;
        DE_in.CorrectShiftSwt = CorrectShiftSwt;
        DE_in.shiftkey = shiftkey;
        DE_in.Ref = Ref;
        DE_in.statusStr = statusStr;
        DE_in.handles = handles;

        startFile = 1;
        for cFileIdx = startFile : NumFilesToThres
            if ChunkLimitTracker <= ChunkLimit
                CurrentFile = ExpNameList{cFileIdx, 1};
                if cFileIdx == startFile
                    ff = imread(CurrentFile,1); dims = size(ff);
                end
                if cFileIdx == 1 || cFileIdx == 4
                    mip = [];
                end
                tifidces = strfind(CurrentFile,'.tif');
                coreFn = CurrentFile(1:tifidces(1)-1);
                if ~isempty(CurrentFile)
                    if NumThresholds == 1
                        if ~isempty(strfind(computer, 'WIN'))
                            if ~exist(strcat(AllDirs{1}, '\', coreFn, ThresholdDataFolderStr),'dir')
                                ThresDataFolder = strcat(AllDirs{1}, '\',coreFn, ThresholdDataFolderStr);
                                mkdir(ThresDataFolder);
                                addpath(ThresDataFolder);
                            else
                                ThresDataFolder = strcat(AllDirs{1}, '\', coreFn, ThresholdDataFolderStr);
                                addpath(ThresDataFolder);
                            end
                        elseif ~isempty(strfind(computer, 'MAC'))
                            if ~exist(strcat(AllDirs{1}, '/', coreFn, ThresholdDataFolderStr),'dir')
                                ThresDataFolder = strcat(AllDirs{1}, '/',coreFn, ThresholdDataFolderStr);
                                mkdir(ThresDataFolder);
                                addpath(ThresDataFolder);
                            else
                                ThresDataFolder = strcat(AllDirs{1}, '/', coreFn, ThresholdDataFolderStr);
                                addpath(ThresDataFolder);
                            end
                        end
                        
                        ThresholdDataFolderList{cFileIdx,cDirIdx} = ThresDataFolder;
                    elseif NumThresholds > 1
                        if ~isempty(strfind(computer, 'WIN'))
                            if ~exist(strcat(AllDirs{1}, '\', coreFn, '_ThresholdOptions_', ThresOptsID),'dir')
                                ThresDataFolder = strcat(AllDirs{1}, '\',coreFn, '_ThresholdOptions_', ThresOptsID);
                                mkdir(ThresDataFolder);
                                addpath(ThresDataFolder);
                            else
                                ThresDataFolder = strcat(AllDirs{1}, '\', coreFn,'_ThresholdOptions_', ThresOptsID);
                                addpath(ThresDataFolder);
                            end
                        elseif ~isempty(strfind(computer, 'MAC'))
                            if ~exist(strcat(AllDirs{1}, '/', coreFn, '_ThresholdOptions_', ThresOptsID),'dir')
                                ThresDataFolder = strcat(AllDirs{1}, '/',coreFn, '_ThresholdOptions_', ThresOptsID);
                                mkdir(ThresDataFolder);
                                addpath(ThresDataFolder);
                            else
                                ThresDataFolder = strcat(AllDirs{1}, '/', coreFn,'_ThresholdOptions_', ThresOptsID);
                                addpath(ThresDataFolder);
                            end
                        end
                    end
                    if strcmp(CorrectShiftSwt, 'On')
                        dispstat(strcat('"',num2str(ChunkTracker), '" down..."', num2str(TotalChunkNum), '" to go! (shift correction being applied)'))
                        if ~isempty(statusStr)
                            if numel(E_valThresholdList) == 1
                                statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                                    '"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "',CellName,'"');
                            else
                                statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold(1)),...
                                    '-',num2str(E_valThreshold(end)),'"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "',...
                                    CellName,'"');
                            end
                            statusStr = strcat(statusStr, '\n"',num2str(ChunkTracker), '" down..."',...
                                num2str(TotalChunkNum), '" to go! (shift correction being applied)');
                            pause(1); handles.text_status.String = sprintf(statusStr);

                        end
                    else
                        dispstat(strcat('"',num2str(ChunkTracker), '" down..."', num2str(TotalChunkNum), '" to go!'))
                        if ~isempty(statusStr)
                            if numel(E_valThresholdList) == 1
                                statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                                    '"; Amplitude threshold =  "', num2str(AmplThreshold),'" for "',...
                                    CellName,'"');
                            else
                                statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold(1)),...
                                    '-',num2str(E_valThreshold(end)),'"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "',...
                                    CellName,'"');
                            end
                            statusStr = strcat(statusStr, '\n"',num2str(ChunkTracker), '" down..."',...
                                num2str(TotalChunkNum), '" to go!');
                            pause(1); handles.text_status.String = sprintf(statusStr);

                        end
                    end

                    DE_in.CurrentFile = CurrentFile;
                    DE_in.dims = dims;
                    DE_in.coreFn = coreFn;
                    DE_in.ChunkTracker = ChunkTracker;
                    DE_in.ThresDataFolder = ThresDataFolder;
                    DE_in.ChunkLimitTracker = ChunkLimitTracker;
                    DE_in.CellName = CellName;
                    DE_in.E_valThreshold = E_valThreshold;
                    DE_in.AmplThreshold = AmplThreshold;
                    DE_in.suffixStrs = suffixStrs;
                    DE_in.mip = mip;

                    [handles, ChunkTracker,ChunkLimitTracker,mip] = ApplyCrossCorr(DE_in);
                end
            end
        end
    end
else
    disp('AARG prep cancelled!')
end  


function [handles, ChunkTracker, ChunkLimitTracker,mip] = ApplyCrossCorr(DE_in)

CurrentFile = DE_in.CurrentFile;
dims = DE_in.dims;
ExpNameList = DE_in.ExpNameList;
coreFn = DE_in.coreFn;
E_valThresholdList = DE_in.E_valThresholdList;
E_valThreshold = DE_in.E_valThreshold;
AmplThreshold = DE_in.AmplThreshold;
TotalChunkNum = DE_in.TotalChunkNum;
ChunkTracker = DE_in.ChunkTracker;
ThresDataFolder = DE_in.ThresDataFolder;
ChunkLimitTracker = DE_in.ChunkLimitTracker;
NumThresholds = DE_in.NumThresholds;
ChunkLimit = DE_in.ChunkLimit;
CurrentDir = DE_in.CurrentDir;
CorrectShiftSwt = DE_in.CorrectShiftSwt;
shiftkey = DE_in.shiftkey;
Ref = DE_in.Ref;
statusStr = DE_in.statusStr;
CellName = DE_in.CellName;
handles = DE_in.handles;
suffixStrs = DE_in.suffixStrs;
mip = DE_in.mip;

ROI_sl = str2double(handles.edit_ROIsl.String);

cd(CurrentDir)
InfoImage = imfinfo(CurrentFile); 
NumberImages = length(InfoImage);
    
if NumThresholds > 1
    NumChunks = ChunkLimit;
elseif NumThresholds == 1
    NumChunks = round(NumberImages/500);
end


for cChunk = 1 : NumChunks

    filestr = strcat(coreFn,'_Chunk_', num2str(cChunk),'_E_val_',num2str(E_valThreshold),'_Ampl_',num2str(AmplThreshold),'Thresholds.mat');
    try
        cd(CurrentDir)
        load(filestr, 'suprathres_Evals','fmedian','E_valThreshold', 'AmplThreshold','ChunkSize','pixelBinningFactor');
        filestr_loc = which(filestr);
        if ~isempty(strfind(computer, 'WIN'))
            if ~strcmp(filestr_loc,strcat(ThresDataFolder,'\',filestr))
%CG: if the current file is not in the ThresDataFolder directory, then it
%is copied to this location. 
                Source = filestr_loc;
                copyfile(Source, ThresDataFolder);
            end
        elseif ~isempty(strfind(computer, 'MAC'))
            if ~strcmp(filestr_loc,strcat(ThresDataFolder,'/',filestr))
%CG: if the current file is not in the ThresDataFolder directory, then it
%is copied to this location. 
                Source = filestr_loc;
                copyfile(Source, ThresDataFolder);
            end
        end
    catch
       
        cd(CurrentDir)
        ChunkStarts = linspace(1, NumberImages, round(NumberImages/500) + 1);
        StartChunkStep = round(ChunkStarts(cChunk));
        EndChunkStep = round(ChunkStarts(cChunk+1)-1);
        tifStack = zeros(dims(1),dims(2),(EndChunkStep-StartChunkStep)+1);
        ii = 0;
        for jj = StartChunkStep : EndChunkStep
            ii = ii + 1;
            tifStack(:,:,ii) = imread(CurrentFile,jj);
        end
     
        tifStack=double(tifStack);
%*******************testing
%         if isempty(mip)
%             mip = max(tifStack, [], 3);
%         else
%             mip = cat(3,mip,max(tifStack, [], 3));
%         end
%        
%         if ChunkTracker == 23
%             mip = max(mip, [], 3);
%             save(strcat('maxIntProj_', suffixStrs{1}, '.mat'), 'mip')
%         elseif ChunkTracker == 47
%             mip = max(mip, [], 3);
%             save(strcat('maxIntProj_', suffixStrs{2}, '.mat'), 'mip')
%         end
%*******************testing
        
        if strcmp(CorrectShiftSwt,'On') &&...
                strcmp(ExpNameList{1,1}, CurrentFile) &&...
                    StartChunkStep == 1
%CG: user has applied the LatShiftCheck function and found the lateral shift
%to be unacceptably high. ApplyCrossCorr will call myShiftCorrect2D to
%correct lateral shift in the raw data file before applying correlation
%thresholds to detect events.
            [~,CellName,~] = fileparts(CurrentDir);
            sk = load(strcat(CellName,'_shiftkey.mat'));
            shiftkey = sk.shiftkey;
            Ref = tifStack(round(shiftkey(2):shiftkey(2)+shiftkey(4)),round(shiftkey(1):shiftkey(1)+shiftkey(3)),1:3);
            Ref = mean(Ref-mySmooth2D_All3_fft(Ref,7,7),3);
        end
        if strcmp(CorrectShiftSwt,'On')
            In = mySmooth2D_All3_fft(tifStack(round(shiftkey(2):shiftkey(2)+shiftkey(4)),round(shiftkey(1):shiftkey(1)+shiftkey(3)),:),1,1);
            In = In-mySmooth2D_All3_fft(In,7,7);
            Out = myShiftDetect2D(In, Ref);
            tifStack = myShiftCorrect2D(tifStack, Out, 'fft');
        end
        
        fmedian = median(tifStack,3); fmedian=myReSize3_bin(fmedian,str2double(handles.edit_pixelBinningFactor.String));
        
        dd=myReSize3_bin(tifStack,str2double(handles.edit_pixelBinningFactor.String));
        clear tifStack
        dd=myBleachingL3(dd);
        dd=mySmooth2D_All3_fft(dd,ROI_sl,ROI_sl);
        AA=myGauss1D_Allx(dd);
        ddN = myNorm_G(dd, AA, [-1 0 0]);
        ddS = myRemove_BK(ddN);

        mm = myNB3(ddS);
        clear ddN
        [Y, X]=hist(mm(:),0:0.001:1); 
        [jj, jm]=max(Y); tt=-20:20; jm=jm+tt; tt(tt<1)=[];tt(tt>numel(jj))=[];
        jm_l = jm > 0; jm = jm(jm_l); A=myGauss1D_fit2(X(jm),Y(jm)); th=A(1)+3*A(2);

        ll=[7 30]; tt=[170 75 45 30 20 14 10 7 5 3 1.7 0.7]; tmp = myNorm_mv(myEXP_Curves(ll, tt));
        th = [th(1) 0];
        ii=find(mm>th(1)&mean(ddS,3)>th(2)); 
        ss=size(ddS); 
        ddSx=reshape(ddS,[],ss(end)); 
        ddSx=ddSx(ii,:);
        ddSx(ddSx<0)=0;
        kk=ll(1)+ll(2); cc=myVar3(ddSx,kk); cc0=cc(:,[kk:end 1:(kk-1)]);
        for jj=1:numel(tt); %fprintf('Now calculate %0.1f\n',tt(jj));    
            cv=myCros_00_d3(ddSx,tmp(jj,:)); cc=cc0; cc(cv<0)=0; ii=find(cc>0); cc(ii)=cv(ii)./sqrt(cc(ii)).*cv(ii)/kk;
            if jj==1; cc1=cc;cc2=cc*0+tt(jj); else ii=find(cc1<cc); cc1(ii)=cc(ii);cc2(ii)=tt(jj); end    
        end; cc=cc1; clear('cv','ii','cc1','cc0');
        cc(:,(1:kk)+end-kk)=0;cc2(:,(1:kk)+end-kk)=0;
        cc=cc(:,[(1:ll(1))+end-ll(1),(1:(end-ll(1)))]);cc2=cc2(:,[(1:ll(1))+end-ll(1),(1:(end-ll(1)))]);
        ii=find(mm>th(1)&mean(ddS,3)>th(2));cc1=zeros(ss(1)*ss(2),ss(3));cc1(ii,:)=cc; cc=reshape(cc1,ss);
        cc1=zeros(ss(1)*ss(2),ss(3));cc1(ii,:)=cc2; cc2=reshape(cc1,ss);
        clear('cc1','ii','jj','ddSx','kk');

        if ~isempty(strfind(computer, 'WIN'))
            md=myAVG3(ddS,4)>AmplThreshold; 
        elseif ~isempty(strfind(computer, 'MAC'))
            md=myAvg3(ddS,4)>AmplThreshold; 
        end
        m1 = cc >= E_valThreshold;
        E_val_array = cc*0;
        E_val_array(m1&md) = 1;

        clear ('cc', 'ntifStackChunk');

        Sz_E_val_array= size(E_val_array);
        ChunkSize = Sz_E_val_array(3);

        conncomps = bwconncomp(E_val_array);
        suprathres_Evals = conncomps.PixelIdxList;
        cd(ThresDataFolder); pixelBinningFactor = num2str(handles.edit_pixelBinningFactor.String);
        save(strcat(coreFn,'_Chunk_', num2str(cChunk),'_E_val_',num2str(E_valThreshold),'_Ampl_',num2str(AmplThreshold),'Thresholds.mat'),...
            'suprathres_Evals','fmedian','E_valThreshold', 'AmplThreshold','ChunkSize','CorrectShiftSwt',...
            'pixelBinningFactor');
%         save(strcat(coreFn,'_Chunk_', num2str(cChunk),'_E_val_',num2str(E_valThreshold),'_Ampl_',num2str(AmplThreshold),'Thresholds.mat'),...
%             'CorrectShiftSwt','-append');

        clear E_val_array;
    end
    ChunkTracker = ChunkTracker + 1;
    ChunkLimitTracker = ChunkLimitTracker + 1;
    dispstat(strcat('"',num2str(ChunkTracker), '" down..."', num2str(TotalChunkNum), '" to go!'))
    if ~isempty(statusStr)
        if numel(E_valThresholdList) == 1
            statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold),...
                '"; Amplitude threshold =  "', num2str(AmplThreshold), '" for "',CellName,'"');
        else
            statusStr = strcat('Detecting events...E value threshold =  "', num2str(E_valThreshold(1)),...
                '-',num2str(E_valThreshold(end)),'"; Amplitude threshold =  "', num2str(AmplThreshold),'" for "',...
                CellName,'"');
        end
        statusStr = strcat(statusStr, '\n"',num2str(ChunkTracker), '" down..."',...
            num2str(TotalChunkNum), '" to go!');
        pause(1); handles.text_status.String = sprintf(statusStr);

    end
end











    