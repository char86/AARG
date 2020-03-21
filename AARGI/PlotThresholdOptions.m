function PlotThresholdOptions(PTOin)

% Outline 
% PlotThresholdOptions handles the plotting of detection with a background
% of the current cell for each of the threshold combinations set by the
% user. 10 results are plotted per figure. 

% Author: Charlie J Gilbride
% version: 1.0.0

%CG: When prompted, the user should select the experiment folders
%containing the threshold option folders having the string ThresOptsID in
%their name. 

if nargin < 1
    options = 'write';
end

ThresOptsID = PTOin.ThresOptsID;
ThresSet = PTOin.ThresSet;
cDir = PTOin.cDir;
specifiedConditionIdx = PTOin.specifiedConditionIdx;
suffixStrs = PTOin.suffixStrs;
options = PTOin.options;
% suffixStrs = {'a_baseline'};

% ThresOptsID = 'E_val_5_Ampl3.5-4';

DisplayBy = 'Threshold'; 
% DisplayBy = 'Condition';

%CG: PTO will display upto a maximum of 10 threshold combinations per
%figure. If DisplayBy is set to 'Threshold', PTO will preferentially put
%different thresholds for the same cell together and always store the results for
%different experiments/cells in different figures. DisplayBy = 'Condition'
%will preferential put images together that belong to the same
%experiment/cell, but represent different treatments/conditions.

scrsz=get(0,'ScreenSize');
% folders = uipickfiles('prompt', 'Select data folders with data to plot');
folders = {};
if nargin < 1
    folders{1}='/Users/cjgilbride/Documents/DecJan2016_17/MATLAB/AARG_pub_datatest/CG0706173';
else
    folders{1} = cDir;
end

Sz_suffixStrs = size(suffixStrs);
TotalConditions = Sz_suffixStrs(1)*Sz_suffixStrs(2);

EvalArray = {}; AmplArray = {}; CountCheck = 0;

dispstat('', 'init')

if iscell(folders)
    
    StrLen = length(suffixStrs{specifiedConditionIdx});
    
    NumFolders = size(folders, 2); 
    
    CondiNameList = {};
    CellNameList = {};
        
    for cCell = 1 : TotalConditions
        for cIdx = 1 : NumFolders
            
            tla = zeros(1,2); cDir = folders{cIdx}; cd(cDir);
            [~, CellName, ~] = fileparts(cDir); 
            if ~isempty(strfind(computer, 'WIN'))
                Slashes = strfind(folders{cIdx}, '\');
            elseif ~isempty(strfind(computer, 'MAC'))
                Slashes = strfind(folders{cIdx}, '/');
            end
            ParentDir = cDir(1:Slashes(end)-1); IDString = CellName(1:2);

            numAnalysisFiles = 0;
            EvalArray = zeros(1, NumFolders, TotalConditions);
            AmplArray = zeros(1, NumFolders, TotalConditions);

            dirInfo = dir;                                                                
            numItems = size(dirInfo,1);                                                 
            for cItem = 1 : numItems                                                         
                currItemName = dirInfo(cItem).name;
                Sz_currItemName = size(currItemName);
                if dirInfo(cItem).isdir 
                    if Sz_currItemName(2) > StrLen
                        if (~isempty(strfind(currItemName, suffixStrs{cCell})) &&...
                                ~isempty(strfind(currItemName, strcat('ThresholdOptions_', ThresOptsID))))||...
                                    (~isempty(strfind(currItemName, suffixStrs{cCell})) &&...
                                        ~isempty(strfind(currItemName, 'ThresholdOptions')) && isempty(ThresOptsID))

                            cd(strcat(cDir, '/', currItemName))

                            dirInfo_sub = dir;                                                                
                            numFiles = size(dirInfo_sub,1);                                                 

                            for cFile = 1 : numFiles
                                cFileName = dirInfo_sub(cFile).name;
                                Sz_cFileName = size(cFileName);
                                if ~isempty(strfind(cFileName, suffixStrs{cCell}))
                                    IDStart = strfind(cFileName, IDString);
                                    IDFinish = strfind(cFileName, suffixStrs{cCell});

                                    FN_length = numel(cFileName(IDStart:IDFinish));
                                else
                                    FN_length = numel(cFileName);
                                end
                                if Sz_cFileName(2) > FN_length 
                                    if strcmp(cFileName(end-2:end),'mat') &&...
                                            ~isempty(strfind(cFileName,'Threshold')) &&...
                                                ~isempty(strfind(cFileName, 'Chunk')) 

                                        CountCheck = CountCheck + 1;
                                        E_Thres = '(?<=\w*_E_val_)\d*(?=_Ampl_)';
                                        E_Val = regexp(cFileName, E_Thres, 'match');
                                        if isempty(E_Val)
%CG: if there is a decimal point in the threshold value for Eval, then the
%expression E_Thres needs to be changed in order to detect it. 
                                            E_Thres = '(?<=\w*_E_val_)\d.\d*(?=_Ampl_)';
                                            E_Val = regexp(cFileName, E_Thres, 'match');
                                        end
                                        AmplThres = '(?<=\w*_Ampl_)\d*(?=Thr)';
                                        AmplVal = regexp(cFileName, AmplThres, 'match');
                                        if isempty(AmplVal)
                                            AmplThres = '(?<=\w*_Ampl_)\d.\d*(?=Thr)';
                                            AmplVal = regexp(cFileName, AmplThres, 'match');
                                        end

                                        if isempty(find(tla,1)) 
                                            tla(1,1) = str2double(E_Val);
                                            tla(1,2) = str2double(AmplVal);
                                        else
                                            tla(end+1, 1) = str2double(E_Val);
                                            Sz_tla = size(tla);
                                            tla(Sz_tla(1), 2) = str2double(AmplVal);
                                        end

                                        lgc_tla_ActThres = tla == str2double(E_Val);
                                        lgc_tla_SizeThres = tla == str2double(AmplVal);

                                        SumLogic = lgc_tla_ActThres + lgc_tla_SizeThres;
                                        CrossSum = sum(SumLogic, 2);

                                        if numel(find(CrossSum == 2)) || numel(find(CrossSum == 4))
%CG: the condition statement exists in order to prevent more chunks of data
%from any one condition being plotted. numel(find(CrossSum == 4))
%accommodates situations where E_Val = AmplVal.
                                            numAnalysisFiles = numAnalysisFiles + 1;
                                            CellNameList{1, cIdx, cCell} = CellName;
                                            if strcmp(DisplayBy, 'Threshold') 
                                                CondiNameList(numAnalysisFiles, cIdx, cCell) = strcat('E_val: ',E_Val, '; ', 'Ampl: ', AmplVal);
                                            elseif strcmp(DisplayBy, 'Condition')
                                                CondiNameList(numAnalysisFiles, cIdx, cCell) = strcat('; E_val: ',E_Val, '; ', 'Ampl: ', AmplVal);
                                            end

                                            CondiFiles{numAnalysisFiles, cIdx, cCell} = cFileName; 
                                            if numAnalysisFiles == 1
                                                EvalArray(1, cIdx, cCell) = str2double(E_Val{1});
                                                AmplArray(1, cIdx, cCell) = str2double(AmplVal{1});
                                            elseif numAnalysisFiles > 1
                                                EvalArray(end+1, cIdx, cCell) = str2double(E_Val{1});
                                                AmplArray(end+1, cIdx, cCell) = str2double(AmplVal{1});
                                            end
                                        end
                                    end
                                end
                            end
                        end

                    end
                end   
            end
        end
    end
    
    if strcmp(options, 'write') || strcmp(options, 'readonly')
    
%CG: CondiNameList can be empty if ~isempty(specifiedConditionIdx) 
   
        CondiNameListEmptiesOut = CondiNameList(~cellfun('isempty',CondiNameList));
        uCondiNameListEmptiesOut = unique(CondiNameListEmptiesOut);
        Sz_uCondiNameListEmptiesOut = size(uCondiNameListEmptiesOut);
        uCondis = Sz_uCondiNameListEmptiesOut(1)*Sz_uCondiNameListEmptiesOut(2);
        Sz_CondiNameList = size(CondiNameList);

        if strcmp(DisplayBy, 'Threshold')
    %CG: the max number of options that can comfortably fit onto one figure is
    %10. So, if there are more than 10 options, we need to create more figures.
            NumOfChunks = (Sz_CondiNameList(1)/uCondis);
            FileNameSelect = 1:NumOfChunks:numAnalysisFiles;
        elseif strcmp(DisplayBy, 'Condition')
            NumOfChunks = (Sz_CondiNameList(1)/uCondis);
            FileNameSelect = 1:NumOfChunks:numAnalysisFiles;
        end

        if ~isempty(EvalArray)
            lgcEvalArray = EvalArray == 0; EvalArray(lgcEvalArray) = [];
            lgcAmplArray = AmplArray == 0; AmplArray(lgcAmplArray) = [];
    %CG: When DisplayBy = 'Threshold', EvalArray and AmplArray will have null element 
    %values. These need to be removed before calling the sort function,
    %otherwise sort will order the null elements before the non-null elements
    %and this will violate assumptions made later. 
            [~, ActArray_sortIdx] = sort(EvalArray, 'ascend');
            [~, SzArray_sortIdx] = sort(AmplArray, 'ascend');
            CondiNameList_sorted = {};CondiFiles_sorted = {};
            for cCell = 1 : TotalConditions
                for cIdx = 1 : NumFolders
                    for cThres = 1 : numAnalysisFiles

                        if strcmp(ThresSet, 'E_val')
                            IdxSet = ActArray_sortIdx(cThres); 
                        elseif strcmp(ThresSet, 'Ampl')
                            IdxSet = SzArray_sortIdx(cThres);   
                        end
                        if ~isempty(CondiNameList{IdxSet, cIdx, cCell})
                            CondiNameList_sorted{cThres, cIdx, cCell} = CondiNameList{IdxSet, cIdx, cCell};
                            CondiFiles_sorted{cThres, cIdx, cCell} = CondiFiles{IdxSet, cIdx, cCell};
                        end

                    end
                end
            end                

            SumOfAll = TotalConditions*NumFolders*Sz_CondiNameList(1);
            wbStep = floor(SumOfAll/100);
            wbCounter = 0;
            if strcmp(DisplayBy, 'Threshold')

                for cIdx = 1 : NumFolders

                    quitWhile = 0; cCell = 0;
                    while cCell <= TotalConditions && quitWhile == 0

                        cCell = cCell + 1;
                        if ~isempty(specifiedConditionIdx)
                            cCell = specifiedConditionIdx; quitWhile = 1;
                        end

                        CellName = CellNameList{1, cIdx, cCell}; handleArray = {};
                        handleArray{1} = figure('OuterPosition',[40,40,scrsz(3)-80,scrsz(4)-80],'Name',... % create figures
                            strcat(strcat(CellName, suffixStrs{cCell}),': Threshold Options'),'NumberTitle','off','Visible','on');

                        CondiCount = 0;
                        for condistep = 1 : uCondis

                            if mod(condistep-1,10) == 0 && condistep-1 ~= 0
                                handleArray{end+1}=figure('OuterPosition',[40,40,scrsz(3)-80,scrsz(4)-80],'Name',... % create figures
                                    strcat(strcat(CellName, suffixStrs{cCell}),': Threshold Options'),'NumberTitle','off','Visible','on');
                                CondiCount = 0;
                            end
                            CondiCount = CondiCount + 1;
                            [~, CellName, ~] = fileparts(folders{cIdx});
                            IDStrSta = strfind(CellName, IDString);

                            if IDStrSta > 1
                                CellName = CellName(IDStrSta:end);
                            end
                            cd(strcat(folders{cIdx}, '/', CellName, suffixStrs{cCell}, '_ThresholdOptions_', ThresOptsID))

                            for cChunk = 1 : NumOfChunks

                                wbCounter = wbCounter + 1;

                                if cChunk == 1
                                    condiName = CondiNameList_sorted{FileNameSelect(condistep) + (cChunk-1), cIdx, cCell};
                                    results = load(CondiFiles_sorted{FileNameSelect(condistep) + (cChunk-1), cIdx, cCell});
                                    fmedian = results.fmedian;
                                    suprathres_Evals=results.suprathres_Evals;
                                    numEvents = size(suprathres_Evals,2);
                                    ChunkSize = results.ChunkSize;
                                else
                                    results = load(CondiFiles_sorted{FileNameSelect(condistep) + (cChunk-1), cIdx, cCell});
                                    fmedian_par = results.fmedian;
                                    suprathres_Evals_par=results.suprathres_Evals;
                                    numEvents_par = size(suprathres_Evals_par,2);
    %                                 fmedian = mean(3,cat(3, fmedian, fmedian_par),'omitnan');
                                    suprathres_Evals = cat(2,suprathres_Evals,suprathres_Evals_par);
                                    numEvents = numEvents + numEvents_par;
                                    if ChunkSize < results.ChunkSize
                                        ChunkSize = results.ChunkSize;
                                    end
                                end
                                if wbStep ~= 0
                                    if ~mod(wbCounter, wbStep) || wbCounter == SumOfAll
                                        dispstat(strcat('~', num2str((wbCounter/SumOfAll)*100), '% complete'))
                                    end
                                else
                                    dispstat(strcat('~', num2str((wbCounter/SumOfAll)*100), '% complete'))
                                end
                            end
                            ss=size(suprathres_Evals,2);
                            ActOverview=false(size(fmedian));
                            for cComp = 1:ss
                                [rr,cc,~] = ind2sub([size(fmedian,1), size(fmedian,2), ChunkSize],suprathres_Evals{cComp});
                                eventidx = sub2ind([size(fmedian)], rr, cc);
                                ActOverview(eventidx) = 1;
                            end

                            ActOverview = logical(ActOverview);
                            logFluoImage=log10(fmedian);                                                
                            logFluoImage=logFluoImage-min(logFluoImage(:));                             
                            logFluoImage=logFluoImage/max(logFluoImage(:));                             
                            logFluoImage=uint8(logFluoImage*256);                                       
                            cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);        
                            colormap(cyanColorMap);                                                     
                            fmedianRGB=ind2rgb(logFluoImage,cyanColorMap);

                            figure(handleArray{end});
                            if condistep == 1 || (mod(condistep-1,10) == 0 && condistep-1 ~= 0); hold on; end

                            leftP = mod(CondiCount-1,5)/5;
                            bottomP = (ceil(CondiCount/5)-1)/2;

                            subplot('Position',[leftP, bottomP, 1/5, 1/2]);
                            [RowClr, ColClr, ~] = ind2sub(size(fmedianRGB), find(ActOverview));

                            FirstDim = ones(1, numel(RowClr));
                            Idces_1D = sub2ind(size(fmedianRGB),RowClr, ColClr, FirstDim');
    %                         Idces_1D = ActOverview == 1;
                            SecDim = ones(1, numel(RowClr)).*2;
                            Idces_2D = sub2ind(size(fmedianRGB),RowClr, ColClr, SecDim');
                            ThirdDim = ones(1, numel(RowClr)).*3;
                            Idces_3D = sub2ind(size(fmedianRGB),RowClr, ColClr, ThirdDim');

                            fmedianRGB(Idces_1D) = 1;
                            fmedianRGB(Idces_2D) = 1;
                            fmedianRGB(Idces_3D) = 1;
                            imagesc(fmedianRGB);axis square; 

                            title(strcat(condiName,' Number of Events:',num2str(numEvents))); 

                            set(gca, 'XTickLabel', '', 'TickLength', [0 0]);
                            set(gca, 'YTickLabel', '', 'TickLength', [0 0]);
                            set(gca, 'box', 'off')

                        end

                    end
                    cd(ParentDir)
                    mkdir(strcat('ThresholdRecords_', CellName))
                    cd(strcat(ParentDir, '/', strcat('ThresholdRecords_', CellName)))

                    dirInfo = dir;                                                                
                    numItems = size(dirInfo,1);

                    cItemV = 0;
                    Version = [];
                    while cItemV < numItems
                        cItemV = cItemV + 1;
                        ci = dirInfo(cItemV).name;
                        if numel(dirInfo(cItemV).name) >= numel(ThresOptsID)
                            if ~isempty(strfind(dirInfo(cItemV).name, ThresOptsID))
                                r = regexp(ci, '(?<=_v)\d*', 'match');
                                if isempty(r) && isempty(Version)
                                    Version = 1;
                                elseif ~isempty(r)
                                    Version = str2double(r) + 1;
                                end
                            end
                        end

                    end

                    if ~isempty(Version) && strcmp(options, 'readonly')
                        ThresOptsID = strcat(ThresOptsID, '_v', num2str(Version));
                    end
                    if ~strcmp(options, 'readonly')
                        mkdir(ThresOptsID)
                        cd(strcat(ParentDir, '/', strcat('ThresholdRecords_', CellName), '/', ThresOptsID))
                        saveas(handleArray{end}, strcat(CellName, suffixStrs{cCell}), 'fig')
                        if ~isempty(Version)
                            save('Version.mat', 'Version')
                        end

                        cd(cDir)
                    end
                end
           end

        elseif isempty(EvalArray)
            disp(strcat('An error occurred for the current data folder', ': "', CellName,...
                '". Perhaps the variable "ThresOptsID" is not set correctly'))
        end
        if ~isempty(EvalArray)
            set(handleArray{end},'Visible','on');
        end
    elseif strcmp(options, 'read')
        slashType = PTOin.slashType;
        cd(strcat(ParentDir, slashType, 'ThresholdRecords_', CellName, slashType, ThresOptsID))
        open(strcat(CellName, suffixStrs{specifiedConditionIdx},'.fig'))
    end
    
end
