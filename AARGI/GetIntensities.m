function GetIntensities(handles,varargin)

% Outline
% GetIntensities captures raw fluorescence intensity data within ROIs
% established by AARG_Core. In each frame, the mean fluorescence intensity
% is taken for each ROI so there is a single, average, fluorescence
% intensity value for each ROI in each frame. 

% Author: Charlie J Gilbride
% version: 1.0.0

p = inputParser;
addParameter(p,'suffixStrs',{'a_baseline','b_TTAP2'})
addParameter(p,'partnamesuff','-file00')
addParameter(p,'DataType', 'Raw') % CG: previouly used 'Raw','170' or 'both'.
%CG: DataType controls how intensity data is collected. Strings Raw and 170
%refer to the pixel resolution. 512x512 is the pixel resolution used for
%data acquisition in the initial study. This resolution is reduced during
%the data processing stage to 170x170. When DataType is set to '170' then
%fluorescence intensity will be measured from data that is normalized with
%background subtraction applied. With 'Raw' measurements are applied to the
%raw data files. 
addParameter(p, 'osType', '')
parse(p,varargin{:}); param = p.Results; 

suffixStrs = param.suffixStrs; partnamesuff = param.partnamesuff; DataType = param.DataType;
osType = param.osType;
dispstat('', 'init')

if nargin == 0 
    Dirs = uipickfiles('Prompt', 'Select data folders to collect fluorescence intensity data');
    handles = [];
else
    Dirs = {}; Dirs{1} = handles.cDir; 
end

if iscell(Dirs)
    
    NumOfDirs = size(Dirs,2); CA_wbcount = {}; NumConditions = size(suffixStrs,2);
    
    for cDirIdx = 1 : NumOfDirs
        cDir = Dirs{cDirIdx}; [~, CellName, ~] = fileparts(cDir);        
        FluoIntFolder = strcat(CellName, '_fIntensity');
        
        if ~isempty(strfind(osType, 'MAC')) 
            fIntensity_path = strcat(cDir, '/', FluoIntFolder);        
        else 
            fIntensity_path = strcat(cDir, '\', FluoIntFolder);        
        end
        cd(cDir); 
        lists = load(strcat(CellName, '_FileNameLists.mat')); ExpNameList = lists.ExpNameList;
        TotalPartitionNum = size(ExpNameList,1); NumStrFilesList = lists.NumStrFilesList;
        
        SkipCount = 0; cCondition = []; wbEndPoint = 0; tallycount = 0;
        for cParti = 1 : TotalPartitionNum 
            cItem = ExpNameList{cParti}; tifindex = strfind(cItem, '.tif');
            
            if exist(FluoIntFolder,'dir') == 7 
                cd(fIntensity_path)
                if isempty(cCondition)
                    cCondition = NumStrFilesList(cParti);
                    CoreType = strcat(CellName,suffixStrs{cCondition}, '_IntensityFiles.mat');
                    if exist(CoreType, 'file') == 2
                        SkipCount = 1;
                    else
                        SkipCount = 0;
                    end
                else
                    if NumStrFilesList(cParti) ~= cCondition
                        cCondition = NumStrFilesList(cParti);
                        CoreType = strcat(CellName,suffixStrs{cCondition}, '_IntensityFiles.mat');
                        if exist(CoreType, 'file') == 2
                            SkipCount = 1;
                        else
                            SkipCount = 0;
                        end
                    end
                end
            else
                SkipCount = 0; 
            end 
              
            if SkipCount == 0
                cFolderTD = strcat(cItem(1:tifindex-1), '_ThresholdData');
                cd(strcat(cDir, '/', cFolderTD)) 
                lists  = load(strcat(cFolderTD, '_FileOrder.mat'));
                TotalChunkSize = lists.TotalChunkSize;
                tallycount = tallycount + 1;
                wbEndPoint = wbEndPoint + TotalChunkSize;
            end
        end
        wbEndPoint = wbEndPoint + tallycount;
        wbStep = floor(wbEndPoint/100);
        CA_wbcount{cDirIdx,1} = wbEndPoint; CA_wbcount{cDirIdx,2} = wbStep;
    end
    
    for cDirIdx = 1 : NumOfDirs
        cDir = Dirs{cDirIdx}; [~, CellName, ~] = fileparts(cDir);
        
        FluoIntFolder = strcat(CellName, '_fIntensity');
        fIntensity_path = strcat(cDir, '/', FluoIntFolder);
        
        cd(cDir); 
        lists = load(strcat(CellName, '_FileNameLists.mat')); ExpNameList = lists.ExpNameList;
        NumStrFilesList = lists.NumStrFilesList;
        TotalPartitionNum = size(ExpNameList,1); 
        
        wbEndPoint = CA_wbcount{cDirIdx,1}; wbStep = CA_wbcount{cDirIdx,2};
        
        timestep = 1; ChunkSizeTotal = 0; TDFolders_sorted = {}; foldercount = 0;
        PrimeTDFolders_sorted = {};
        for cParti = 1 : TotalPartitionNum 
            cItem = ExpNameList{cParti}; tifindex = strfind(cItem, '.tif');
            cFolderTD = strcat(cItem(1:tifindex-1), '_ThresholdData');
            cd(strcat(cDir, '/', cFolderTD)) 
            TDFolders_sorted{cParti,1} = strcat(cDir, '/', cFolderTD);
            if isempty(strfind(cFolderTD,partnamesuff))
                foldercount = foldercount + 1;
                PrimeTDFolders_sorted{foldercount,1} = strcat(cDir, '/', cFolderTD);
                cd(strcat(cDir, '/', cFolderTD))
                CAs = load(strcat('CAs_Core_',cItem(1:tifindex-1)));
                PrimeTDFolders_sorted{foldercount,2} = size(CAs.CA_ROIs);
                cd(cDir)
            end
        end
        
        ROI_map = []; ROI_Idx_map_Raw = []; 
        cCondition = []; SkipFluoInt = 0; ChunkSizeTotal = 0;
        for cParti = 1 : TotalPartitionNum
            cItem = ExpNameList{cParti}; cStr = suffixStrs{NumStrFilesList(cParti)};
            if cParti > 1; lastStr = suffixStrs{NumStrFilesList(cParti-1)}; else lastStr = ''; end
            
            if isempty(cCondition)
                cCondition = NumStrFilesList(1);
                
                if cParti == 1
                    cd(PrimeTDFolders_sorted{cCondition,1}); 
                    ROIs = load(strcat('CAs_Core_', CellName, suffixStrs{cCondition}, '.mat')); 
                    CA_ROIs = ROIs.CA_ROIs; cd(cDir)
                    try 
                        addpath(fIntensity_path); cd(fIntensity_path)
                    catch
                        mkdir(FluoIntFolder); addpath(fIntensity_path); 
                    end
                end
                
                CoreType = strcat(CellName,suffixStrs{cCondition}, '_IntensityFiles.mat');
                if exist(CoreType,'file') == 2
                    SkipFluoInt = 1; 
                else
                    SkipFluoInt = 0;
                    CA_IntensityRaw = cell(PrimeTDFolders_sorted{cCondition,2});
                    CA_Intensity170 = cell(PrimeTDFolders_sorted{cCondition,2});
                end
            else
                if NumStrFilesList(cParti) ~= cCondition
                    cd(fIntensity_path); CoreType = strcat(CellName,suffixStrs{cCondition}, '_IntensityFiles.mat');
                    if SkipFluoInt == 0
                        try
                            ints = load(CoreType);
                            if ~isempty(ints)
                                if strcmp(DataType,'Raw') 
                                    save(CoreType, 'CA_IntensityRaw','-append')
                                elseif strcmp(DataType,'170')
                                    save(CoreType, 'CA_Intensity170','-append')
                                elseif strcmp(DataType,'both')
                                    save(CoreType, 'CA_Intensity170','CA_IntensityRaw','-append')
                                end
                            end
                        catch
                            if strcmp(DataType,'Raw') 
                                save(CoreType, 'CA_IntensityRaw')
                            elseif strcmp(DataType,'170')
                                save(CoreType, 'CA_Intensity170')
                            elseif strcmp(DataType,'both')
                                save(CoreType, 'CA_Intensity170','CA_IntensityRaw')
                            end
                        end
                    end
                    cCondition = NumStrFilesList(cParti);
                    cd(PrimeTDFolders_sorted{cCondition,1}); 
                    ROIs = load(strcat('CAs_Core_', CellName, suffixStrs{cCondition}, '.mat')); 
                    CA_ROIs = ROIs.CA_ROIs; cd(cDir)
                    try 
                        addpath(fIntensity_path); cd(fIntensity_path)
                    catch
                        mkdir(FluoIntFolder); addpath(fIntensity_path); 
                    end
                    
                    CoreType = strcat(CellName,suffixStrs{cCondition}, '_IntensityFiles.mat');
                    if exist(CoreType,'file') == 2
                        SkipFluoInt = 1; ChunkSizeTotal = ChunkSizeTotal + (wbEndPoint*0.5); 
                    else
                        SkipFluoInt = 0;
                        CA_IntensityRaw = cell(PrimeTDFolders_sorted{cCondition,2});
                        CA_Intensity170 = cell(PrimeTDFolders_sorted{cCondition,2});
                    end
                    
                elseif NumStrFilesList(cParti) == cCondition
                    
                end
            end
            
            if SkipFluoInt == 0  
                
                tdf = TDFolders_sorted{cParti};
                if isempty(strfind(cItem,partnamesuff))
                    fmedian = []; cfile = 0; cd(tdf); dirInfo = dir;
                    while isempty(fmedian) 
                        cfile = cfile + 1; currentfile = dirInfo(cfile).name;
                        if ~isempty(strfind(currentfile, 'Chunk_1'))
                            data = load(currentfile); fmedian = data.fmedian;
                        end
                    end
                end
                
                cd(cDir)
                xdim512 = size(imread(cItem,1),1); ydim512 = size(imread(cItem,1),1);
                
                Slashes = strfind(tdf,'/');
                fileorder = load(strcat(TDFolders_sorted{cParti}, '/', tdf(Slashes(end)+1:end),'_FileOrder.mat'));

                fii = struct;

                fii.CA_ROIs = CA_ROIs;
                fii.CA_IntensityRaw = CA_IntensityRaw;
                fii.CA_Intensity170 = CA_Intensity170;
                fii.DataType = DataType;
                fii.ExpNameList = ExpNameList;
                fii.ChunkSizeTotal = ChunkSizeTotal;
                fii.fmedian = fmedian;
                fii.cParti = cParti;
                fii.cDir = cDir;
                fii.TotalPartitionNum = TotalPartitionNum;
                fii.fileorder = fileorder;
                fii.xdim512 = xdim512;
                fii.ydim512 = ydim512;
                fii.ROI_map = ROI_map;
                fii.ROI_Idx_map_Raw = ROI_Idx_map_Raw;
                fii.ChunkSizeTotal = ChunkSizeTotal;
                fii.wbEndPoint = wbEndPoint;
                fii.wbStep = wbStep;
                fii.handles = handles;
                fii.CellName = CellName;
                fii.cCondition = cCondition;
                fii.NumConditions = NumConditions;
                fii.suffixStrs = suffixStrs;
                fii.NumStrFilesList = NumStrFilesList;
                fii.fIntensity_path = fIntensity_path;
                fii.cStr = cStr;
                fii.lastStr = lastStr;
                
                [fio] = FluoInt(fii);

                CA_Intensity170 = fio.CA_Intensity170;
                CA_IntensityRaw = fio.CA_IntensityRaw;
                CA_ROIs = fio.CA_ROIs;
                ROI_map = fio.ROI_map;
                ROI_Idx_map_Raw = fio.ROI_Idx_map_Raw;
                ChunkSizeTotal = fio.ChunkSizeTotal;
            end

        end
    end
end


function [fio] = FluoInt(fii)

%CG: FluoInt captures raw fluorescence intensity data from the raw tif
%files if the DataType variable is set to 'raw' (this is the default
%setting).

%CG: FluoInt can also normalize the imaging data (if DataType variable is
%set to '170' or 'both') from the tif files containing the raw data (just
%like in the DetectEvents function). It then collects the average
%fluorescence intensity, frame-by-frame, for each ROI. A new CA_...
%variable is then created, CA_Intensity170, which will store the
%fluoescence intensity data for each ROI in the corresponding cell. At the
%same time, FluoInt rescales ROIs by a factor of 3 (to the original 512*512
%dimensions) and collects fluoescence intensity from the raw data.

CA_ROIs = fii.CA_ROIs;
CA_IntensityRaw = fii.CA_IntensityRaw;
CA_Intensity170 = fii.CA_Intensity170;
DataType = fii.DataType;
ExpNameList = fii.ExpNameList;
fmedian = fii.fmedian;
cParti = fii.cParti;
fileorder = fii.fileorder;
xdim512 = fii.xdim512;
ydim512 = fii.ydim512;
ROI_map = fii.ROI_map;
ROI_Idx_map_Raw = fii.ROI_Idx_map_Raw;
ChunkSizeTotal = fii.ChunkSizeTotal;
wbEndPoint = fii.wbEndPoint;
wbStep = fii.wbStep;
handles = fii.handles;
CellName = fii.CellName;
cCondition = fii.cCondition;
NumConditions = fii.NumConditions;
suffixStrs = fii.suffixStrs;
NumStrFilesList = fii.NumStrFilesList;
fIntensity_path = fii.fIntensity_path;
cStr = fii.cStr;
lastStr = fii.lastStr;

Sz_CA_ROIs = size(CA_ROIs);
xdim170 = size(fmedian,1); ydim170 = size(fmedian,1); bw512 = zeros(xdim170,ydim170);

TotalCCSList = fileorder.TotalCCSList;
PartiMorphList = TotalCCSList(1):TotalCCSList(end);

if isempty(lastStr) || ~strcmp(cStr, lastStr)
%CG: if lastStr is '' then cParti = 1 and we must set parameters like
%ROI_map for the first condition. If ~strcmp(cStr, lastStr), this must
%mean we have started a new condition and must set the parameters for this
%condition.
    ROI_map = zeros(1,3); ROISwitch = 0;
    CA_Intensity170 = cell(Sz_CA_ROIs); CA_IntensityRaw = cell(Sz_CA_ROIs);

    if numel(Sz_CA_ROIs) == 2
        TotalMapSize = Sz_CA_ROIs(1)*Sz_CA_ROIs(2);
    elseif numel(Sz_CA_ROIs) == 3
        TotalMapSize = Sz_CA_ROIs(1)*Sz_CA_ROIs(2)*Sz_CA_ROIs(3);
    end
    for cCellIdx = 1 : TotalMapSize
        if ~isempty(CA_ROIs{cCellIdx}) && ~ischar(CA_ROIs{cCellIdx})
            new_ROI_sl = numel(CA_ROIs{cCellIdx}); ROI_Idx_map_170 = zeros(1,new_ROI_sl);
            ROI_Idx_map_Raw = zeros(1,new_ROI_sl*new_ROI_sl);
        end
    end
    cellcounter = 0;
    for cCellIdx = 1 : TotalMapSize
        
        if ~isempty(CA_ROIs{cCellIdx}) && ROISwitch == 0
            cellcounter = cellcounter + 1;

            [r, c, z] = ind2sub(Sz_CA_ROIs, cCellIdx);
            ROI_Idx_map_170(cellcounter, :) = CA_ROIs{cCellIdx};

            ROI_map(cellcounter,1) = r; ROI_map(cellcounter,2) = c;
            ROI_map(cellcounter,3) = z; ROISwitch = 1;
%CG: find the index values for each ROI based on the data with the original
%size. 
            [rr,cc] = ind2sub([xdim170,ydim170], CA_ROIs{cCellIdx}); rr=rr.*3; cc=cc.*3;
            new_ROI_sl = numel(CA_ROIs{cCellIdx}); ROI_sl = sqrt(numel(CA_ROIs{cCellIdx}));
            rr_new = (rr(1)-1:1:rr(ROI_sl)+1); rr_new = repmat(rr_new,1,new_ROI_sl);
            unique_cc = unique(cc); exp_unique_cc = (unique_cc(1)-1:1:unique_cc(end)+1);
            NumExpUniqueEls = numel(exp_unique_cc); cc_new = [];
            for cEl = 1 : NumExpUniqueEls
                El = exp_unique_cc(cEl); cc_new_ones = ones(1,new_ROI_sl);
                cc_new_ones = cc_new_ones.*El;
                if isempty(cc_new); cc_new = cc_new_ones;
                else cc_new = cat(2,cc_new,cc_new_ones); end
            end
            Idces512 = sub2ind([xdim512, ydim512], rr_new, cc_new); ROI_Idx_map_Raw(cellcounter, :) = Idces512; 
            bw512(rr_new,cc_new)=1;

        elseif ~isempty(CA_ROIs{cCellIdx}) && ROISwitch == 1
            cellcounter = cellcounter + 1;

            [r, c, z] = ind2sub(Sz_CA_ROIs, cCellIdx);
            ROI_Idx_map_170(cellcounter, :) = CA_ROIs{cCellIdx};

            ROI_map(cellcounter,1) = r; ROI_map(cellcounter,2) = c;
            ROI_map(cellcounter,3) = z;

            [rr,cc] = ind2sub([170,170], CA_ROIs{cCellIdx}); rr=rr.*3; cc=cc.*3;
            new_ROI_sl = numel(CA_ROIs{cCellIdx}); ROI_sl = sqrt(numel(CA_ROIs{cCellIdx}));
            rr_new = (rr(1)-1:1:rr(ROI_sl)+1); rr_new = repmat(rr_new,1,new_ROI_sl);
            unique_cc = unique(cc); exp_unique_cc = (unique_cc(1)-1:1:unique_cc(end)+1);
            NumExpUniqueEls = numel(exp_unique_cc); cc_new = [];
            for cEl = 1 : NumExpUniqueEls
                El = exp_unique_cc(cEl); cc_new_ones = ones(1,new_ROI_sl);
                cc_new_ones = cc_new_ones.*El;
                if isempty(cc_new); cc_new = cc_new_ones;
                else cc_new = cat(2,cc_new,cc_new_ones); end
            end
            Idces512 = sub2ind([xdim512, ydim512], rr_new, cc_new); ROI_Idx_map_Raw(cellcounter, :) = Idces512; 
            bw512(rr_new,cc_new)=1;
        end
    end
end
NumDataChunks = numel(TotalCCSList);
for cChunk = 2 : NumDataChunks

    StartChunkStep = TotalCCSList(cChunk-1); 
    if cChunk == NumDataChunks; EndChunkStep = TotalCCSList(cChunk);
    else EndChunkStep = TotalCCSList(cChunk)-1; end
    ChunkSize = TotalCCSList(cChunk) - TotalCCSList(cChunk-1);
    QS512 = zeros(xdim512,ydim512,(EndChunkStep-StartChunkStep)+1);
    morphStartChunkStep = find(PartiMorphList == StartChunkStep);
    morphEndChunkStep = find(PartiMorphList == EndChunkStep);
    ii = 0;
    for jj = morphStartChunkStep : morphEndChunkStep
        ii = ii + 1; QS512(:,:,ii) = imread(ExpNameList{cParti},jj);
    end
    
    if strcmp(DataType,'170') || strcmp(DataType,'both')
        dd=myReSize3_bin(QS512,3); dd=myBleachingL3(dd); dd=mySmooth2D_All3_fft(dd,1.7,1.7);
        AA=myGauss1D_Allx(dd); ddN = myNorm_G(dd, AA, [-1 0 0]); QS170 = myRemove_BK(ddN);
        clear('dd','AA','ddN')
    end

    Sz_ROI_map = size(ROI_map);

    for cRow = 1 : Sz_ROI_map(1)
        if strcmp(DataType,'Raw') 
            [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));
            numberOfPixels = size(rr,2);
            sQS512 = NaN(numberOfPixels,size(QS512,3));
            for cPixel = 1 : numberOfPixels
                sQS512(cPixel, :) = QS512(rr(cPixel),cc(cPixel),:);
            end
            mQS512 = mean(sQS512,1);
            if StartChunkStep == 1
                CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} = reshape(mQS512,1,numel(mQS512));
            else
                CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} =...
                    cat(2,CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)}, reshape(mQS512,1,numel(mQS512)));
            end
        end

%         if StartChunkStep == 1 
%             if strcmp(DataType,'170') 
%                 [rr,cc,~] = ind2sub(size(QS170),ROI_Idx_map_170(cRow,:));sQS170 = QS170(rr,cc,:); mQS170 = mean(mean(sQS170,1),2); 
%                 CA_Intensity170{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} = reshape(mQS170,1,numel(mQS170));
%             elseif strcmp(DataType,'Raw') 
%                 [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));
%                 numberOfPixels = size(rr,2);
%                 sQS512 = NaN(numberOfPixels,size(QS512,3));
%                 for cPixel = 1 : numberOfPixels
%                     sQS512(cPixel, :) = QS512(rr(cPixel),cc(cPixel),:);
%                 end
%                 mQS512 = mean(sQS512,1);
% %                 [rr,cc,~] = ind2sub(size(QS170),ROI_Idx_map_170(cRow,:));sQS170 = QS170(rr,cc,:); mQS170 = mean(mean(sQS170,1),2);
%                 CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} = reshape(mQS512,1,numel(mQS512));
%             elseif strcmp(DataType,'both')
%                 [rr,cc,~] = ind2sub(size(QS170),ROI_Idx_map_170(cRow,:));sQS170 = QS170(rr,cc,:); mQS170 = mean(mean(sQS170,1),2); 
%                 CA_Intensity170{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} = reshape(mQS170,1,numel(mQS170));
%                 
% %                 [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));sQS512 = QS512(rr,cc,:); mQS512 = mean(mean(sQS512,1),2); 
%                 CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} = reshape(mQS512,1,numel(mQS512));
%             end
%         else
%             if strcmp(DataType,'170') 
%                 [rr,cc,~] = ind2sub(size(QS170),ROI_Idx_map_170(cRow,:));sQS170 = QS170(rr,cc,:); mQS170 = mean(mean(sQS170,1),2); 
%                 CA_Intensity170{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} =...
%                     cat(2,CA_Intensity170{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)}, reshape(mQS170,1,numel(mQS512)));
%             elseif strcmp(DataType,'Raw') 
%                 [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));
%                 numberOfPixels = size(rr,2);
%                 sQS512 = NaN(numberOfPixels,size(QS512,3));
%                 for cPixel = 1 : numberOfPixels
%                     sQS512(cPixel, :) = QS512(rr(cPixel),cc(cPixel),:);
%                 end
%                 mQS512 = mean(sQS512,1);
% %                 [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));sQS512 = QS512(rr,cc,:); mQS512 = mean(mean(sQS512,1),2);
%                 CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} =...
%                     cat(2,CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)}, reshape(mQS512,1,numel(mQS512)));
%             elseif strcmp(DataType,'both')
%                 [rr,cc,~] = ind2sub(size(QS170),ROI_Idx_map_170(cRow,:));sQS170 = QS170(rr,cc,:); mQS170 = mean(mean(sQS170,1),2); 
%                 CA_Intensity170{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} =...
%                     cat(2,CA_Intensity170{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)}, reshape(mQS170,1,numel(mQS170)));
%                 [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));
%                 numberOfPixels = size(rr,2);
%                 sQS512 = NaN(numberOfPixels,size(QS512,3));
%                 for cPixel = 1 : numberOfPixels
%                     sQS512(cPixel, :) = QS512(rr(cPixel),cc(cPixel),:);
%                 end
% %                 [rr,cc,~] = ind2sub(size(QS512),ROI_Idx_map_Raw(cRow,:));sQS512 = QS512(rr,cc,:); mQS512 = mean(mean(sQS512,1),2);
%                 CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)} =...
%                     cat(2,CA_IntensityRaw{ROI_map(cRow, 1), ROI_map(cRow, 2), ROI_map(cRow, 3)}, reshape(mQS512,1,numel(mQS512)));
%             end
%         end 
    end

    ChunkSizeTotal = ChunkSizeTotal + ChunkSize;
    
    if ~isempty(handles)
        statusStr = strcat('saving fluorescence intensity files...',...
            num2str((ChunkSizeTotal/wbEndPoint)*100), '% complete for "',CellName,'" at condition "',...
        num2str(cCondition),'/',num2str(NumConditions),'"');
        handles.text_status.String = sprintf(statusStr);
    else
        dispstat(strcat('saving fluorescence intensity files...',...
        num2str((ChunkSizeTotal/wbEndPoint)*100), '% complete for "', CellName, '" at condition "',...
        num2str(cCondition),'/',num2str(NumConditions),'"'))
    end
end

if cParti == 1
    cc512 = bwconncomp(bw512);
    NumROIs = cc512.NumObjects;
    if NumROIs ~= cellcounter
        errordlg('ERROR: problem resizing ROIs - resized ROIs now overlap');
    end
end
SaveSwitch = 0;
if cParti == numel(NumStrFilesList)
    SaveSwitch = 1;
elseif cParti < numel(NumStrFilesList)
    if NumStrFilesList(cParti) ~= NumStrFilesList(cParti+1)
        SaveSwitch = 1;
    end
end
if SaveSwitch == 1
%CG: in this case, FluoInt is on the last partitioned file in the current
%condition and it is time to save the data
    cd(fIntensity_path); CoreType = strcat(CellName,suffixStrs{cCondition},'_IntensityFiles.mat');
    try
        ints = load(CoreType);
        if ~isempty(ints)
            if strcmp(DataType,'Raw') 
                save(CoreType, 'CA_IntensityRaw','-append')
            elseif strcmp(DataType,'170')
                save(CoreType, 'CA_Intensity170','-append')
            elseif strcmp(DataType,'both')
                save(CoreType, 'CA_Intensity170','CA_IntensityRaw','-append')
            end
        end
    catch
        if strcmp(DataType,'Raw') 
            save(CoreType, 'CA_IntensityRaw')
        elseif strcmp(DataType,'170')
            save(CoreType, 'CA_Intensity170')
        elseif strcmp(DataType,'both')
            save(CoreType, 'CA_Intensity170','CA_IntensityRaw')
        end
    end
end

fio = struct;

fio.CA_Intensity170 = CA_Intensity170;
fio.CA_IntensityRaw = CA_IntensityRaw;
fio.CA_ROIs = CA_ROIs;
fio.ROI_map = ROI_map;
fio.ROI_Idx_map_Raw = ROI_Idx_map_Raw;
fio.ChunkSizeTotal = ChunkSizeTotal;
fio.fmedian = fmedian;
fio.wbStep = wbStep;

