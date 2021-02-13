function BranchTracks(datadirs, keywordStrs, expDirIndex, branchIDs,...
    customStartPointDistances, configStruct)
%CG: BranchTracks organizes ROIs according to their location in preparation
%for graphing and analyzing the data (assuming the user wishes to use ROI
%location as an interesting variable). BranchTracks creates
%'..._BranchOrderedROIs.mat' file containing:
%
% - CA_BranchGroups: cell array with each cell containing a structure
% holding data about each roi on a single branch. The structure will have the following
% fields: 'roi'; 'amps'; 'growths'; 'decays'; 'distance'. e.g. 'amps' is a
% numeric array containing the peak amplitude measurements for all events
% not rejected by the user. 
%
% - BranchGroupMtx: 2D matrix containing the branch ID numbers (see
% CA_Bones; '..._CellBonesCA.mat'for definition). Connected branches are aligned in rows
% with primary (i.e. more proximal) branches in the first column and
% distal-most branches in the last column. e.g. BranchGroupMtx(1,:) = [11,
% 12, 14]; In this case, branch number 14 branches off branch number 12,
% which in turn branches from no. 11. Branch 11 is directly connected to
% the cell body.
%
% - BranchPntDistMtx: contains the maximum length of each branch from the
% cell body. Both BranchGroupMtx and BranchPntDistMtx are used to index the
% data in CA_BranchGroups, so we can find the max distance of a branch with
% ID stored in BranchGroupMtx(2,2) by going to BranchPntDistMtx(2,2). 
%
%The aim of BranchTracks is organize the ROIs such that each cell in
%CA_BranchGroups contains all rois on a single branch and all structures in
%the same row of CA_BranchGroups will be linearly connected. BranchTracks
%also organizes the data such that proximal branches are aligned with
%distal branches to make the longest tracks.

% Author: Charlie J Gilbride
% version: 1.0.0

% Updates:
% version: 1.0.1 - names of experiments (or cells) contributing to the
% n-values are saved to the parent directory. This enables users to know
% contribution from independent cell cultures. These names are saved to the
% nested cell array Ncells

if nargin == 0
    load('directories.mat','datadirs')
    %datadirs = uipickfiles('prompt','Select data folder(s)');
    %datadirs = {}; datadirs{1} = '/Volumes/No2/MATLAB/AARG_Data/ReAnalysis/TTAP2/CG0408172';
    %datadirs{1} = '/Volumes/No2/MATLAB/AARG_Data/ReAnalysis/TTAP2/CG2508172';
    suffixStrs = {'a_baseline','b_TTAP2'}; 
    
    binWidth = 20; %CG: length of dendrite (in µm) along which rois will be binned.

    resc = 3; %CG: if the resolution was not set correctly in CellBones, use resc to 
    %rescale to the correct value. e.g. the raw data is a 512*512 image. This
    %has a resolution of 6.25pixels/µm. If the images are resized to 170*170,
    %then the resolution is 2.0833pixels/µm. If the distance calculation is
    %made without taking into account the resizing, then an appropriate value
    %for resc would be 3. If the distance calculation already takes into
    %account resizing, then resc should be 1. If resc is not 1 then a warning
    %will appear. 
    
    resc_dir = [1 2]; %CG: resc_dir is the directory or directories which will be rescaled
    %using the value resc.

    if resc ~= 1; disp('WARNING: distance values are being rescaled.'); end
    
    customStartPointDistances = [83.2, 146.1, 138]; expDirIndex = [2, 4, 4];
    branchIDs = [1, 7, 1];
    
else
    suffixStrs = keywordStrs;
%CG: rescaling needs to be applied in the case of CG0408172 and CG0408173 for the TTAP2 group.    
%     if strcmp(suffixStrs{end},'b_TTAP2'); resc_dir = [1 2]; resc = 3; else resc_dir = []; resc = 1; end
    binWidth = 20; %CG: default value for configStruct.lengthOfDendritePerGraph
%     if resc ~= 1; disp('WARNING: distance values are being rescaled.'); end
end


%CG: customStartPointDistances is the start point distance from the soma. This value
%should be found using the Fiji's segmented line tool. expDirIndex is
%the data directory to which the corresponding value in customStartPointDistances
%applies. e.g. if customStartPointDistances = [84.2 51.9]; and expDirIndex = [3 7];
%this will mean that the 3rd directory in the list (datadirs) has a start
%point of 84.2 microns and the 7th directory in the list has a start point
%of 51.9 microns. branchIDs specifies precisely which branch this applies
%to. branchIDs is useful if there is more than one start point for an
%experiment where the cell body is not visible. For example: if 
%customStartPointDistances = [84.2 51.9 75.3]; expDirIndex = [3 7 7];
%branchIDs = [1 1 2] would mean the experiment in directory 3 has one start
%point with a non-zero start point value and there is a single branch
%starting from this non-zero start point with a branch ID of one. The
%experiment in the 7th directory has two branches starting at non-zero
%start points - one with a branch ID of 1 and the other with a branch ID of
%2. Branch 1 starts 51.9 unit distances from the soma, while branch 2
%starts 75.3 unit distances from the soma. 


if sum(~isnan(customStartPointDistances)) > 0 || sum(~isnan(expDirIndex)) > 0
    if numel(customStartPointDistances) ~= numel(expDirIndex)
        disp('variables inappropriately set. Analysis cancelled')
    else
        disp('start points distances not set to 0')
    end
elseif sum(isnan(customStartPointDistances)) ~= sum(isnan(expDirIndex))
    disp('variables inappropriately set. Analysis cancelled')
end
    
% plotMeanOfMedians = 0; 
% %CG: plot mean of medians should be made an input argument and set by the
% %user in the configureGraph GUI. 

if iscell(datadirs) && ~isempty(datadirs)
    NumDirs = size(datadirs,2); dirdata = struct; grand_max = []; NumConditions = size(suffixStrs,2);
    for cDirIdx = 1 : NumDirs
        cd(datadirs{cDirIdx}); [~, CellName, ~] = fileparts(datadirs{cDirIdx});

        for cCondIdx = 1 : NumConditions
            try
%                 BreakTry
                branchorders = load(strcat(CellName,suffixStrs{cCondIdx},'_BranchOrderedROIs.mat'));
                CA_BranchGroups = branchorders.CA_BranchGroups;
                BranchPntDistMtx = branchorders.BranchPntDistMtx;
                BranchGroupMtx = branchorders.BranchGroupMtx;
                grand_max = branchorders.grand_max;

            catch
                bones = load(strcat(CellName,'_CellBonesCAs.mat'));
                CA_ROIc = bones.CA_ROIc; NumTotalROIs = size(CA_ROIc,1);
                branchpool = []; 
                max_contig = 0;
                for cROIIdx = 1 : NumTotalROIs
%CG: CA_ROIc{cROIIdx,8} contains the branch numbers that connect the
%current ROI to a start point. The largest value in branchpool after this
%for-loop will indicate the total number of branch levels for the current cell.
                    if ~isempty(CA_ROIc{cROIIdx,8})
                        branchpool = [branchpool,CA_ROIc{cROIIdx,8}];
                        prov_max_contig = numel(CA_ROIc{cROIIdx,8});
                        if prov_max_contig > max_contig
                            max_contig = prov_max_contig;
                        end
                    end
                end
%CG: max_contig also indicates the most distal branch.
                blkmtx = zeros(NumTotalROIs,1); CA_BranchOrg = {}; cntr = 1;
%CG: blkmtx prevents the double for-loop below looking up the same ROI
%more than once. CA_BranchOrg assists with organising the branch data.
                for cBLidx = 1 : max_contig
                    BranchOrgMtx = []; bpd = [];
                    for cROIIdx = 1 : NumTotalROIs
                        
                        if numel(CA_ROIc{cROIIdx,8}) == cBLidx && blkmtx(cROIIdx,1) ~= 1
                            branchLineage = CA_ROIc{cROIIdx,8}; startBranch = branchLineage(end);
                            blkmtx(cROIIdx,1) = 1; raw_dist = CA_ROIc{cROIIdx,9};
%                             if ~isempty(find(resc_dir == cDirIdx,1))
%                                 raw_dist = raw_dist*resc;
%                             end

                            if ~isempty(find(expDirIndex == cDirIdx,1))
                                aa = expDirIndex == cDirIdx;
                                bb = branchIDs == startBranch;
                                cc = sum(cat(1,aa,bb),1) == 2;
%CG: find which start point-connected branch the current ROI is directly or
%indirectly connected to.
                                raw_dist = raw_dist+customStartPointDistances(cc);
                            end
                            
                            if isempty(BranchOrgMtx)
                                BranchOrgMtx = CA_ROIc{cROIIdx,8};
                                bpd = raw_dist;
                            else
                                BranchOrgMtx = cat(1,BranchOrgMtx,CA_ROIc{cROIIdx,8});
                                bpd = cat(1,bpd,raw_dist);
                            end

                        end
                    end
%CG: identify unique combinations of branches.
                    Sz_BranchOrgMtx = size(BranchOrgMtx); blkmtx2 = zeros(Sz_BranchOrgMtx); 
                    uBranchOrgMtx = []; max_ubpd = []; cntr = 1;

                    while ~isempty(find(blkmtx2==0,1))
                        cCombo = BranchOrgMtx(cntr,:); %branchpnts = bpd(cntr,:);
                        idces = find(blkmtx2(:,1) == 0); switch1 = 0; ubpd = [];

                        for cidx = 1 : numel(idces) 
                            if sum(BranchOrgMtx(idces(cidx),:) - cCombo) == 0 && switch1 == 0
                                switch1 = 1; uBranchOrgMtx = [uBranchOrgMtx;cCombo];
                                ubpd = [ubpd; bpd(idces(cidx),:)]; blkmtx2(idces(cidx),:) = 1;
                            elseif sum(BranchOrgMtx(idces(cidx),:) - cCombo) == 0 && switch1 == 1
                                blkmtx2(idces(cidx),:) = 1; 
                                ubpd = [ubpd; bpd(idces(cidx),:)];
                            end
                        end
                        
                        if ~isempty(ubpd)
                            if size(ubpd,1) == 1 && size(max_ubpd,2) == size(ubpd,2)
                                max_ubpd = [max_ubpd; ubpd];
                            else 
                                max_ubpd = [max_ubpd; max(ubpd)];
                            end
                        end
                        cntr = cntr + 1;
                    end
                    CA_bpd{cBLidx} = max_ubpd;
                    CA_BranchOrg{cBLidx} = uBranchOrgMtx;
                end

                CA_BranchGroups = {};

                NumContigGrps = size(CA_BranchOrg,2); BranchGroupMtx = [];
                BranchPntDistMtx = [];
%CG: number of groups of contiguous branches. Maximum number of contiguous
%branches is usually 3 and is unlikely to exceed 4. 
                for cGrpIdx = 1 : NumContigGrps
                    cGrp = CA_BranchOrg{NumContigGrps-(cGrpIdx-1)}; Sz_cGrp = size(cGrp);
                    cbpGrp = CA_bpd{NumContigGrps-(cGrpIdx-1)}; 
                    for cRow = 1 : Sz_cGrp(1)
                        for cCol = 1 : Sz_cGrp(2)
                            if isempty(BranchGroupMtx)
                                BranchGroupMtx = cGrp(cRow,Sz_cGrp(2)-(cCol-1));
                                BranchPntDistMtx = cbpGrp(cRow,Sz_cGrp(2)-(cCol-1));
                            else
                                if (cRow == 1 && cCol == 1) || (cRow > 1 && cCol == 1)
                                    BranchGroupMtx = cat(1,BranchGroupMtx,zeros(1,max_contig));
                                    BranchGroupMtx(end,cCol) = cGrp(cRow,Sz_cGrp(2)-(cCol-1));

                                    BranchPntDistMtx = cat(1,BranchPntDistMtx,zeros(1,max_contig));
                                    BranchPntDistMtx(end,cCol) = cbpGrp(cRow,Sz_cGrp(2)-(cCol-1));
                                else
                                    BranchGroupMtx(end,cCol) = cGrp(cRow,Sz_cGrp(2)-(cCol-1));
                                    BranchPntDistMtx(end,cCol) = cbpGrp(cRow,Sz_cGrp(2)-(cCol-1));
                                end

                                if cCol == Sz_cGrp(2) 
%CG: if about to enter a new group with more that one branch, a blank layer
%should be added to the matrix. This can be used in case there are ROIs
%found above the branch point, but do not belong to any other contiguous
%group of branches. 
                                end
                            end
                        end
                    end
                end

                lglidx = BranchGroupMtx > 0;

                Sz_BranchGroupMtx = size(BranchGroupMtx); NumCombos = Sz_BranchGroupMtx(1);
                CAIdces = cell2mat(CA_ROIc(:,3)); CA_dists = CA_ROIc(:,10);
                empties = cellfun('isempty',CA_dists); CA_dists(empties) = {NaN};
                DistVals = cell2mat(CA_dists);
                data = load(strcat(CellName,suffixStrs{cCondIdx},'_MeasuresRAW.mat'));
                parameters = data.parameters;
                for cComboIdx = 1 : NumCombos
                    cCombo = BranchGroupMtx(cComboIdx,:); cCombo = cCombo(cCombo>0); IdxList = [];
                    for cROIIdx = 1 : NumTotalROIs
                        if numel(CA_ROIc{cROIIdx,8}) == numel(cCombo>0)
                            if sum(CA_ROIc{cROIIdx,8} - cCombo) == 0
                                IdxList = [IdxList, cROIIdx];
                            end
                        end
                    end
                    cIdces = CAIdces(IdxList); Sz_cIdces = size(cIdces); bstr = struct;
                    NumROIs = Sz_cIdces(1)*Sz_cIdces(2);
                    for cItemIdx = 1 : NumROIs
                        roi_lgic = parameters(1,:) == cIdces(cItemIdx);
                        AllAmps = parameters(2,:); amps = AllAmps(roi_lgic);
                        AllGrowths = parameters(3,:); growths = AllGrowths(roi_lgic);
                        AllDecays = parameters(4,:); decays = AllDecays(roi_lgic);
                        distance = DistVals(CAIdces == cIdces(cItemIdx)); 
                        AllXCos = parameters(5,:); xcos = AllXCos(roi_lgic);
                        AllYCos = parameters(5,:); ycos = AllYCos(roi_lgic);
                        
%                         if ~isempty(find(resc_dir == cDirIdx,1))
%                             distance = distance*resc;
%                         end
                        
                        if ~isempty(find(expDirIndex == cDirIdx,1))
                            aa = expDirIndex == cDirIdx;
                            bb = branchIDs == startBranch;
                            cc = sum(cat(1,aa,bb),1) == 2;
                            distance = distance+customStartPointDistances(cc);
                        end

                        bstr(cItemIdx).roi = cIdces(cItemIdx);
                        bstr(cItemIdx).amps = amps;
                        bstr(cItemIdx).growths = growths;
                        bstr(cItemIdx).decays = decays;
                        bstr(cItemIdx).distance = distance;
                        bstr(cItemIdx).xcos = xcos;
                        bstr(cItemIdx).ycos = ycos;
                    end
%CG: when cCombo = [11 12 14]; this means the ROIs, which were just
%evaluated, are located on branch 14. Data from these ROIs should therefore
%be stored (at least temporarily) in the CA_BranchGroups cell array, which
%has the same coordinates as the ID no. 14 in BranchGroupMtx.
                    CA_BranchGroups{cComboIdx,numel(cCombo)} = bstr;
                end
                trackarray = (max_contig:-1:1);
%CG: sort the data in CA_BranchGroups according to the branch point
%distance data in BranchPntDistMtx and the branch ID number in
%BranchGroupMtx.
                for cColIdx = 1 : numel(trackarray)
%CG: all ROIs in the most distal branches will already be stored in the
%correct location in CA_BranchGroups, so there is no need to process them
%further. Hence, cColIdx = 2 : max_contig rather than cColIdx = 1 :
%max_contig.
                    colnum = trackarray(cColIdx);

%CG: evaluate columns in reverse order.

                    lgl = BranchGroupMtx > 0; cComboIdx = 1;
                    while cComboIdx <= NumCombos
                        if BranchGroupMtx(cComboIdx,colnum) > 0 
%CG: data are organised such that single branches with no branch points are
%put in the lower rows of the BranchGroupMtx and BranchPntDistMtx. These
%single branches do not need to be organised any further.
                            cCombo = BranchGroupMtx(cComboIdx,:); %cCombo_s=cCombo(cCombo>0); 
                            cbranchID = cCombo(colnum);
                            
                            IDidx = find(BranchGroupMtx == cbranchID); 
                            [MaxBranchPnt,~] = max(BranchPntDistMtx(IDidx));
                            if cCombo(end) == 0 %&& numel(cCombo(cCombo>0)) > 1
%CG: the coordinates of the arrays that must be evaluated will always have 0-valued 
%elements next to them in BranchGroupMtx.
                                maxbp_idx = find(BranchPntDistMtx==MaxBranchPnt);
                                bstr = CA_BranchGroups{cComboIdx,colnum}; Sz_bstr = size(bstr);
                                NumROIsOnBranch = Sz_bstr(1)*Sz_bstr(2); new_bstr = struct;
                                new_cntr = 1;
                                if numel(IDidx) > 1
                                    delArray = [];
                                    for cROIIdx = 1 : NumROIsOnBranch
                                        cROIdist = bstr(cROIIdx).distance;
                                        if cROIdist < MaxBranchPnt+0.5
        %CG: providing the ROI is within 500nm of the maximum branch length, it
        %will also be included. This is to avoid ROIs being kicked off a branch
        %just because the connection point between two branch lines falls just
        %below the connection point between the green line (spine line) and the
        %white line (for the dendrite). 
                                            new_bstr(new_cntr).roi = bstr(cROIIdx).roi;
                                            new_bstr(new_cntr).amps = bstr(cROIIdx).amps;
                                            new_bstr(new_cntr).growths = bstr(cROIIdx).growths;
                                            new_bstr(new_cntr).decays = bstr(cROIIdx).decays;
                                            new_bstr(new_cntr).distance = bstr(cROIIdx).distance;
                                            new_bstr(new_cntr).xcos = bstr(cROIIdx).xcos;
                                            new_bstr(new_cntr).ycos = bstr(cROIIdx).ycos;
                                            new_cntr = new_cntr + 1;
                                            delArray = [delArray,cROIIdx];
                %                             bstr(cROIIdx).roi = []; bstr(cROIIdx).amps = [];
                %                             bstr(cROIIdx).growths = []; bstr(cROIIdx).decays = [];
                %                             bstr(cROIIdx).distance = [];
                                        end

                                    end
                                    NumEntries2Del = size(delArray,2);
                                    if cComboIdx == 14
                                        stophere = 1;
                                    end
                                    for cEntry = 1 : NumEntries2Del
                                        bstr(delArray(cEntry)).roi = []; bstr(delArray(cEntry)).amps = [];
                                        bstr(delArray(cEntry)).growths = []; bstr(delArray(cEntry)).decays = [];
                                        bstr(delArray(cEntry)).distance = []; bstr(delArray(cEntry)).xcos = [];
                                        bstr(delArray(cEntry)).ycos = [];
                                    end
                                    if ~isempty(delArray)
                                        roi_out = {bstr(~cellfun(@isempty,{bstr.roi})).roi};
        %CG: some ROIs may not have been excluded, but may not contain measurable data. 
        %So, there will be a mismatch between the cells in roi_out and amps_out
        %unless the same bstr.roi are indexed out for the bstr.amps field. 
                                        amps_out = {bstr(~cellfun(@isempty,{bstr.roi})).amps};
                                        growths_out = {bstr(~cellfun(@isempty,{bstr.roi})).growths};
                                        decays_out = {bstr(~cellfun(@isempty,{bstr.roi})).decays};
                                        distance_out = {bstr(~cellfun(@isempty,{bstr.roi})).distance};
                                        xcos_out = {bstr(~cellfun(@isempty,{bstr.roi})).xcos};
                                        ycos_out = {bstr(~cellfun(@isempty,{bstr.roi})).ycos};

                                        sz_out = size(roi_out); OutSize = sz_out(1)*sz_out(2);

        %CG: if all parts of the roi field are empty, then make the appropriate
        %cell in CA_BranchGroups empty. 
                                        if isempty(roi_out)
                                            CA_BranchGroups(cComboIdx,colnum) = {[]};
                                        else
                                            bstr = struct;
                                            for cCell = 1 : OutSize
                                                bstr(cCell).roi = roi_out{cCell};
                                                bstr(cCell).amps = amps_out{cCell};
                                                bstr(cCell).growths = growths_out{cCell};
                                                bstr(cCell).decays = decays_out{cCell};
                                                bstr(cCell).distance = distance_out{cCell};
                                                bstr(cCell).xcos = xcos_out{cCell};
                                                bstr(cCell).ycos = ycos_out{cCell};
                                            end
                                            CA_BranchGroups{cComboIdx,colnum} = bstr;
                                        end
                                    end
                                    if NumEntries2Del > 0
        %CG: if there is more than one branch point with a maximum branch point value, then the index at the 
        %top of the maxbp_idx array should be taken. This is appropriate because it
        %is desirable to have the largest track possible and then larger the track
        %the higher up (i.e. smaller row number) it is placed in the
        %CA_BranchGroups cell array. 
        %                                 CA_BranchGroups{maxbp_idx(1)} = new_bstr;
        %CG: More branch points ~= greater distance. Must go through each of the
        %indexed rows and find the one with the greatest maximum distance
                                        cmax = [];
                                        [rr,~] = ind2sub(size(BranchPntDistMtx),maxbp_idx);
                                    
                                        for cItemIdx = 1 : numel(maxbp_idx)
                                            branchpnts = BranchPntDistMtx(rr(cItemIdx),:);
                                            bpnts_idx = find(branchpnts > 0); cmax = [cmax;branchpnts(bpnts_idx(end))]; 
                                        end
                                        [~,gmax_idx] = max(cmax);
                                        CA_BranchGroups{maxbp_idx(gmax_idx)} = new_bstr;

                                    end
                                end
                                
                            end
                            grand_max = [grand_max,MaxBranchPnt];
                        end
                        cComboIdx = cComboIdx + 1;
                    end
                end
                save(strcat(CellName,suffixStrs{cCondIdx},'_BranchOrderedROIs.mat'),'CA_BranchGroups',...
                    'BranchPntDistMtx', 'BranchGroupMtx', 'grand_max')
%                 save(strcat(CellName,suffixStrs{cCondIdx},'_BranchOrderedROIs.mat'),'grand_max', '-append')
            end
            dirdata(cDirIdx,cCondIdx).bgrps = CA_BranchGroups;
            dirdata(cDirIdx,cCondIdx).distmtx = BranchPntDistMtx;
            dirdata(cDirIdx,cCondIdx).bgrpmtx = BranchGroupMtx;
        end
        
    end
    grand_max = max(grand_max); distanceBins = configStruct.lengthOfDendritePerGraph;
    if ischar(distanceBins)
        if strcmp(distanceBins,'default')
            distanceBins = binWidth;
        end
    end
    [CA_Bins,SizeBins] =...
        CompileTracks(dirdata,NumDirs,distanceBins,grand_max,datadirs,suffixStrs);
    
    if configStruct.dontDisplayMedianPlots == 1 && configStruct.dontDisplayHistograms == 0
        [CA_BinMedians] = AddMedianField(CA_Bins,suffixStrs,NumDirs);
        GraphMedianData(CA_BinMedians,suffixStrs,datadirs,SizeBins,configStruct)
    elseif configStruct.dontDisplayMedianPlots == 0 && configStruct.dontDisplayHistograms == 1
        GraphBinROIs(CA_Bins,suffixStrs,datadirs,SizeBins,configStruct);
    elseif configStruct.dontDisplayMedianPlots == 0 && configStruct.dontDisplayHistograms == 0
        GraphBinROIs(CA_Bins,suffixStrs,datadirs,SizeBins,configStruct);
        [CA_BinMedians] = AddMedianField(CA_Bins,suffixStrs,NumDirs);
        GraphMedianData(CA_BinMedians,suffixStrs,datadirs,SizeBins,configStruct)
    end
end

function [CA_Bins,SizeBins] =...
    CompileTracks(dirdata,NumDirs,distanceBins,grand_max,datadirs,suffixStrs)

rem = mod(ceil(grand_max), distanceBins); LastEdge = ceil(grand_max)-rem;
SizeBins = (0:distanceBins:LastEdge+distanceBins); NumBins = numel(SizeBins);
NumConditions = size(suffixStrs,2); CA_Bins = cell(NumConditions,NumBins-1);
for cDirIdx = 1 : NumDirs
    cd(datadirs{cDirIdx}); [~, CellName, ~] = fileparts(datadirs{cDirIdx});
    for cCondIdx = 1 : NumConditions
        CA_BranchGroups = dirdata(cDirIdx,cCondIdx).bgrps; 
        BranchPntDistMtx = dirdata(cDirIdx,cCondIdx).distmtx;
        BranchGroupMtx = dirdata(cDirIdx,cCondIdx).bgrpmtx;

        Total_R = size(CA_BranchGroups,1); Total_C = size(CA_BranchGroups,2); 
        for cRow = 1 : Total_R
            for cCol = 1 : Total_C
                if ~isempty(CA_BranchGroups{cRow,cCol})
                    cstr = CA_BranchGroups{cRow,cCol}; NumROIs = size(cstr,2);
                    for cROIIdx = 1 : NumROIs
                        
                        cdist = cstr(cROIIdx).distance; roibinfound = 0; 
                        cSizeIdx = 2; nContribution = [];
                        while roibinfound == 0
                            upperb = SizeBins(cSizeIdx); lowerb = SizeBins(cSizeIdx-1); 
%CG: find current ROI's distance in range of lowerb ? x < upperb
                            if lowerb <= cdist && cdist < upperb
                                roibinfound = 1; 
                                nContribution = [cDirIdx BranchGroupMtx(cRow,cCol)];
                            else
                                cSizeIdx = cSizeIdx + 1;
                            end
                        end
                        if isempty(CA_Bins{cCondIdx,cSizeIdx-1})
                            binstr = struct;
                            binstr.roi = cstr(cROIIdx).roi;
                            binstr.amps = cstr(cROIIdx).amps;
                            binstr.growths = cstr(cROIIdx).growths;
                            binstr.decays = cstr(cROIIdx).decays;
                            binstr.distance = cstr(cROIIdx).distance;
                            binstr.contribution = nContribution;
                            %CG: v1.0.1
                            binstr.cellName = CellName;
                        elseif ~isempty(CA_Bins{cCondIdx,cSizeIdx-1})
                            binstr = CA_Bins{cCondIdx,cSizeIdx-1};
                            NumPreROIs = size(binstr,2);
                            binstr(NumPreROIs+1).roi = cstr(cROIIdx).roi;
                            binstr(NumPreROIs+1).amps = cstr(cROIIdx).amps;
                            binstr(NumPreROIs+1).growths = cstr(cROIIdx).growths;
                            binstr(NumPreROIs+1).decays = cstr(cROIIdx).decays;
                            binstr(NumPreROIs+1).distance = cstr(cROIIdx).distance;
                            binstr(NumPreROIs+1).contribution = nContribution;
                            %CG: v1.0.1
                            binstr(NumPreROIs+1).cellName = CellName;
                        end
                        CA_Bins{cCondIdx,cSizeIdx-1} = binstr;
                    end
                end
            end
        end
    end
end

function GraphBinROIs(CA_Bins,suffixStrs,datadirs,SizeBins,configStruct)

NumConditions = size(suffixStrs,2); 
screenSize = get(groot,'Screensize'); fontSize = 16;
positionScaled = 0.75; %CG: reduce size of the graph figure (easier to manipulate figure after presentation on some platforms).

legendLabels = cell(size(suffixStrs));
for cCondIdx = 1 : NumConditions
    %CG: any underscore in the condition name messes up the presentation of
    %the word when it is put into a graph, so it is better to have legend
    %labels that don't have an underscore.
    defaultConditionName = suffixStrs{cCondIdx};
    underScorePosition = strfind(defaultConditionName, '_');
    if ~isempty(underScorePosition)
        legendLabels{cCondIdx} = defaultConditionName(underScorePosition+1:end);
    else
        legendLabels{cCondIdx} = defaultConditionName;
    end
end
if iscell(datadirs) && ~isempty(datadirs)
    
    xAxisLabel = ''; conditionString = strcat('. Condition: "', suffixStrs{cCondIdx},'"');
    if configStruct.dontDisplayHistograms == 0 && configStruct.displayPeakAmplitude == 1
        fh_hist = figure('Name',strcat('Distribution of individual event peak amplitudes separated according to distance from soma',conditionString));
        xAxisLabel = 'Peak amplitude (a.u.)';
    elseif configStruct.dontDisplayHistograms == 0 && configStruct.displayDecayTimeConstant == 1
        fh_hist = figure('Name',strcat('Distribution of individual event decay time constants separated according to distance from soma', conditionString));
        xAxisLabel = 'tau decay constant (s)';
    end
    
    fh_hist.NumberTitle = 'off'; fh_hist.Position = screenSize.*positionScaled; 
    set(findall(fh_hist, '-property', 'Units'), 'Units', 'Normalized')
    set(fh_hist, 'Resize','on')
   
    if isempty(configStruct.maximumDistance)
        NumBins = size(CA_Bins,2); 
    elseif ~isempty(configStruct.maximumDistance)
        NumBins = size(CA_Bins,2); cBinIdx = 1;
        %CG: if maximum distance is non-empty, then the user has
        %specified that s/he only wants to see graphs with distances up
        %to a value of configStruct.maximumDistance
        while cBinIdx <= NumBins
            numberOfStructs = size(CA_Bins{cCondIdx,cBinIdx},2);
            distances = CA_Bins{cCondIdx,cBinIdx}; structIdx = 1;
            while structIdx <= numberOfStructs 
                if isempty(find(distances(structIdx).distance < configStruct.maximumDistance,1))
                    NumBins = cBinIdx-1; structIdx = numberOfStructs;
                end
                structIdx = structIdx + 1;
            end
            cBinIdx = cBinIdx + 1;
        end
    end

    ansx = NumBins/2; 
    
    Ncells = cell(NumBins,NumConditions);
    Ndirs = cell(NumBins,NumConditions); 
    Nbranches = cell(NumBins,NumConditions);
    ROICounts = cell(NumBins,NumConditions);
    
    for cCondIdx = 1 : NumConditions
            
        for cBinIdx = 1 : NumBins

            titleStr = strcat(num2str(SizeBins(cBinIdx)),...
                '-',num2str(SizeBins(cBinIdx+1)),'µm');
                
            binstr = CA_Bins{cCondIdx,cBinIdx}; NumROIs = size(binstr,2); 
            if configStruct.displayPeakAmplitude == 1; amps = []; end
            if configStruct.displayDecayTimeConstant == 1; decays = []; end
            if configStruct.displayNValues == 1; nContributions = []; nCellNames = {}; end
            ROICount = 0; 
            for cROIIdx = 1 : NumROIs
%                 if cROIIdx == 53
%                     stophere = 1;
%                 end
                if configStruct.displayPeakAmplitude == 1;
                    amps = [amps,binstr(cROIIdx).amps];
                    if configStruct.displayNValues == 1 && sum(~isnan(binstr(cROIIdx).amps)) > 0;
                        nContributions = [nContributions; binstr(cROIIdx).contribution];
                        if isempty(nCellNames); nCellNames{1} = binstr(cROIIdx).cellName; else nCellNames{end+1,1} = binstr(cROIIdx).cellName; end
                        ROICount = ROICount + 1;
                    end
                end
                if configStruct.displayDecayTimeConstant == 1;
                    decays = [decays,binstr(cROIIdx).decays];
                    if configStruct.displayNValues == 1 && sum(~isnan(binstr(cROIIdx).decays)) > 0;
                        nContributions = [nContributions; binstr(cROIIdx).contribution];
                        if isempty(nCellNames); nCellNames{1} = binstr(cROIIdx).cellName; else nCellNames{end+1,1} = binstr(cROIIdx).cellName; end
                        ROICount = ROICount + 1;
                    end
                end
                
            end
            ROICounts{cBinIdx,cCondIdx} = ROICount;
            if configStruct.displayDecayTimeConstant == 1;
                decaysLogicFilter = decays == 5000000000000000;
                decays(decaysLogicFilter) = NaN;
                %CG: Some events with very slow decay times may have curve fits
                %with exponential decays that are extremely large. In such cases, the
                %most appropriate action would be to discard the event or
                %accept only its amplitude measurement. The decay of such
                %events may be distorted by a second, overlapping event, but
                %this may not have been apparent when classifying the events.
                %It seems any such events have decay constants of
                %"5000000000000000". 
            end
%CG: plot histograms
            
            notEnoughData = 0;
            if configStruct.displayPeakAmplitude == 1
                if sum(~isnan(amps)) == 0; notEnoughData = 1; end
            end
            if configStruct.displayDecayTimeConstant == 1;
                if sum(~isnan(decays)) == 0; notEnoughData = 1; end
            end
            
            if notEnoughData == 0
            
                if configStruct.displayNValues == 1;
                    unique_nCellNames = {}; %CG: v1.0.1
                    uniqueNDirs = unique(nContributions(:,1));
                    Ndirs{cBinIdx,cCondIdx} = numel(uniqueNDirs); NbranchesArray = [];
                    for cDirIdx = 1 : Ndirs{cBinIdx,cCondIdx}
                        cDir = uniqueNDirs(cDirIdx); directorySelection = nContributions(:,1);
                        branchIDs = nContributions(:,2);
                        targetIndices = directorySelection == cDir;
                        NbranchesArray = [NbranchesArray, numel(unique(branchIDs(targetIndices)))];
                        
                        % CG: also retrieve names for all contibuting
                        % cells for update v1.0.1
                        targetNames = nCellNames(targetIndices);
                        unique_nCellNames(cDirIdx,1) = targetNames(1);
                       
                    end
                    Nbranches{cBinIdx,cCondIdx} = sum(NbranchesArray);
                    Ncells{cBinIdx,cCondIdx} = {unique_nCellNames};
                end


                if configStruct.dontDisplayHistograms == 0 
                    figure(fh_hist); subplot(2,ceil(ansx),cBinIdx);
                    if cCondIdx == 1
                        if configStruct.displayPeakAmplitude == 1; 
    %                         h1 = histogram(amps, 'DisplayStyle', 'stairs'); 
                            ampsNaNLogic = isnan(amps); amps(ampsNaNLogic) = [];
                            h1 = histogram(amps, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
                        end
                        if configStruct.displayDecayTimeConstant == 1; 
    %                         h1 = histogram(decays, 'DisplayStyle', 'stairs');
                            decaysNaNLogic = isnan(decays); decays(decaysNaNLogic) = [];
                            h1 = histogram(decays, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
                        end

                        if cCondIdx == 1; h1.EdgeColor = [0 0 1]; end

                    elseif cCondIdx > 1
                        ax1 = gca; hold(ax1,'on');

                        if configStruct.displayPeakAmplitude == 1; 
    %                         h1 = histogram(amps, 'DisplayStyle', 'stairs'); 
                            ampsNaNLogic = isnan(amps); amps(ampsNaNLogic) = [];
                            h1 = histogram(amps, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
                        end
                        if configStruct.displayDecayTimeConstant == 1; 
    %                         h1 = histogram(decays, 'DisplayStyle', 'stairs'); 
                            decaysNaNLogic = isnan(decays); decays(decaysNaNLogic) = [];
                            h1 = histogram(decays, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
                        end

                        if configStruct.displayNValues == 1;
                            firstDir = datadirs{1}; slashes = strfind(datadirs{1}, '/'); parentDir = firstDir(1:slashes(end)); 
                            cd(parentDir) 
                            save('NamesOfExperimentsForN.mat', 'Ncells')
                            
                            labelNdirs = strcat('N = "',num2str(Ndirs{cBinIdx,1}),'/', num2str(Ndirs{cBinIdx,cCondIdx}), '" experiments');
                            labelNbranches = strcat('N = "',num2str(Nbranches{cBinIdx,1}), '/', num2str(Nbranches{cBinIdx,cCondIdx}), '" branches');
                            labelNrois = strcat('N = "',num2str(ROICounts{cBinIdx,1}), '/', num2str(ROICounts{cBinIdx,cCondIdx}), '" ROIs');

                            xlims1 = (ax1.XLim); ylims1 = (ax1.YLim);
                            t1xy = [xlims1(2)*0.4 ylims1(2)*0.7];
                            text(t1xy(1),t1xy(2),labelNdirs, 'FontSize', fontSize)
                            t2xy = [xlims1(2)*0.4 ylims1(2)*0.6];
                            text(t2xy(1),t2xy(2),labelNbranches, 'FontSize', fontSize)
                            t3xy = [xlims1(2)*0.4 ylims1(2)*0.5];
                            text(t3xy(1),t3xy(2),labelNrois, 'FontSize', fontSize)
    %                         t4xy = [xlims1(2)*0.4 ylims1(2)*0.5];
    %                         text(t4xy(1),t4xy(2),labelRanksum, 'FontSize', fontSize)
    %                         t5xy = [xlims1(2)*0.4 ylims1(2)*0.4];
    %                         text(t5xy(1),t5xy(2),labelZVal, 'FontSize', fontSize)
                        end

                        if cCondIdx ~= 1; h1.EdgeColor = [1 0 0]; end
    %                     if cCondIdx ~= 1; h1.EdgeColor = [0.6 0.6 0.6]; end
    %                     if cCondIdx ~= 1; h1.EdgeColor = [1 0 0]; end

                        title(titleStr)
                        if cBinIdx == 1 || cBinIdx == (NumBins/2)+1; 
                            ylabel('Probability');
                        else
                            ylabel('')
                            if configStruct.hideUnnecessaryTickLabels == 1
                                ax1.YTickLabel = ''; %ax1.TickLength = [0 0];
                            end
                        end
                        if cBinIdx > (NumBins/2); 
                            xlabel(xAxisLabel);
                        else 
                            xlabel('')
                            if configStruct.hideUnnecessaryTickLabels == 1
                                ax1.XTickLabel = ''; %ax1.TickLength = [0 0];
                            end
                        end
                        if cBinIdx == 1; 
                            legend(legendLabels, 'Location', 'northwest')
                        end
                    end
                    ax1.TickDir = 'out'; ax1.LineWidth = 2; h1.LineWidth = 2;
                    ax1.XTickLabelRotation = 45;
                    if ischar(configStruct.binWidth)
                        if strcmp(configStruct.binWidth, 'default') && configStruct.displayPeakAmplitude == 1
    %                         h1.BinWidth = 200;
                        elseif strcmp(configStruct.binWidth, 'default') && configStruct.displayDecayTimeConstant == 1
                            h1.BinWidth = 1;
                        end
                    elseif isnumeric(configStruct.binWidth) && h1.BinWidth ~= configStruct.binWidth
                        h1.BinWidth = configStruct.binWidth;
                    end


                    ax1 = gca; ax1.Box = 'off';
                    ax1.FontSize = fontSize;

                    if ~ischar(configStruct.xAxisMaxHistogram)
                        ax1.XLim = [configStruct.xAxisMinHistogram configStruct.xAxisMaxHistogram];
                    else
                        if strcmp(configStruct.xAxisMaxHistogram, '?') && ~ischar(configStruct.xAxisMinHistogram)
                            defaultValues = ax1.XLim;
                            ax1.XLim = [defaultValues(1) defaultValues(2)];
                        end
                    end 
                    if ~ischar(configStruct.yAxisMaxHistogram)
                        ax1.YLim = [configStruct.yAxisMinHistogram configStruct.yAxisMaxHistogram];
                    else
                        if strcmp(configStruct.yAxisMaxHistogram, '?') && ~ischar(configStruct.yAxisMinHistogram)
                            defaultValues = ax1.YLim;
                            ax1.YLim = [defaultValues(1) defaultValues(2)];
                        end
                    end                    
                end
            elseif notEnoughData == 1
                disp(strcat('Not enough data error: ned error1'))
            end
        end
    end
end