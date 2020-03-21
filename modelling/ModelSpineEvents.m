function ModelSpineEvents(varargin)

p = inputParser;

addParameter(p,'KeepROIs', 0)
%CG: Switch 'KeepROIs' to '1' if seeking to model effect of ROIs
%established from a previous condition. When KeepROIs = 1, the modelling
%script will generate model ROIs and enter them into randomly selected
%cells with a 3D cell array (the same kind of cell array that is generated
%by CoreOne and CoreTwo) - the Pre-Established SLM_CAs. This is useful for
%testing the ability of CoreThree and CoreFour to appropriately expand
%their cell arrays when a cell containing a ROI established by CoreOne or
%CoreTwo is encountered.

%CG: another important feature that will be activated if KeepROIs = 1 is
%that the modelling script will give priority to model events that are
%present before the first model frame is processed (all such events are
%entered into the SLM_CAs before the rest of the model data is generated.
%These events, called Distinguished Events in the modelling script, are
%intended to model baseline data and the ROIs established in the baseline
%data take precedence over the ideal ROIs of model events in a subsequent
%data group. Some relevant variables: DistinguishedFixedLocs,
%SLM_CA_ROIs_PreEst.
addParameter(p,'DefString', 'ModelDataTest')
addParameter(p,'SaveString', 'Save') %CG: 'Save' or 'NoSave'
addParameter(p,'LoadOrSaveIt', 'SaveIt') %CG: debugging feature. 
addParameter(p,'turnONimages', 0) 
addParameter(p,'EventNum', 4000) 
addParameter(p,'NumReps', 500) 
%CG: Total number of events (EventNum) and the number of events that will
%be RepeaterEvents(NumReps).
addParameter(p,'NumRepeaterLocs', 200) 
addParameter(p,'LargeEventSize',8)
addParameter(p,'StdEventSize_Max',3)
addParameter(p,'ROI_sl',3)
%CG: ROI_sl is the side length of the ROIs that will be modelled. 
addParameter(p,'ModelMatrixDim',[512 512 500])

parse(p,varargin{:}); param = p.Results;
KeepROIs = param.KeepROIs; DefString = param.DefString; SaveString = param.SaveString; 
LoadOrSaveIt = param.LoadOrSaveIt; turnONimages = param.turnONimages; EventNum = param.EventNum;
NumReps = param.NumReps; NumRepeaterLocs = param.NumRepeaterLocs;
LargeEventSize = param.LargeEventSize; StdEventSize_Max = param.StdEventSize_Max;
ROI_sl = param.ROI_sl; ModelMatrixDim = param.ModelMatrixDim;

SLM = zeros(ModelMatrixDim(1), ModelMatrixDim(2), ModelMatrixDim(3));

NumDistinguishedRepLocs = round(NumRepeaterLocs*(0.9));
%CG: Number of locations across the matrix that are available for the
%LargeEvents and RepeaterEvents to be located.
NumLargeEvents = NumRepeaterLocs; 
%CG: Every Repeater location will contain one and only one LargeEvent.
NeighborhoodSize = ROI_sl^2;

if KeepROIs == 1
    Sz_SLM = size(SLM);
    MaxNumPEComps = 100;
    MaxNumPEtimesteps = Sz_SLM(3);
    ZDim_PreEst = 2;
elseif KeepROIs == 0
    MaxNumPEComps = 0;
    MaxNumPEtimesteps = 0;
    ZDim_PreEst = 0;
end
%CG: MaxNumPEtimesteps, MaxNumPEComps and ZDim_PreEst determine the
%dimensions of the Pre-Established SLM_CAs. In real conditions the
%MaxNumPEComps (meaning: Maximum number of pre-established components')
%would be determined by the largest number of components (i.e. events) that
%occurred in any one frame. This means that CA_ROIs would have at least one
%row with that number of non-empty elements. There does not have to be at
%least one such row in the case of the SLM_CA_ROIs that is fed into the
%modelling script. For the modelling script to be still validating AARG it
%is sufficient that a scalable number of cells are non-empty in the CA_ROIs
%cell array and that the same cells in the CA_EventSize array each have a
%randomly generated event size. Similarly, with the ZDim_PreEst variable:
%this determines how many layers the Pre_Established SLM_CAs have. 

ROI_plusEZ_sl = ROI_sl*3;

%CG: An efficient way to ensure no unwanted touching or overlap between the
%ROIs is to create an exclusion zone around the ROI such that as long as
%the centroid of the current ROI does not overlap with the exclusion zone,
%there is no possibility that the two ROIs would be overlapping or
%touching. The exclusion zone is aimed at preventing the ROIs of FreeEvents
%erroneously overlapping with the ROIs of RepeaterEvents. Repeater events
%themselves will be randomly assigned a centroid value anywhere within the
%centroid neighbourhood of the designated fixed location centroid indices
%(i.e. one of the values in the FixedLocs array). This means that the
%exclusion zone must include any element around the centroid where the
%centroid neighbourhood of the RepeaterEvent COULD finally be. In addition
%to this, I re-wrote then modelling script slightly so that if a FreeEvent
%has to be reassigned to another location only the centroid index is
%checked against the ForbiddenIndices (see below), rather than the entire
%centroid neighbourhood. I suppose this would make the script a bit faster.

if KeepROIs == 1; DefString = strcat(DefString, '_KR'); end
if nargin < 1; LoadOrSaveIt = 'SaveIt'; end

if strcmp(LoadOrSaveIt, 'LoadIt')
    EventTypeRec = {};
    FEOrderCount_rec = [];
%CG: 'LoadIt' was used as a debugging feature.
    DesDir = uigetdir('', 'Select model to re-load');
    addpath(DesDir);
    cd(DesDir);

    load(strcat('RepeaterEvents_', DefString, '.mat'), 'RepeaterEvents')
%     load(strcat('REs2DLocs_', DefString, '.mat'), 'REs2DLocs')
    NumRepeaterEvents = numel(RepeaterEvents);
    load(strcat('FixedLocs_', DefString, '.mat'), 'FixedLocs');
    
    if KeepROIs == 0
        DistinguishedFixedLocs = [];
        SLM_CA_ROIs_PreEst = {}; SLM_CA_EventSize_PreEst = {}; PE_timestepList = [];
        PE_componentList = []; PreEst_EventSizeMtx = []; PreEst_cDFL = []; PreEst_cCell = [];
    end
    load(strcat('DistinguishedFixedLocs_', DefString, '.mat'), 'DistinguishedFixedLocs');
    load(strcat('PE_SLM_CAs_', DefString, '.mat'), 'SLM_CA_ROIs_PreEst', 'SLM_CA_EventSize_PreEst',...
      'PE_timestepList', 'PE_componentList', 'PreEst_EventSizeMtx', 'PreEst_cDFL', 'PreEst_cCell');
    
    NumRepeaterLocs = numel(FixedLocs);
    NumLargeEvents = NumRepeaterLocs; 
    load(strcat('LargeEvents_', DefString, '.mat'), 'LargeEvents');
    load(strcat('FE_Sizes_', DefString, '.mat'), 'LargeEventSize', 'StdEventSize_Max', 'NeighborhoodSize');
    [xdim, ydim, timesteps] = size(SLM);
  
    load(strcat('RNHoods_lgcMtx_', DefString, '.mat'), 'RepeaterNhoods_Mtx_log')
    
    load(strcat('FEOrderArray_', DefString,'.mat'), 'FEOrderArray_sc','FEOrderArray_fl');
    load(strcat('REOrderArray_', DefString,'.mat'), 'REOrderArray_sc','REOrderArray_fl');
    load(strcat('LEOrderArray_', DefString,'.mat'), 'LEOrderArray_sc','LEOrderArray_fl');
    FEOrderCounter_sc = 1;
    REOrderCounter_sc = 1;
    LEOrderCounter_sc = 1;
    load(strcat('RESizeList_', DefString,'.mat'), 'RESizeList_sc')
    load(strcat('FESizeList_', DefString,'.mat'), 'FESizeList_sc')
    RESizeCounter_sc = 1;
    FESizeCounter_sc = 1;

    load(strcat('FCLocs_sorted_', DefString,'.mat'), 'FCLocs_sorted');
    
%%%%%%%testing
    load(strcat('SLM_', DefString,'.mat'), 'SLM')
    SLM_done = SLM; SLM = SLM.*0;
%%%%%%%testing

elseif strcmp(LoadOrSaveIt, 'SaveIt')   
%CG: When the user wish to generate new or rewrite old modelling data. 

    DesDir = uigetdir('', 'Save model files to (just select parent folder. Folder containing model data created automatically...');

    mkdir(DesDir, DefString)
    addpath(strcat(DesDir, '/', DefString))
    
    cd(strcat(DesDir, '/', DefString));

    save(strcat('SLM_', DefString, '.mat'), 'SLM');
    [xdim, ydim, timesteps] = size(SLM);
    RepeaterEvents = randperm(EventNum, NumReps);
    save(strcat('RepeaterEvents_', DefString, '.mat'), 'RepeaterEvents')

    NumRepeaterEvents = numel(RepeaterEvents);
    FixedLocs = randperm(xdim*ydim, NumRepeaterLocs);
    LargeEvents = randsample(RepeaterEvents, NumLargeEvents);

    save(strcat('NumLargeEvents_', DefString, '.mat'), 'NumLargeEvents');
    save(strcat('LargeEvents_', DefString, '.mat'), 'LargeEvents');
    save(strcat('FE_Sizes_', DefString, '.mat'), 'LargeEventSize', 'StdEventSize_Max', 'NeighborhoodSize');

    LEOrderArray_sc = zeros(1,NumLargeEvents); LEOrderArray_fl = zeros(1,NumLargeEvents);
    REOrderArray_sc = zeros(1,NumRepeaterEvents); REOrderArray_fl = zeros(1,NumRepeaterEvents);
    FEOrderArray_sc = zeros(1,EventNum); FEOrderArray_fl = zeros(1,EventNum);
    
    FEOrderCounter_sc = 1;
    REOrderCounter_sc = 1;
    LEOrderCounter_sc = 1;
%CG: the exact centroid locations for each event need to be saved when strcmp(LoadOrSaveIt, 'SaveIt')  
%(including centroids of free events) otherwise events can end up touching
%or overlapping if the data is loading during code development. 
    RESizeList_sc = zeros(1,NumRepeaterEvents);
    FESizeList_sc = zeros(1,(EventNum));
    RESizeCounter_sc = 1;
    FESizeCounter_sc = 1;
%CG: the OrderArrays and SizeLists record the order in which locations were
%selected from the FixedLocs array - in the case of LargeEvents and
%RepeaterEvents. Sizes of the RepeaterEvents and FreeEvents are also stored
%when the 'SpecialSave' argument is called. This is the essential task that
%'SpecialSave' does. The size of Non-LargeEvents and the locations of
%LargeEvents and RepeaterEvents are randomly selected, which would make it
%otherwise impossible to reproduce exactly the same model. FreeEvent
%locations also vary randomly, but this is taken care of later.

end

dispstat('','init')

AccessCtrl = 1;
AccessCtrl2 = 1;
AccessCtrl3 = 1;
ErrorCount = 1;
ErrorCount3 = 1;

%SLM = Synthetic Logical Matrix; This matrix will contain circular 2D
%events of a random size within an expected range (to represent detected
%calcium elevations in spines) distributed across random (FCLocs) and
%non-random (FixedLocs) locations in the matrix.
if NumRepeaterLocs < NumLargeEvents
    errordlg('Error: The number of large events should not exceed the number of Repeater (fixed) locations');
    BreakVariable = NonExistentVariable;
end

LELimit = round(2*(NumLargeEvents/timesteps));
if LELimit == 0
    LELimit = 1;
end

%CG: LELimit sets the threshold for how many LEs can occur in any given
%timestep. This seems to help limit the number of times the modelling
%script needs to be reset. 

SLM_Compiled = zeros(xdim, ydim);
SLM_TouchPrevention = zeros(xdim, ydim);
SLM_TouchPrevention_ii = zeros(xdim, ydim);

%CG: in addition to ensuring that there are no overlapping events (in terms
%event dimensions) it is also important to ensure that the ROIs of Free
%Events do not overlap spatially or temporally 
SLM_ROI_TP = zeros(xdim, ydim);
SLM_ROI_TP_ii = zeros(xdim, ydim);

FECounterRecord = zeros(1, (timesteps-1));

NhoodSize_plusEZ = ROI_plusEZ_sl^2;
%CG: '_plusEZ' means 'plus exclusion zone'. The exclusion zone is a
%ROI_sl-element thick edge around the ROI_sl-by-ROI_sl ROI that should also
%be included in the ForbiddenIndices array in order to avoid Free Events
%touching (and therefore overlapping with) the established ROIs of Repeater
%Events.
ForbiddenIndices = zeros(numel(NumRepeaterLocs)*NhoodSize_plusEZ, 1);
%CG: It is important to ensure that no repeater event centroids are too
%close to the edge of the matrix. If the distance of the centroid to the
%closest edge of the matrix is less than the radius of the LargeEvents then
%a LargeEvent in the modelling script cell array for EventSize will have a
%different size compared to the EventSize cell array produced by
%AARG

NRL_Counter = 0;
if NumRepeaterLocs <= 1

    NRL_Counter = 1;

end

CNhood_rad = ceil(sqrt(NeighborhoodSize)/2);

%CG: AARG generates a Window Matrix and overlays the centroid of the ideal
%ROI of the current event with the centre point of the Window Matrix.
%Determining how close to the edge any event can be depends on if the
%diameter of the LargeEventSize or the side length of the Window Matrix is
%larger. LargeEventSize+CNhood_rad needs to be used because when a
%LargeEvent or RepeaterEvent is placed at a FixedLoc location, the exact
%centroid coordinates of this model event is determined randomly from a
%list of all CNhood values for the current FixedLoc index. This makes the
%script model synaptic events more closely.

WinMtx_sl_hl = ceil((3*ROI_sl)/2);
if strcmp(LoadOrSaveIt, 'SaveIt') 

    if (LargeEventSize+CNhood_rad) >= WinMtx_sl_hl
%CG: MinRval, MaxRval, etc: these variables are used to ensure that
%RepeaterEvents are not too close to the edge of the model.
        MinRval = (LargeEventSize + CNhood_rad);
        MinCval = (LargeEventSize + CNhood_rad);
        MaxRval = (xdim) - (LargeEventSize + CNhood_rad);
        MaxCval = (ydim) - (LargeEventSize + CNhood_rad);

    elseif (LargeEventSize+CNhood_rad) < WinMtx_sl_hl

        MinRval = (WinMtx_sl_hl);
        MinCval = (WinMtx_sl_hl);
        MaxRval = (xdim) - (WinMtx_sl_hl);
        MaxCval = (ydim) - (WinMtx_sl_hl);

    end

    RepeaterNhoods_plusEZ_Mtx = zeros(xdim, ydim);
    RepeaterNhoods_plusEZ_MemMtx = RepeaterNhoods_plusEZ_Mtx;
    RepeaterNhoods_Mtx = zeros(xdim, ydim);
    RepeaterNhoods_MemMtx = RepeaterNhoods_Mtx;
    RepeaterARL = 0;
    Error5Flag = 0;

    for RepeaterLoc = 1 : NumRepeaterLocs

%CG: The current for-loop finds all the CNhood indices for the Repeater
%centroids. None of the FCLocs indices should overlap with these.

%CG: It is important to ensure that no repeater event centroids are too
%close to the edge of the matrix. If the distance of the centroid to the
%closest edge of the matrix is less than the radius of the LargeEvents then
%a LargeEvent in the SLM cell array for EventSize will have a
%different size compared to the EventSize cell array produced by
%AARG
        PosiCtrl = 1;
        while PosiCtrl == 1 && RepeaterARL < 10000

            RepeaterARL = RepeaterARL + 1; dispstat(num2str(RepeaterARL))
            [R_RL, C_RL, ~] = ind2sub(size(SLM), FixedLocs(RepeaterLoc)); 

            if R_RL >= MinRval && R_RL <= MaxRval && C_RL >= MinCval && C_RL <= MaxCval    

                Row_CNhoodDiag_plusEZ = R_RL-(floor(ROI_sl/2)+ROI_sl):1:R_RL+(floor(ROI_sl/2)+ROI_sl);
                Col_CNhoodDiag_plusEZ = C_RL-(floor(ROI_sl/2)+ROI_sl):1:C_RL+(floor(ROI_sl/2)+ROI_sl);

                Row_CNhoodDiag = R_RL-floor(ROI_sl/2):1:R_RL+floor(ROI_sl/2);
                Col_CNhoodDiag = C_RL-floor(ROI_sl/2):1:C_RL+floor(ROI_sl/2);
%CG: Row_CNhood gives the row values for the centroid neighborhood.
                for cEl = 1 : ROI_sl
                    Row_CNhood(cEl:ROI_sl:ROI_sl^2) = Row_CNhoodDiag(cEl);
                end
%CG: Row_CNhood_plusEZ gives the row values for the centroid neighborhood plus the
%(ROI_sl)-element thick exclusion zone surrounding the CNhood. 
                for cEl = 1 : ROI_plusEZ_sl
                    Row_CNhood_plusEZ(cEl:ROI_plusEZ_sl:ROI_plusEZ_sl^2) = Row_CNhoodDiag_plusEZ(cEl);
                end 

                for cEl = 1 : ROI_sl
                    if cEl == 1
                        Col_CNhood(cEl:1:ROI_sl) = Col_CNhoodDiag(cEl);
                    elseif cEl > 1
                        Col_CNhood(((cEl-1)*ROI_sl)+1:1:ROI_sl*cEl) = Col_CNhoodDiag(cEl);
                    end
                end

                for cEl = 1 : ROI_plusEZ_sl
                    if cEl == 1
                        Col_CNhood_plusEZ(cEl:1:ROI_plusEZ_sl) = Col_CNhoodDiag_plusEZ(cEl);
                    elseif cEl > 1
                        Col_CNhood_plusEZ(((cEl-1)*ROI_plusEZ_sl)+1:1:ROI_plusEZ_sl*cEl) = Col_CNhoodDiag_plusEZ(cEl);
                    end
                end

                Ind_RepeaterNhood_plusEZ = sub2ind(size(RepeaterNhoods_plusEZ_Mtx), Row_CNhood_plusEZ, Col_CNhood_plusEZ);
                Ind_RepeaterNhood = sub2ind(size(RepeaterNhoods_Mtx), Row_CNhood, Col_CNhood);

                RepeaterNhoods_plusEZ_Mtx(Ind_RepeaterNhood_plusEZ) = 1;
                RepeaterNhoods_Mtx(Ind_RepeaterNhood) = 1;

                CC_plusEZ = bwconncomp(RepeaterNhoods_plusEZ_Mtx);
                NumRepeaterComps = CC_plusEZ.NumObjects;
%CG: if the NumRepeaterComps does not match the value of RepeaterLoc, this
%could mean that the events are overlapping. This is more likely for more
%complex matrices and would cause problems. 
                if NumRepeaterComps == RepeaterLoc
                    RepeaterNhoods_plusEZ_MemMtx = RepeaterNhoods_plusEZ_Mtx;
                    RepeaterNhoods_MemMtx = RepeaterNhoods_Mtx;
                    RepeaterARL = 0;
                    PosiCtrl = 0;
                elseif NumRepeaterComps ~= RepeaterLoc
                    RepeaterNhoods_plusEZ_Mtx = RepeaterNhoods_plusEZ_MemMtx;
                    RepeaterNhoods_Mtx = RepeaterNhoods_MemMtx;
                    FixedLocs(RepeaterLoc) = randperm(xdim*ydim, 1);
                end
            else 
                FixedLocs(RepeaterLoc) = randperm(xdim*ydim, 1);
            end
        end

        if RepeaterARL >= 10000
            Error5Flag = 1;
            break 
        end

        StartIdx = (RepeaterLoc*NhoodSize_plusEZ) - NhoodSize_plusEZ + 1;
        EndIdx = (RepeaterLoc*NhoodSize_plusEZ);
        ForbiddenIndices(StartIdx : EndIdx) = Ind_RepeaterNhood_plusEZ;    

    end

    if Error5Flag == 1
        errordlg(strcat('Error5: A valid repeater location could not be found',...
            '. Try reducing the number of repeater locations to avoid this error'));
        BreakVariable = NonExistentVariable;
    end

    RepeaterNhoods_Mtx_log = logical(RepeaterNhoods_Mtx);
    save(strcat('RNHoods_lgcMtx_', DefString, '.mat'), 'RepeaterNhoods_Mtx_log')

end

SLM_ROI_TP(:) = 0;
SLM_ROI_TP(RepeaterNhoods_Mtx_log) = 1;
SLM_ROI_TP_ii(:) = 0;
SLM_ROI_TP_ii(RepeaterNhoods_Mtx_log) = 1;

FixedLocs_4LEs = FixedLocs;
% FixedLocs_Mem = FixedLocs_4LEs;

if strcmp(LoadOrSaveIt, 'SaveIt') && KeepROIs == 1
    DistinguishedFixedLocs = randsample(FixedLocs, NumDistinguishedRepLocs);
    DistinguishedFixedLocs_copy = DistinguishedFixedLocs;
%CG: DistinguishedFixedLocs represent the ROIs established from previous
%conditions. If the modelling script has been written correctly, then
%this function should not overwrite these ROIs, even if the current event
%is found to be a LargeEvent. 

        SLM_CA_ROIs_PreEst = cell(MaxNumPEtimesteps, MaxNumPEComps, ZDim_PreEst);
        SLM_CA_EventSize_PreEst = cell(MaxNumPEtimesteps, MaxNumPEComps, ZDim_PreEst);

        TotalPreEstCellNum = MaxNumPEtimesteps*MaxNumPEComps*ZDim_PreEst;
        NonEmpCells = randperm(TotalPreEstCellNum, NumDistinguishedRepLocs);

        PreEst_cDFL = zeros(NumDistinguishedRepLocs, 1);
        PreEst_cCell = zeros(NumDistinguishedRepLocs, 1);
%CG: PreEst_ROI_idx stores the DistinguishedFixedLoc value in the first
%column and in the second column, the respective cell index in which it was
%used when SLM_CA_ROIs_PreEst was filled in. 

        PE_timestepList = [];
        PE_componentList = [];

        PreEstCounter = 0;

        PreEst_EventSizeMtx = zeros(xdim, ydim);

        for cCell = 1 : TotalPreEstCellNum

            if ~isempty(find(NonEmpCells == cCell,1))

                PreEstCounter = PreEstCounter + 1;

                if numel(DistinguishedFixedLocs_copy) > 1
                    cDFL = randsample(DistinguishedFixedLocs_copy, 1);
                else
                    cDFL = DistinguishedFixedLocs_copy;
                end

                lgcIdx = DistinguishedFixedLocs_copy == cDFL;
                DistinguishedFixedLocs_copy(lgcIdx) = [];

                PreEst_cDFL(PreEstCounter, 1) = cDFL;
                PreEst_cCell(PreEstCounter, 1) = cCell;

                [R_RL, C_RL, Z_RL] = ind2sub(size(SLM), cDFL);
                Row_CNhoodDiag = R_RL-floor(ROI_sl/2):1:R_RL+floor(ROI_sl/2);
                Col_CNhoodDiag = C_RL-floor(ROI_sl/2):1:C_RL+floor(ROI_sl/2);

                for cEl = 1 : ROI_sl
                    Row_CNhood(cEl:ROI_sl:ROI_sl^2) = Row_CNhoodDiag(cEl);
                end
                for cEl = 1 : ROI_sl
                    if cEl == 1
                        Col_CNhood(cEl:1:ROI_sl) = Col_CNhoodDiag(cEl);
                    elseif cEl > 1
                        Col_CNhood(((cEl-1)*ROI_sl)+1:1:ROI_sl*cEl) = Col_CNhoodDiag(cEl);
                    end
                end

                PE_CNhood = sub2ind(size(SLM(:,:,1)), Row_CNhood, Col_CNhood);

                SLM_CA_ROIs_PreEst{cCell} = PE_CNhood;
                sesize_PreEst = randperm(StdEventSize_Max, 1);

                se = strel('disk', sesize_PreEst, 0);

                FEM = se.getnhood;
%CG: FEM = Fake Event Matrix 
                [xdim_se, ydim_se] = size(FEM);

                T = round(((xdim_se)/2));
                S = round(((ydim_se)/2));

                CInd_FEM = sub2ind(size(FEM), T, S);
%CG: CInd_FEM is the centroid location of the fake event matrix (FEM).
                AllPosInds_FEM = find(FEM);

                [RowVec_FEM, ColVec_FEM, ZVec_FEM] = ind2sub(size(FEM), AllPosInds_FEM);
                [RowVec_CIndFEM, ColVec_CIndFEM, ZVec_CIndFEM] = ind2sub(size(FEM), CInd_FEM);

                SubtVec_R = ones(size(RowVec_FEM)).*RowVec_CIndFEM;
                SubtVec_C = ones(size(ColVec_FEM)).*ColVec_CIndFEM;
                SubtVec_Z = ones(size(ZVec_FEM)).*ZVec_CIndFEM;

                PosiVec_R = RowVec_FEM - SubtVec_R;
                PosiVec_C = ColVec_FEM - SubtVec_C;
                PosiVec_Z = ZVec_FEM - SubtVec_Z;

                TempVec_RowDLocs = ones(size(RowVec_FEM)).*R_RL;
                TempVec_ColDLocs = ones(size(ColVec_FEM)).*C_RL;
                TempVec_ZDLocs = ones(size(ZVec_FEM)).*Z_RL;

                FinalRowVector = TempVec_RowDLocs + PosiVec_R;
                FinalColVector = TempVec_ColDLocs + PosiVec_C;
                FinalZVector = TempVec_ZDLocs + PosiVec_Z;

                Ind2Convert_PE = sub2ind(size(SLM), FinalRowVector(:), FinalColVector(:), FinalZVector(:));
                Ind2Convert_PE = Ind2Convert_PE';
               
                SLM_CA_EventSize_PreEst{cCell} = Ind2Convert_PE;
                PreEst_EventSizeMtx(Ind2Convert_PE) = 1;

                [cTimestep, cComp] = ind2sub(size(SLM_CA_ROIs_PreEst), cCell);
                if isempty(PE_timestepList)
                    PE_timestepList = cTimestep;
                    PE_componentList = cComp;
                else
                    PE_timestepList(end+1) = cTimestep;
                    PE_componentList(end+1) = cComp;
                end   
            end
        end
elseif strcmp(LoadOrSaveIt,'SaveIt') && KeepROIs == 0   
        DistinguishedFixedLocs = [];
        SLM_CA_ROIs_PreEst = {}; SLM_CA_EventSize_PreEst = {}; PE_timestepList = [];
        PE_componentList = []; PreEst_EventSizeMtx = []; PreEst_cDFL = []; PreEst_cCell = [];   
end

if ~strcmp(LoadOrSaveIt, 'LoadIt')
    save(strcat('DistinguishedFixedLocs_', DefString, '.mat'), 'DistinguishedFixedLocs');
    save(strcat('PE_SLM_CAs_', DefString, '.mat'), 'SLM_CA_ROIs_PreEst', 'SLM_CA_EventSize_PreEst',...
        'PE_timestepList', 'PE_componentList', 'PreEst_EventSizeMtx', 'PreEst_cDFL', 'PreEst_cCell');
end
if strcmp(LoadOrSaveIt, 'SaveIt')
    save(strcat('FixedLocs_', DefString, '.mat'), 'FixedLocs');
end

%CG: The ForbiddenIndices list originally incorporated only the Repeater
%event locations, but it is expanded to include any index that is
%within one StdEventSize_Max+CNhood_rad or WinMtx_sl_hl from any edge.
%Without this there would be a risk that standard events are recorded as
%having different sizes by the Model and AARG. WinMtx_sl_hl
%refers to the the half-length of the window matrix. The window matrix is
%an array created in the AARG_spines function. It looks at the relevant
%section of OnsetMatrix (see AARG Cores) in such a way that ROIs that
%appear smaller than their actual size in the window matrix will be
%non-overlapping events and all touching or overlapping events will appear
%to have the expected size. The central element of the window matrix
%overlays the centroid of the current event. The window matrix doesn't show
%the ideal ROI of the current event, only the ROIs that would overlap with
%or touch it. AARG automatically ignores events that lie closer than
%WinMtx_sl_hl to any of the 4 borders. This has to happen otherwise
%AARG would not be able to form the window matrix.

%CG: Indices that are too close to the matrix edge (specifically for
%FreeEvents) are incorporated into the ForbiddenIndices (with the
%MoreFIdces variable.

if WinMtx_sl_hl <= (StdEventSize_Max+CNhood_rad)

    SLMcopy = logical(SLM(:,:,1));
    SLMcopy(:) = 1;
    SLMcopy((StdEventSize_Max+CNhood_rad)+1 : xdim - (StdEventSize_Max+CNhood_rad),...
        (StdEventSize_Max+CNhood_rad)+1 : ydim - (StdEventSize_Max+CNhood_rad)) = 0;

elseif WinMtx_sl_hl > (StdEventSize_Max+CNhood_rad)

    SLMcopy = logical(SLM(:,:,1));
    SLMcopy(:) = 1;
    SLMcopy(WinMtx_sl_hl : (xdim - (WinMtx_sl_hl + 1)), WinMtx_sl_hl : (ydim - (WinMtx_sl_hl + 1))) = 0;

end

MoreFIdces = find(SLMcopy);
ForbiddenIndices = cat(1, ForbiddenIndices, MoreFIdces);    
LastIndex = sub2ind(size(SLM), xdim, ydim, timesteps);

if strcmp(LoadOrSaveIt, 'SaveIt')

    FCLocs = randperm(LastIndex, EventNum);
    Error4Flag = 0;
    Sz_FCLocs = size(FCLocs);
    NumFCLocs = Sz_FCLocs(1)*Sz_FCLocs(2);

    SLM2ndCopy = logical(SLMcopy);
    SLM2ndCopy(MoreFIdces) = 0;

    for CurrentEl = 1 : NumFCLocs    
%CG: The current for-loop ensures that none of the FCLocs are equal to any
%of the CNhoods for the FixedLocs.

        [Row_FCLoc, Col_FCLoc, ~] = ind2sub(size(SLM), FCLocs(CurrentEl));  
        FCLocs2D = sub2ind(size(SLM(:, :, 1)), Row_FCLoc, Col_FCLoc);

        if ~isempty(find(ForbiddenIndices(:) == FCLocs2D, 1))

            CleanFCLocs = 0;
            RepeaterAV = 0;
            while CleanFCLocs == 0 && Error4Flag == 0


                RepeaterAV = RepeaterAV + 1;
                AlterValue = randperm(LastIndex, 1);
                [AlterRow_FCLoc, AlterCol_FCLoc, ~] = ind2sub(size(SLM), AlterValue);
                AlterValue2D = sub2ind(size(SLM(:, :, 1)), AlterRow_FCLoc, AlterCol_FCLoc);

                if ~isempty(find(ForbiddenIndices(:) == AlterValue2D, 1))
%CG: the ForbiddenIndices list contains all the indices that the centroid
%of the current event would have to avoid in order to prevent having its
%CNhood overlapping or touching the CNhood of a previously established
%event. 
                    CleanFCLocs = 0;
                    RepeaterAV = 0;
                elseif isempty(find(ForbiddenIndices(:) == AlterValue2D, 1))
                    FCLocs(CurrentEl) = AlterValue;
                    SLM2ndCopy(AlterValue2D) = 1;
                    CleanFCLocs = 1;
                end   

                if RepeaterAV >= 10000
                    Error4Flag = 1;
                    break
                end
            end
        end

        if Error4Flag == 1
            errordlg('Error4: 200,000 attempts were made to relocate one of the Model Event centroids. It is likely that an acceptable location will not be found with the settings you are using. Try increasing the time dimension of your model matrix or decreasing the number of model events');
            break
            BreakVariable = NonExistentVariable;
        end
    end    

    FCLocs_sorted = sort(FCLocs);
%CG: FCLocs contain a list of all the centroid indices that could be used
%as a central point for a FreeEvent. Sorting this array is one of the steps
%that is necessary to ensure that the appropriate index value is given for
%the current frame (or 'timestep'). 
end

EventNum4TS = 0;
timestepsDone = 0;
EventsDone = 0;

PosLELocs = zeros(NumRepeaterEvents, ROI_sl^2);
DefLELocs = zeros(NumLargeEvents, ROI_sl^2);
MasterFEs = zeros(NumLargeEvents, 1);
timestepMFE = zeros(NumLargeEvents, 1);
MasterTally = zeros(NumLargeEvents, 1);
zMFE = zeros(NumLargeEvents, 1);

RepeaterFEs = zeros(NumRepeaterEvents, 1);
timestepRFE = zeros(NumRepeaterEvents, 1);
zRFE = zeros(NumRepeaterEvents, 1);

if KeepROIs == 0
    SLM_CA_ROIs = {};
    SLM_CA_EventSize = {};
    SLM_CA_EventNum = {};
    SLM_CA_FrameIEI = {};
elseif KeepROIs == 1
    SLM_CA_ROIs = SLM_CA_ROIs_PreEst;
    SLM_CA_EventSize = SLM_CA_EventSize_PreEst;
    SLM_CA_EventNum = cell(MaxNumPEtimesteps, MaxNumPEComps, ZDim_PreEst);
    SLM_CA_FrameIEI = cell(MaxNumPEtimesteps, MaxNumPEComps, ZDim_PreEst);
    SLM_CA_EventType = cell(MaxNumPEtimesteps, MaxNumPEComps, ZDim_PreEst);
end

counterLE = 0;
counterNLE = 0;
LargeEventCounter = 0;

CurrentFCLocsIdx = 0;

RepeaterE = 0;
RepeaterF = 0;

Error2Flag = 0;
Error3Flag = 0;

NewNumLargeEvents = NumLargeEvents;

CheckingFEOverlaps = 0;

CA_Debug = {};

CA_Debug{1,1} = 'timestep';
CA_Debug{1,2} = 'NumLEs';
CA_Debug{1,3} = 'NumREs';
CA_Debug{1,4} = 'NumFEs';
CA_Debug{1,5} = 'TotalEventNum';
CA_Debug{1,6} = 'FixedLocsTaken';

wb = waitbar(0, strcat('Generating model: "', DefString, '"'));
wbStep = floor(timesteps/100);

% CompDetMtx = zeros(Sz_SLM(1), Sz_SLM(2));

%%%%testing 
% fh = figure;
% fh_new = figure('Name', 'New', 'NumberTitle','Off'); fh_done = figure('Name', 'Done', 'NumberTitle','Off');
%%%%testing
REOrderCounter_scList = [];

for timestep = 1 : timesteps

    RepeaterK = 0;
    RepeaterC = 0;

    IdxStart = (xdim*ydim)*(timestep - 1) + 1;
    IdxEnd = (xdim*ydim)*timestep;

    if Error2Flag == 1 || Error3Flag == 1        
        break       
    end

    EventsDone = EventsDone + EventNum4TS;    
    timestepsDone = timestepsDone + 1;

    if timestep == 1
        CurrentFCLocsRange = (1 : xdim*ydim);
    elseif timestep > 1    
        CurrentFCLocsRange = (IdxStart : IdxEnd);
    end 

    LastIdx = CurrentFCLocsRange(end);

    if timestep == 1
        CurrentFCLocsIdx = 0;
        ElValue_FCLocs = FCLocs_sorted(1);
    end
    LastIdxCtrl = 0;

%CG: this while-loop finds the fake centroid locations for the current
%timestep (Current_FCLocs variable). The values within this array may be
%altered at a later stage if inappropriate overlap is found with another
%event. 

    FELimitCtrl = 1;
  
    while FELimitCtrl == 1 && CurrentFCLocsIdx <= EventNum - 1

        RepeaterE = RepeaterE + 1;
        CurrentFCLocsIdx = CurrentFCLocsIdx + 1;
        LastIdxCtrl = LastIdxCtrl + 1;

        if LastIdxCtrl == 1
            Ind_RightOnes = zeros(1);
            Ind_RightOnes(1) = CurrentFCLocsIdx;
        elseif LastIdxCtrl > 1
            Ind_RightOnes(end+1) = CurrentFCLocsIdx;
        end

        if CurrentFCLocsIdx+1 > EventNum
            FELimitCtrl = 0;
        else
            ElValue_FCLocs = FCLocs_sorted(CurrentFCLocsIdx);
            ElValue_FCLocs_List(RepeaterE) = FCLocs_sorted(CurrentFCLocsIdx);
        end

        if ElValue_FCLocs > LastIdx
            FELimitCtrl = 0;
            Ind_RightOnes(end) = [];
            CurrentFCLocsIdx = CurrentFCLocsIdx - 1;
            RepeaterE = RepeaterE - 1;
        end
    end

    Current_FCLocs = FCLocs_sorted(Ind_RightOnes);

    if timestep == 1

        FakeEventStart = 1;
        FakeEventEnd = numel(Current_FCLocs);

    elseif timestep > 1

        FakeEventStart = FakeEventEnd + 1;
        FakeEventEnd = (numel(Current_FCLocs) + FakeEventStart) - 1;

    end

    FakeEventCounter = 0;
    FCLocsReassign = 0;
    RepeaterEventCounter = 0;

    if timestep > 1
        FECounterRecord(timestep-1) = FreeEventCounter;
    end
    FreeEventCounter = 0;

    SLM_TouchPrevention(:) = 0;
    SLM_TouchPrevention_ii(:) = 0;

    RECtrl_Counter = 0;

    Str_FCLocs = struct([]);
%CG: Empty structure array to hold data for sorting at the end of each
%FakeEvent loop. Sorting happens prior to entering the data into the cell
%arrays. 

    FixedLocs_TakenThisRound = [];
    if FakeEventEnd >= FakeEventStart && FakeEventStart <= EventNum
        for FakeEvent = FakeEventStart : FakeEventEnd
            
            if FakeEventStart > FakeEventEnd 
                SLM_CA_ROIs(timestep) = {[]};
                SLM_CA_EventNum(timestep) = {[]};
                SLM_CA_EventSize(timestep) = {[]};
                SLM_CA_FrameIEI(timestep) = {[]};

            elseif FakeEventStart ~= FakeEventEnd
     
                if Error3Flag == 1 
                   break
                end
                generalRep = 0;
                FakeEventCounter = FakeEventCounter + 1;

                FakeEvents_ThisTS = (FakeEventStart : FakeEventEnd);

                if ~isempty(find(RepeaterEvents(:) == FakeEvents_ThisTS(FakeEventCounter), 1)) 
                    RECtrl_Counter = RECtrl_Counter + 1;
                end

                if NumRepeaterLocs < RECtrl_Counter

                    errordlg('Error2: In one frame, the number of Repeater Events exceeded the NumRepeaterLocs. Try calling this function again. To decrease the chances of this error re-occurring, increase the value for NumRepeaterLocs');

                    Error2Flag = 1;
                    BreakVariable = NonExistentVariable;
                    break 

                end

                if ~isempty(find(RepeaterEvents(:) == FakeEvent,1)) 
                    RepeaterEventCounter = RepeaterEventCounter + 1;
                elseif isempty(find(RepeaterEvents(:) == FakeEvent,1))
                    FreeEventCounter = FreeEventCounter + 1;
                end
                
                if FakeEvent == FakeEventStart 

                    for CurrentLE = 1 : NewNumLargeEvents
                        Ind_SingleLE = find(FakeEvents_ThisTS(:) == LargeEvents(CurrentLE));

                        if ~isempty(Ind_SingleLE) && CurrentLE == 1
                            Ind_LEsThisTS = Ind_SingleLE;
                        elseif ~isempty(Ind_SingleLE) && CurrentLE > 1 && ~exist('Ind_LEsThisTS', 'var')
                            Ind_LEsThisTS = Ind_SingleLE;
                        elseif ~isempty(Ind_SingleLE) && CurrentLE == 1 && exist('Ind_LEsThisTS', 'var')
                            Ind_LEsThisTS = Ind_SingleLE;
                        elseif ~isempty(Ind_SingleLE) && CurrentLE > 1 && exist('Ind_LEsThisTS', 'var')
                            Ind_LEsThisTS(end+1) = Ind_SingleLE;
                        elseif isempty(Ind_SingleLE)
                            Ind_LEsThisTS = Ind_SingleLE;
                        end
                    end
                end
                NumAttempts = 0; 
                TouchingEvents = 1;
                

                while TouchingEvents > 0

                    TouchingEvents = 0; 
%CG: TouchingEvents is an important variable that switches off the while
%loop if all events in the current frame (timestep) are not overlapping or
%touching. If two or more events in the current frame are overlapping or
%touching, AARG would treat them as a single event, but the modelling
%script would still record them as two distinct events. The user can set a
%limit on how many times the function can loop at this point using the
%variable NumAttempts. Depending on the inputs the user tries to set, it
%can be that the function loops indefinitely, if no limit is put in place.
%CG: For some reason, Matlab has no straight forward abort function, so I
%use this 'controlled error' to trigger a forced abort of the function
%call. 
                    NumAttempts = NumAttempts + 1;
                    generalRep = generalRep + 1;

                    if NumAttempts > 1000 && ~isempty(find(RepeaterEvents(:) == FakeEvent,1)) && isempty(find(LargeEvents(:) == FakeEvent,1))
                        errordlg('a place for the current (Repeater) Event cannot be found');
                        BreakVariable = NonExistentVariable;
                        break
                    elseif NumAttempts > 1000 && ~isempty(find(LargeEvents(:) == FakeEvent,1))
                        errordlg('a place for the current (Large) Event cannot be found');
                        BreakVariable = NonExistentVariable;
                        break
                    elseif NumAttempts > 1000 && isempty(find(RepeaterEvents(:) == FakeEvent,1))
                        errordlg('a place for the current (Free) Event cannot be found');
                        BreakVariable = NonExistentVariable;
                        break
                    end

                    if exist('Ind_LEsThisTS', 'var') && FakeEvent == FakeEventStart

                        NumLEs_ThisTS = numel(FakeEvents_ThisTS(Ind_LEsThisTS));

%CG: The number of LEs in this timestep, should not exceed the number set
%by the scalar variable 'LELimit'. LELimit helps to distribute the
%LargeEvents more evenly across the model and reduce the number of
%controlled errors. 
                        if NumLEs_ThisTS > LELimit
                            LargeEvents_ThisTS = FakeEvents_ThisTS(Ind_LEsThisTS);
                            ExcessLEs = NumLEs_ThisTS - LELimit;

                            LEs2Remove = randsample(LargeEvents_ThisTS, ExcessLEs);
                            NumLEs2Remove = numel(LEs2Remove);

                            for Remover = 1 : NumLEs2Remove
                                Ind_EventRemover = find(LargeEvents_ThisTS == LEs2Remove(Remover));

                                if ~isempty(Ind_EventRemover) && Remover == 1
                                    Ind_EventRemoverSum = Ind_EventRemover;
                                elseif ~isempty(Ind_EventRemover) && Remover > 1 && ~exist('Ind_EventRemover', 'var')
                                    Ind_EventRemoverSum = Ind_EventRemover;                    
                                elseif ~isempty(Ind_EventRemover) && Remover > 1 && exist('Ind_EventRemover', 'var')
                                    Ind_EventRemoverSum(end+1) = Ind_EventRemover;
                                end

                            end

                            LargeEvents_ThisTS(Ind_EventRemoverSum) = [];

                            Ind_RestrLargeEvents = LargeEvents > FakeEventEnd;
                            RestrLargeEvents = LargeEvents(Ind_RestrLargeEvents);

                            Ind_RestrRepeaterEvents = RepeaterEvents > FakeEventEnd;
                            RestrRepeaterEvents = RepeaterEvents(Ind_RestrRepeaterEvents);

                            LEsUnique = 1;
                            while LEsUnique == 1
                                LEsUnique = 0;
                                LargeEvents2Add = randsample(RestrRepeaterEvents, NumLEs2Remove);

                                LargeEvents = horzcat(RestrLargeEvents, LargeEvents2Add, LargeEvents_ThisTS);

                                CheckUniqueness = unique(LargeEvents);
                                TK = numel(CheckUniqueness);
                                UL = numel(LargeEvents);
                                if TK ~= UL
                                    LEsUnique = 1;
                                end
                            end

                            NewNumLargeEvents = numel(LargeEvents);
                        end
%CG: 'LargeEvents' now excludes all the Repeater events below FakeEventEnd
%(that is the last FakeEvent to appear in the current timestep), with the
%exception of LELimit number of LEs that fall into the current timestep.
                    end

                    RepeaterF = RepeaterF + 1;

                    if Error3Flag == 1
                       break
                    end

%CG: If, at the end of this while-loop, the current FakeEvent is found to
%overlap with or touch another event in the same frame, it is necessary to
%reassign the centroid index for this FakeEvent. It is also important to
%ensure that this newly assigned index does not overlap with any of the
%CNhood values for the FixedLocs. The condition statement 'if
%FCLocsReassign == 1' and while-loop 'while ~isempty(IndexForbidden)' take
%care of this.
                    if FCLocsReassign == 1
%CG: it is necessary to check the CNhood for the possible event, to see if
%any of these indices overlap with the elements of established CNhood
%elements. In the AARG_spine function, the element of each ROI of each
%event is screened against all the ROIs currently stored in the CA_ROI cell
%array.
                        if isempty(find(RepeaterEvents(:) == FakeEvent,1));

                            ReassignmentCtrl = 1;

                            while ReassignmentCtrl == 1
                                ReassignmentCtrl = 0;
                                PosAlternative = randsample(CurrentFCLocsRange, 1); 
                                [RowSubs, ColSubs, ~] = ind2sub(size(SLM), PosAlternative);
                                PosAlternative2D = sub2ind(size(SLM(:,:,1)), RowSubs, ColSubs);
                                if ~isempty(find(ForbiddenIndices(:) == PosAlternative2D, 1))
                                    ReassignmentCtrl = 1; 
                                elseif isempty(find(ForbiddenIndices(:) == PosAlternative2D, 1))
                                    FCLocs_sorted(FakeEvent) = PosAlternative;
                                end
                            end
                        end  
                    end

                    if ~isempty(find(RepeaterEvents(:) == FakeEvent,1))... 
                            && ~isempty(find(LargeEvents(:) == FakeEvent,1));
                        sesize = LargeEventSize;
                        if strcmp(LoadOrSaveIt, 'LoadIt')
                            RESizeCounter_sc = RESizeCounter_sc + 1;
                            FESizeCounter_sc = FESizeCounter_sc + 1;
                        end
                    elseif ~isempty(find(RepeaterEvents(:) == FakeEvent,1))...
                            && isempty(find(LargeEvents(:) == FakeEvent,1));
                        if strcmp(LoadOrSaveIt, 'LoadIt')                   
                            sesize = RESizeList_sc(RESizeCounter_sc);
                            RESizeCounter_sc = RESizeCounter_sc + 1;
                            FESizeCounter_sc = FESizeCounter_sc + 1;
                        else
                            sesize = randi(StdEventSize_Max);
                        end

                    elseif isempty(find(RepeaterEvents(:) == FakeEvent,1))...
                            && isempty(find(LargeEvents(:) == FakeEvent,1));
                        if strcmp(LoadOrSaveIt, 'LoadIt')
                            sesize = FESizeList_sc(FESizeCounter_sc);
                            FESizeCounter_sc = FESizeCounter_sc + 1;
                        else
                            sesize = randi(StdEventSize_Max);
                        end
                    end
                    
                    se = strel('disk', sesize, 0);

                    FEM = se.getnhood;
%CG: FEM = Fake Event Matrix 
                    [xdim_se, ydim_se] = size(FEM);

                    T = round(((xdim_se)/2));
                    S = round(((ydim_se)/2));

                    CInd_FEM = sub2ind(size(FEM), T, S);
%CG: CInd_FEM is the centroid location of the fake event matrix (FEM).
                    AllPosInds_FEM = find(FEM);

                    [RowVec_FEM, ColVec_FEM, ZVec_FEM] = ind2sub(size(FEM), AllPosInds_FEM);
                    [RowVec_CIndFEM, ColVec_CIndFEM, ZVec_CIndFEM] = ind2sub(size(FEM), CInd_FEM);

                    SubtVec_R = ones(size(RowVec_FEM)).*RowVec_CIndFEM;
                    SubtVec_C = ones(size(ColVec_FEM)).*ColVec_CIndFEM;
                    SubtVec_Z = ones(size(ZVec_FEM)).*ZVec_CIndFEM;

                    PosiVec_R = RowVec_FEM - SubtVec_R;
                    PosiVec_C = ColVec_FEM - SubtVec_C;
                    PosiVec_Z = ZVec_FEM - SubtVec_Z;

%CG: '0' in the positioning vector (PosiVec...) should line up with FCLocs(x),
%where x is any linear index value. Only the Centroid subscripts should
%have all 0 values. 

                    if ~isempty(find(RepeaterEvents(:) == FakeEvent,1));

                        if strcmp(LoadOrSaveIt, 'LoadIt')
                            
                            if ~isempty(find(LargeEvents(:) == FakeEvent,1)) 

                                loaded_FixedLoc = LEOrderArray_fl(LEOrderCounter_sc);
                                Ind_FixedLoc = LEOrderArray_sc(LEOrderCounter_sc);
                                FixedLocsElement = FixedLocs(Ind_FixedLoc);
                                LEOrderCounter_sc = LEOrderCounter_sc + 1;

                                Ind_FixedLoc4sc = Ind_FixedLoc;
                                FixedLocs_4LEs = FixedLocs;

                                REOrderCounter_sc = REOrderCounter_sc + 1;
                                FEOrderCounter_sc = FEOrderCounter_sc + 1;
                            elseif isempty(find(LargeEvents(:) == FakeEvent,1))  && ~isempty(find(RepeaterEvents(:) == FakeEvent,1))

                                loaded_FixedLoc = REOrderArray_fl(REOrderCounter_sc);
                                Ind_FixedLoc = REOrderArray_sc(REOrderCounter_sc);
                                FixedLocsElement = FixedLocs(Ind_FixedLoc);
%CG: zero-value elements are placed in the REOrderArray_sc to represent the
%LargeEvents. Non-Large repeater events must not be zero values. 

                                REOrderCounter_sc = REOrderCounter_sc + 1;
                                FEOrderCounter_sc = FEOrderCounter_sc + 1;
                                
%                             elseif isempty(find(RepeaterEvents(:) == FakeEvent,1))
%                                 Ind_FixedLoc = FEOrderArray_sc(FEOrderCounter_sc);
%                                 FEOrderCounter_sc = FEOrderCounter_sc + 1;
                            end
                        else
%CG: if the current event is a LargeEvent or a RepeaterEvent, then the
%modelling script will select a location for this event from the FixedLocs
%list. If this is found later to overlap with another event (using the
%SLM_TouchPrevention matrix or the SLM_ROI_TP matrix) then another location
%will be randomly assigned. 
                            if NRL_Counter == 1
                                FixedLocsElement = randsample(FixedLocs, 1);
                                Ind_FixedLoc = find(FixedLocs == FixedLocsElement);
                            elseif NRL_Counter == 0 && isempty(find(LargeEvents(:) == FakeEvent,1));   
                                FixedLocsElement = randsample(FixedLocs, 1);
                                Ind_FixedLoc = find(FixedLocs == FixedLocsElement);

                            elseif NRL_Counter == 0 && ~isempty(find(LargeEvents(:) == FakeEvent,1));
                                LargeEventCounter = LargeEventCounter + 1;
                                Sz_FixedLocs_4LEs = numel(FixedLocs_4LEs);     
                                if Sz_FixedLocs_4LEs > 1

                                    FixedLocsElement = randsample(FixedLocs_4LEs, 1);                 
                                    FixedLocs_Mem = FixedLocs_4LEs;
                                    Ind_2Remove = FixedLocs_4LEs == FixedLocsElement;
%CG: The original, unadulterated FixedLocs HAS TO be used here in order to
%find the unique and correct order of indices that must be applied to
%insert the LargeEvents into the SLM, so that LEs do not occur at more than
%one location. 
                                    Ind_FixedLoc = find(FixedLocs(:) == FixedLocsElement);
                                    Ind_FixedLoc4sc = FixedLocs_4LEs == FixedLocsElement;

                                elseif Sz_FixedLocs_4LEs == 1

                                    FixedLocsElement = FixedLocs_4LEs;
                                    FixedLocs_Mem = FixedLocs_4LEs;
                                    Ind_2Remove = FixedLocs_4LEs == FixedLocsElement;

                                    Ind_FixedLoc = find(FixedLocs(:) == FixedLocsElement);
                                    Ind_FixedLoc4sc = FixedLocs_4LEs == FixedLocsElement;

                                elseif Sz_FixedLocs_4LEs < 1

                                    disp('Error: No more FixedLocs elements for LEs!')
                                end
                            end
                        end

                        if ~isempty(find(LargeEvents(:) == FakeEvent,1));
                            if strcmp(LoadOrSaveIt, 'LoadIt')
                                Ind_CRepeat = loaded_FixedLoc; 
                                Ind_CRepeatTS1 = Ind_CRepeat - (xdim*ydim*(timestep-1));
                            elseif strcmp(LoadOrSaveIt, 'SaveIt')
                                Ind_CRepeat = FixedLocs_4LEs(Ind_FixedLoc4sc) + (xdim*ydim*(timestep-1));
                                Ind_CRepeatTS1 = FixedLocs_4LEs(Ind_FixedLoc4sc);
                            end
                        elseif isempty(find(LargeEvents(:) == FakeEvent,1));
                            if strcmp(LoadOrSaveIt, 'LoadIt')
                                Ind_CRepeat = loaded_FixedLoc;
                                Ind_CRepeatTS1 = Ind_CRepeat - (xdim*ydim*(timestep-1));
                            elseif strcmp(LoadOrSaveIt, 'SaveIt')
                                Ind_CRepeat = FixedLocs(Ind_FixedLoc) + (xdim*ydim*(timestep-1));
                                Ind_CRepeatTS1 = FixedLocs(Ind_FixedLoc); 
                            end
                        end
%CG: FixedLocs provides the spatial information about where the current
%event should be placed (approximately) to get a spatially overlapping
%event. To get the correct index value for a timestep that is >1 the
%expression "...+ (xdim*ydim*(timestep-1))" must be used.

                        [XCentroid, YCentroid, ~] = ind2sub(size(SLM), Ind_CRepeat);

%CG: Get the vector coordinates to find the centroid. We also want to get
%the Centroid Neighborhood (CNhood). The code will randomly choose an index
%from the CNhood and then make this the new centroid location for the
%current FakeEvent.
                        Row_CNhoodDiag = XCentroid-floor(ROI_sl/2):1:XCentroid+floor(ROI_sl/2);
                        Col_CNhoodDiag = YCentroid-floor(ROI_sl/2):1:YCentroid+floor(ROI_sl/2);

                        for cEl = 1 : ROI_sl
                            Row_CNhood(cEl:ROI_sl:ROI_sl^2) = Row_CNhoodDiag(cEl);
                        end

                        for cEl = 1 : ROI_sl
                            Col_CNhood(((cEl-1)*ROI_sl)+1:1:ROI_sl*cEl) = Col_CNhoodDiag(cEl);
                        end

                        Z_CTs = ones(numel(Row_CNhood),1);
                        Z_CTs = Z_CTs';
                        Z_CTs = Z_CTs.*timestep;
                        Ind_CentroidNhood = sub2ind(size(SLM), Row_CNhood, Col_CNhood, Z_CTs);
%CG: NewCentroidLoc is a randomly selected point in Ind_CentroidNhood. This is 
%equivalent to the centroid of a real event overlapping with an old event.
                        if isempty(find(LargeEvents(:) == FakeEvent,1));
                            
                            cEventType = 'RepeaterEvent';
                            if strcmp(LoadOrSaveIt, 'SaveIt')
                                NewCentroidLoc = randsample(Ind_CentroidNhood, 1);
                            elseif strcmp(LoadOrSaveIt, 'LoadIt')
                                NewCentroidLoc = Ind_CRepeat;
                            end
                            
                            [RowVec_FCLocs, ColVec_FCLocs, ZVec_FCLocs] = ind2sub(size(SLM(:,:,timestep)), NewCentroidLoc);

%CG: The NewCentroidLoc should be the variable used to create the vector arrays:
%RowVec_FCLocs, ColVec_FCLocs, ZVec_FCLocs for RepeaterEvents. However,
%the CNhood for the NewCentroidLoc must be found for the SLM_CA_ROIs cell
%array :). 
                            [XNewCentroid, YNewCentroid, ~] = ind2sub(size(SLM), NewCentroidLoc);

                            Row_NewCNhoodDiag = XNewCentroid-floor(ROI_sl/2):1:XNewCentroid+floor(ROI_sl/2);
                            Col_NewCNhoodDiag = YNewCentroid-floor(ROI_sl/2):1:YNewCentroid+floor(ROI_sl/2);

                            for cEl = 1 : ROI_sl
                                Row_NewCNhood(cEl:ROI_sl:ROI_sl^2) = Row_NewCNhoodDiag(cEl);
                            end

                            for cEl = 1 : ROI_sl
                                Col_NewCNhood(((cEl-1)*ROI_sl)+1:1:ROI_sl*cEl) = Col_NewCNhoodDiag(cEl);
                            end

                            Ind_NewCentroidNhood = sub2ind(size(SLM), Row_NewCNhood, Col_NewCNhood);

                        elseif ~isempty(find(LargeEvents(:) == FakeEvent,1));

                            cEventType = 'LargeEvent';
                            [RowVec_FCLocs, ColVec_FCLocs, ZVec_FCLocs] = ind2sub(size(SLM(:,:,timestep)), Ind_CRepeat);

                            Ind_CentroidNhood = sub2ind(size(SLM), Row_CNhood, Col_CNhood);
                            Ind_NewCentroidNhood = Ind_CentroidNhood; 
                        end
                    elseif isempty(find(RepeaterEvents(:) == FakeEvent,1));
%CG: If the current FakeEvent does not happen to be in the RepeaterEvent
%array, then it is a FreeEvent and can simply use the the centroid location
%randomly selected and stored in the FCLocs list, but with the corrective
%term for the current timestep: '...+ (xdim*ydim*(timestep-1))'.
                        cEventType = 'FreeEvent';
                        if strcmp(LoadOrSaveIt, 'LoadIt')
                            loaded_FreeLoc = FEOrderArray_sc(FEOrderCounter_sc);
                            FEOrderCount_rec = [FEOrderCount_rec, FEOrderCounter_sc];
                            [RowVec_FCLocs, ColVec_FCLocs, ZVec_FCLocs] = ind2sub(size(SLM), loaded_FreeLoc);
                            FEOrderCounter_sc = FEOrderCounter_sc + 1;
                        elseif strcmp(LoadOrSaveIt, 'SaveIt')
                            [RowVec_FCLocs, ColVec_FCLocs, ZVec_FCLocs] = ind2sub(size(SLM), FCLocs_sorted(FakeEvent));
                        end
                        Row_FreeEventCNhoodDiag = RowVec_FCLocs-floor(ROI_sl/2):1:RowVec_FCLocs+floor(ROI_sl/2);
                        Col_FreeEventCNhoodDiag = ColVec_FCLocs-floor(ROI_sl/2):1:ColVec_FCLocs+floor(ROI_sl/2);

                        for cEl = 1 : ROI_sl
                            Row_FreeEventCNhood(cEl:ROI_sl:ROI_sl^2) = Row_FreeEventCNhoodDiag(cEl);
                        end

                        for cEl = 1 : ROI_sl
                            Col_FreeEventCNhood(((cEl-1)*ROI_sl)+1:1:ROI_sl*cEl) = Col_FreeEventCNhoodDiag(cEl);
                        end
                        Ind_FreeEventCentroidNhood = sub2ind(size(SLM(:, :, 1)), Row_FreeEventCNhood, Col_FreeEventCNhood);

                    end 

%CG: The last condition statement provides the centroid for the current
%event, which is either a 'Repeater' or a 'Free' event. Prior to the last condition 
%statement the current Fake Event was built along with the positioning
%vectors. Now we have the position vectors and information about where the
%event should be. With this information, the code can position the current
%FakeEvent at the appropriate location in the SLM. 
                    TempVec_RowFCLocs = ones(size(RowVec_FEM)).*RowVec_FCLocs;
                    TempVec_ColFCLocs = ones(size(ColVec_FEM)).*ColVec_FCLocs;
                    TempVec_ZFCLocs = ones(size(ZVec_FEM)).*ZVec_FCLocs;

                    FinalRowVector = TempVec_RowFCLocs + PosiVec_R;
                    FinalColVector = TempVec_ColFCLocs + PosiVec_C;
                    FinalZVector = TempVec_ZFCLocs + PosiVec_Z;

                    Ind2Test_TP = sub2ind(size(SLM_TouchPrevention), FinalRowVector(:), FinalColVector(:)); 

                    SLM_TouchPrevention(Ind2Test_TP) = 1;

                    if isempty(find(RepeaterEvents(:) == FakeEvent,1))
                        SLM_ROI_TP(Ind_FreeEventCentroidNhood) = 1;
                    end

                    CC = bwconncomp(SLM_TouchPrevention);
                    NumComponents = CC.NumObjects;

                    CC_ROI = bwconncomp(SLM_ROI_TP);
                    NumComponents_ROI = CC_ROI.NumObjects;
                    if timestep > 1
                        CumulFECount = sum(FECounterRecord(1:end));
                    elseif timestep == 1
                        CumulFECount = 0;
                    end

                    if isempty(FixedLocs_TakenThisRound)
                        TwoIsACrowd = [];
                    elseif ~isempty(FixedLocs_TakenThisRound) && ~isempty(find(RepeaterEvents(:) == FakeEvent,1))
                        if isempty(find(LargeEvents(:) == FakeEvent,1))
                            TwoIsACrowd = find(FixedLocs_TakenThisRound == FixedLocs(Ind_FixedLoc));
                        elseif ~isempty(find(LargeEvents(:) == FakeEvent,1))
                            TwoIsACrowd = find(FixedLocs_TakenThisRound == FixedLocs(Ind_FixedLoc));
                        end
                    end

                   
                    if NumComponents == FakeEventCounter &&...
                            ~isempty(find(RepeaterEvents(:) == FakeEvent,1)) && isempty(TwoIsACrowd)

                        Ind2Convert = sub2ind(size(SLM), FinalRowVector(:), FinalColVector(:), FinalZVector(:));
                        Ind2Convert_Compiled = sub2ind(size(SLM_Compiled), FinalRowVector(:), FinalColVector(:));
%CG: Ind2Convert is an array containing the final index values for the SLM,
%that should be covered by the current FakeEvent. Ind2Covert_Compiled is
%intend for a 2D SLM, when the user wants to see the location of all
%FakeEvents, regardless of their timing.
                        SLM(Ind2Convert) = 1;   
                        SLM_Compiled(Ind2Convert_Compiled) = 1;   
                        SLM_TouchPrevention_ii(Ind2Test_TP) = 1;

                        if isempty(FixedLocs_TakenThisRound)
                            FixedLocs_TakenThisRound(1) = FixedLocs(Ind_FixedLoc);
                        elseif ~isempty(FixedLocs_TakenThisRound)
                            FixedLocs_TakenThisRound(end+1) = FixedLocs(Ind_FixedLoc);
                        end

                        if ~isempty(find(LargeEvents(:) == FakeEvent,1))
                            if ~strcmp(LoadOrSaveIt, 'LoadIt')
                                FixedLocs_4LEs(Ind_2Remove) = [];
%                                 FixedLocs_Mem = FixedLocs_4LEs;
                            end
                            if strcmp(LoadOrSaveIt, 'SaveIt')

                                LEOrderArray_sc(LEOrderCounter_sc) = Ind_FixedLoc;
                                LEOrderArray_fl(LEOrderCounter_sc) = Ind_CRepeat;
                                LEOrderCounter_sc = LEOrderCounter_sc + 1;
                                if isempty(REOrderCounter_scList);
                                    REOrderCounter_scList = REOrderCounter_sc;
                                else
                                    REOrderCounter_scList(end+1) = REOrderCounter_sc;
                                end
                                REOrderCounter_sc = REOrderCounter_sc + 1;
                                FEOrderCounter_sc = FEOrderCounter_sc + 1;
                                RESizeCounter_sc = RESizeCounter_sc + 1;
                                FESizeCounter_sc = FESizeCounter_sc + 1; 
%CG: If the current event is a LargeEvent, then the modelling script needs 
%to by-pass the current value of the counter for the SizeList and
%OrderArray of the FreeEvents and/or the RepeaterEvents.
                            end

                        elseif isempty(find(LargeEvents(:) == FakeEvent,1))

                            if strcmp(LoadOrSaveIt, 'SaveIt')

                                REOrderArray_fl(REOrderCounter_sc) = NewCentroidLoc;
                                REOrderArray_sc(REOrderCounter_sc) = Ind_FixedLoc;
                                REOrderCounter_sc = REOrderCounter_sc + 1;
                                FEOrderCounter_sc = FEOrderCounter_sc + 1;

                                RESizeList_sc(RESizeCounter_sc) = sesize;
                                RESizeCounter_sc = RESizeCounter_sc + 1; 

                                FESizeCounter_sc = FESizeCounter_sc + 1;
                            end
                        end

                    elseif NumComponents == FakeEventCounter && isempty(find(RepeaterEvents(:) == FakeEvent,1))...
                            && NumComponents_ROI == (FakeEventCounter + NumRepeaterLocs - RepeaterEventCounter + CumulFECount)

                        Ind2Convert = sub2ind(size(SLM), FinalRowVector(:), FinalColVector(:), FinalZVector(:));
                        Ind2Convert_Compiled = sub2ind(size(SLM_Compiled), FinalRowVector(:), FinalColVector(:));

                        SLM(Ind2Convert) = 1;   
                        SLM_Compiled(Ind2Convert_Compiled) = 1;    
                        SLM_TouchPrevention_ii(Ind2Test_TP) = 1;

                        SLM_ROI_TP_ii(Ind_FreeEventCentroidNhood) = 1;

                        if strcmp(LoadOrSaveIt, 'SaveIt')
                            FESizeList_sc(FESizeCounter_sc) = sesize;
                            FESizeCounter_sc = FESizeCounter_sc + 1;
                            
                            FEOrderArray_sc(FEOrderCounter_sc) = FCLocs_sorted(FakeEvent);
                            FEOrderCounter_sc = FEOrderCounter_sc + 1;
                        end

                    elseif NumComponents == FakeEventCounter && isempty(find(RepeaterEvents(:) == FakeEvent,1))...
                            && NumComponents_ROI < (FakeEventCounter + NumRepeaterLocs - RepeaterEventCounter + CumulFECount)

%CG: some small events can have a total size that is less than the size of
%their ROIs. In this case, it is necessary to check that the ROIs are not
%touching or overlapping. The number of ROIs visible in the SLM_ROI_TP
%matrix should be exactly (FakeEventCounter + NumRepeaterLocs -
%RepeaterEventCounter + CumulFECount)

                        if isempty(find(RepeaterEvents(:) == FakeEvent,1));
                            RepeaterC = RepeaterC + 1;
                        end
                        TouchingEvents = 1;
                        FCLocsReassign = 1;

                        Ind_ROI_TP = SLM_ROI_TP_ii == 1;
                        SLM_ROI_TP(:) = 0;
                        SLM_ROI_TP(Ind_ROI_TP) = 1;

                        Ind_TP = SLM_TouchPrevention_ii == 1;
                        SLM_TouchPrevention(:) = 0;
                        SLM_TouchPrevention(Ind_TP) = 1;

                    elseif NumComponents < FakeEventCounter || ~isempty(TwoIsACrowd)
%CG: This can only mean that the current Fake Event has touched or
%overlapped a preexisting event in the current timestep. Therefore, the
%current FCLocs index (FCLocs(FakeEventCounter)) should be reselected.

                        if isempty(find(RepeaterEvents(:) == FakeEvent,1));
                            RepeaterC = RepeaterC + 1;
                        end

                        if ~isempty(find(RepeaterEvents(:) == FakeEvent,1));
                            RepeaterK = RepeaterK + 1;
                            if RepeaterK == 1 
                                IdxCheck = zeros(1);
                                IdxCheck(1) = FixedLocsElement;
                                IdxCheck_u = IdxCheck;
                            elseif RepeaterK > 1
                                IdxCheck(end+1) = FixedLocsElement;
                                IdxCheck_u = unique(IdxCheck);
                            end
                            if numel(IdxCheck_u) > NumRepeaterLocs
                                errordlg('Error3: This function has attempted to find a FixedLoc for the current RepeaterEvent, where the current RepeaterEvent will not overlap spatially and temporally with any other event in the current timestep. In this case, it was not possible. The chance of this error recurring can be reduced by increasing the NumRepeaterLocs.');
                                Error3Flag = 1;
                                BreakVariable = NonExistentVariable;
                                break
                            end
                        end
                        TouchingEvents = 1;
                        FCLocsReassign = 1;
%CG: In this case SLM_TouchPrevention should revert back to its condition 
%before the current loop. For this, we can use the copy stored in
%SLM_TouchPrevention_ii.
                        Ind_ROI_TP = SLM_ROI_TP_ii == 1;
                        SLM_ROI_TP(:) = 0;
                        SLM_ROI_TP(Ind_ROI_TP) = 1;
                        Ind_TP = SLM_TouchPrevention_ii == 1;
                        SLM_TouchPrevention(:) = 0;
                        SLM_TouchPrevention(Ind_TP) = 1;

                       if ~isempty(find(LargeEvents(:) == FakeEvent,1));
%                            matVis(SLM(:,:,timestep))
                           FixedLocs_4LEs = FixedLocs_Mem;
                       end

                    elseif NumComponents > (FakeEventCounter + NumRepeaterLocs)

                        figure('Name', 'TouchPrevention Mat'), imshow(SLM_TouchPrevention);
                        figure('Name', 'TouchPrevention_ii Mat'), imshow(SLM_TouchPrevention_ii);
                        disp('error!')

                    end
                    
                end

                if ~isempty(find(RepeaterEvents(:) == FakeEvent,1)) && ~isempty(find(LargeEvents(:) == FakeEvent,1));

                    Sz_Ind2Convert_Compiled = size(Ind2Convert_Compiled);
                    if Sz_Ind2Convert_Compiled(1) > Sz_Ind2Convert_Compiled(2)
                        Ind2Convert_Compiled = Ind2Convert_Compiled';
                    end
                    
                    Str_FCLocs(FakeEventCounter).TypeStr = 'LargeEvent';
                    Str_FCLocs(FakeEventCounter).CIndexStr = Ind_CRepeat;
                    Str_FCLocs(FakeEventCounter).FakeEventCounterStr = FakeEventCounter;
                    Str_FCLocs(FakeEventCounter).EventSizeStr = Ind2Convert_Compiled;
                    Str_FCLocs(FakeEventCounter).CIndexTS1Str = Ind_CRepeatTS1;
                    Str_FCLocs(FakeEventCounter).NCNhoodStr = Ind_NewCentroidNhood;
                    Str_FCLocs(FakeEventCounter).FLElement = FixedLocsElement;
                    Str_FCLocs(FakeEventCounter).CellOverlapTypeStr = [];
%CG: if the FixedLocsElement equals an element value in the DistinguishedFixedLocs
%array, then it is necessary to designate this, because then this large
%event will not be treated like a normal LargeEvent and generate a ROI. 
                    if KeepROIs == 1
                        if ~isempty(find(DistinguishedFixedLocs == FixedLocsElement, 1))
                            Str_FCLocs(FakeEventCounter).DOverlapTypeStr = 'D_Overlap';
                        elseif isempty(find(DistinguishedFixedLocs == FixedLocsElement, 1))
                            Str_FCLocs(FakeEventCounter).DOverlapTypeStr = '';
                        end
                    elseif KeepROIs == 0
                        Str_FCLocs(FakeEventCounter).DOverlapTypeStr = '';
                    end
                    Res_Ind_NewCentroidNhood = Ind_NewCentroidNhood';
                    ForbiddenIndices = cat(1, ForbiddenIndices, Res_Ind_NewCentroidNhood); 

                elseif ~isempty(find(RepeaterEvents(:) == FakeEvent,1)) && isempty(find(LargeEvents(:) == FakeEvent,1));
                    
                    Sz_Ind2Convert_Compiled = size(Ind2Convert_Compiled);
                    if Sz_Ind2Convert_Compiled(1) > Sz_Ind2Convert_Compiled(2)
                        Ind2Convert_Compiled = Ind2Convert_Compiled';
                    end
                    
                    Str_FCLocs(FakeEventCounter).TypeStr = 'RepeaterEvent';
                    Str_FCLocs(FakeEventCounter).CIndexStr = NewCentroidLoc;
                    Str_FCLocs(FakeEventCounter).FakeEventCounterStr = FakeEventCounter;
                    Str_FCLocs(FakeEventCounter).EventSizeStr = Ind2Convert_Compiled;
                    Str_FCLocs(FakeEventCounter).CIndexTS1Str = Ind_CRepeatTS1;
                    Str_FCLocs(FakeEventCounter).NCNhoodStr = Ind_NewCentroidNhood;
                    Str_FCLocs(FakeEventCounter).FLElement = FixedLocsElement;
                    Str_FCLocs(FakeEventCounter).CellOverlapTypeStr = [];
                    if KeepROIs == 1
                        if ~isempty(find(DistinguishedFixedLocs == FixedLocsElement, 1))
                            Str_FCLocs(FakeEventCounter).DOverlapTypeStr = 'D_Overlap';
                        elseif isempty(find(DistinguishedFixedLocs == FixedLocsElement, 1))
                            Str_FCLocs(FakeEventCounter).DOverlapTypeStr = '';
                        end
                    elseif KeepROIs == 0
                        Str_FCLocs(FakeEventCounter).DOverlapTypeStr = '';
                    end
                    Res_Ind_NewCentroidNhood = Ind_NewCentroidNhood';
                    ForbiddenIndices = cat(1, ForbiddenIndices, Res_Ind_NewCentroidNhood);

                elseif isempty(find(RepeaterEvents(:) == FakeEvent,1)) 
                    
                    Sz_Ind2Convert_Compiled = size(Ind2Convert_Compiled);
                    if Sz_Ind2Convert_Compiled(1) > Sz_Ind2Convert_Compiled(2)
                        Ind2Convert_Compiled = Ind2Convert_Compiled';
                    end
                    Str_FCLocs(FakeEventCounter).TypeStr = 'FreeEvent';
                    Str_FCLocs(FakeEventCounter).CIndexStr = FCLocs_sorted(FakeEvent);
                    Str_FCLocs(FakeEventCounter).FakeEventCounterStr = FakeEventCounter;
                    Str_FCLocs(FakeEventCounter).EventSizeStr = Ind2Convert_Compiled;
                    Str_FCLocs(FakeEventCounter).FECNhoodStr = Ind_FreeEventCentroidNhood;
                    Str_FCLocs(FakeEventCounter).DOverlapTypeStr = '';

                    if timestep == 1 
                        CheckingFEOverlaps = CheckingFEOverlaps + 1;
                        if CheckingFEOverlaps == 1
                            Ind_FreeEventCentroidNhood_store = Ind_FreeEventCentroidNhood;
                        else 
                            if Ind_FreeEventCentroidNhood(1) < Ind_FreeEventCentroidNhood_store(1)
                                 Ind_FreeEventCentroidNhood_store = Ind_FreeEventCentroidNhood;
                            end
                        end
                    end
                    Res_Ind_FreeEventCentroidNhood = Ind_FreeEventCentroidNhood';
                    ForbiddenIndices = cat(1, ForbiddenIndices, Res_Ind_FreeEventCentroidNhood);
                end    
            end

        end
        
    end
    
    if FakeEventStart ~= FakeEventEnd && FakeEventEnd >= FakeEventStart && FakeEventStart <= EventNum
        if Error2Flag == 1   
            break 
        end

        StructLimit = numel(FakeEvents_ThisTS);
        SortingArray = zeros(StructLimit, 1);

        for LooperIdxi = 1 : StructLimit
            SortingArray(LooperIdxi) = Str_FCLocs(LooperIdxi).CIndexStr;
        end

%CG: The SortingArray will be used by sort to retrieve the indices, which will 
%tell the script which order to insert the data, contained in the struct
%array, into the cell arrays. 

%CG: Events have to be sorted because AARG will order events detected in
%the same frame according to their linear index values. 

        [~, SortingIdx] = sort(SortingArray, 'ascend');
 
        FakeEventCounter_sorted = 0;

        Sz_SLM_CA = size(SLM_CA_ROIs);

        for LooperIdxii = 1 : StructLimit
 
            SortedIdx = SortingIdx(LooperIdxii);
            FakeEventCounter_sorted = FakeEventCounter_sorted + 1;
            if timestep == 327 && FakeEventCounter_sorted == 4
                stophere = 1;
            end

            lgc_PreEst = [];
            if strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')
                lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;
            end
            
                
            if KeepROIs == 1
                Idx_PE = find(PE_timestepList == timestep);
                if ~isempty(find(Idx_PE))
                    if PE_componentList(Idx_PE) == FakeEventCounter_sorted
                        Str_FCLocs(SortedIdx).CellOverlapTypeStr = 'PE_Cell_OverlapEvent';
                    else
                        Str_FCLocs(SortedIdx).CellOverlapTypeStr = '';
                    end
                else
                    Str_FCLocs(SortedIdx).CellOverlapTypeStr = '';
                end
            elseif KeepROIs == 0
                Str_FCLocs(SortedIdx).CellOverlapTypeStr = '';
            end
            
            if strcmp(Str_FCLocs(SortedIdx).TypeStr, 'LargeEvent') == 1;
%CG: if the current FakeEvent is a large event (LE), then we must record
%its index for the cell arrays. 

                counterLE = counterLE + 1; 

                DefLELocs(counterLE,:) = Str_FCLocs(SortedIdx).NCNhoodStr;
                MasterFEs(counterLE) = FakeEventCounter_sorted;
                timestepMFE(counterLE) = timestep;

                if numel(Sz_SLM_CA) == 2 && timestep <= MaxNumPEtimesteps && FakeEventCounter_sorted <= MaxNumPEComps
                    if isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted})
                        ZDim = 1;   
                    elseif ~isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted})
                        ZDim = 2;
                    end
                elseif numel(Sz_SLM_CA) == 3 && timestep <= MaxNumPEtimesteps && FakeEventCounter_sorted <= MaxNumPEComps
                    ZDim = 1;    
                    if isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted, 1}) && ZDim <= Sz_SLM_CA(3)
%CG: if SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted, 1} is
%empty, it can be that the more than one cell along the 3rd dimension is
%non-empty. The for-loop is used to find the non-empty cell with the
%current (timestep, FakeEventCounter_sorted) coordinates and having the
%largest z-coordinate value. 'ZDim = ZDim + 1' finds the next non-empty
%cell along the z-axis.
                        ZDim = [];
                        for cZDim = 1 : Sz_SLM_CA(3)
                            if ~isempty(SLM_CA_EventSize{timestep, FakeEventCounter_sorted, cZDim})
                                ZDim = cZDim;
                            end
                        end
                        if isempty(ZDim) 
                            ZDim = 1;
                        elseif ~isempty(ZDim)
                            ZDim = ZDim + 1;
                        end
                    elseif ~isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted,1}) && ZDim <= Sz_SLM_CA(3)
%CG: if SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted,1} is
%non-empty, this can mean that SLM_CA_EventSize{timestep,
%FakeEventCounter_sorted,1+n} is also non-empty, where n is any positive
%integer. In this case, the possibility that there is an empty cell between
%two non-empty cells is accommodated. This should not occur with the actual
%data, but there is no logical reason why the model should have to mimic
%this in order to remain a valid model.

%CG: Note that SLM_CA_EventSize rather than SLM_CA_ROIs must be used to
%check cell availablility because a non-empty SLM_CA_EventSize cell will be
%empty in SLM_CA_ROIs if the event was a RepeaterEvent.

                        for cZDim = 1 : Sz_SLM_CA(3)
                            if ~isempty(SLM_CA_EventSize{timestep, FakeEventCounter_sorted, cZDim})
                                ZDim = cZDim;
                            end
                        end
                        if isempty(ZDim) 
                            ZDim = 2;
                        elseif ~isempty(ZDim)
                            ZDim = ZDim + 1;
                        end
                    end 
                elseif timestep > MaxNumPEtimesteps || FakeEventCounter_sorted > MaxNumPEComps
                    ZDim = 1;
                end

                zMFE(counterLE) = ZDim;
%CG: Since the SLM_CAs can expand in the 3rd dimension if a pre-established
%ROI is already occupying the current cell, the next available ROI needs to
%be found. In theory, providing there are only two timelapse sequences per
%experiment, the 3rd dim of the SLM_CAs should not expand beyond 2. 
            elseif strcmp(Str_FCLocs(SortedIdx).TypeStr, 'RepeaterEvent') == 1;

                Ind_PosLELocs = find(PosLELocs(:) == Str_FCLocs(SortedIdx).CIndexTS1Str);

                Mod_Ind_PosLELocs = Ind_PosLELocs(:);
%CG: the current event has been added to the PosLELocs list as the last
%index. For the script to work correctly a list without this index needs to
%be used. 
                counterNLE = counterNLE + 1;
%CG: if the current FakeEvent is a standard event, then we must tally it 
%by increasing the value of the approriate element in MasterTally.
                MasterTally(counterNLE) = +1;
                PosLELocs(counterNLE,:) = Str_FCLocs(SortedIdx).NCNhoodStr;
                RepeaterFEs(counterNLE) = FakeEventCounter_sorted;
                timestepRFE(counterNLE) = timestep;

                if numel(Sz_SLM_CA) == 2 && timestep <= MaxNumPEtimesteps && FakeEventCounter_sorted <= MaxNumPEComps
                    if isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted})
                        ZDim = 1;
                    elseif ~isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted})
                        ZDim = 2;
                    end
                elseif numel(Sz_SLM_CA) == 3 && timestep <= MaxNumPEtimesteps && FakeEventCounter_sorted <= MaxNumPEComps
                    
                    if isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted, 1}) 
    
                        ZDim = [];
                        for cZDim = 1 : Sz_SLM_CA(3)
                            if ~isempty(SLM_CA_EventSize{timestep, FakeEventCounter_sorted, cZDim})
                                ZDim = cZDim;
                            end
                        end
                        if isempty(ZDim) 
                            ZDim = 1;
                        elseif ~isempty(ZDim)
                            ZDim = ZDim + 1;
                        end
                    elseif ~isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted}) && ZDim <= Sz_SLM_CA(3)
                         
                        for cZDim = 1 : Sz_SLM_CA(3)
                           if ~isempty(SLM_CA_EventSize{timestep, FakeEventCounter_sorted, cZDim})
                               ZDim = cZDim;
                           end
                        end
                        if isempty(ZDim) 
                            ZDim = 2;
                        elseif ~isempty(ZDim)
                            ZDim = ZDim + 1;
                        end
                    end 

                elseif timestep > MaxNumPEtimesteps || FakeEventCounter_sorted > MaxNumPEComps
                    ZDim = 1;
                end
                zRFE(counterNLE) = ZDim;

            elseif strcmp(Str_FCLocs(SortedIdx).TypeStr, 'FreeEvent') == 1;

                if numel(Sz_SLM_CA) == 2 && timestep <= MaxNumPEtimesteps && FakeEventCounter_sorted <= MaxNumPEComps
                    if isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted})
                        ZDim = 1;
                    elseif ~isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted})
                        ZDim = 2;
                    end
                elseif numel(Sz_SLM_CA) == 3 && timestep <= MaxNumPEtimesteps && FakeEventCounter_sorted <= MaxNumPEComps
                    ZDim = 1;
                    if isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted, 1}) && ZDim <= Sz_SLM_CA(3)
                        ZDim = [];
                        for cZDim = 1 : Sz_SLM_CA(3)
                            if ~isempty(SLM_CA_EventSize{timestep, FakeEventCounter_sorted, cZDim})
                                ZDim = cZDim;
                            end
                        end
                        if isempty(ZDim) 
                            ZDim = 1;
                        elseif ~isempty(ZDim)
                            ZDim = ZDim + 1;
                        end
                        
                    elseif ~isempty(SLM_CA_EventSize_PreEst{timestep, FakeEventCounter_sorted, 1}) && ZDim <= Sz_SLM_CA(3)
                        for cZDim = 1 : Sz_SLM_CA(3)
                           if ~isempty(SLM_CA_EventSize{timestep, FakeEventCounter_sorted, cZDim})
                               ZDim = cZDim;
                           end
                        end
                        if isempty(ZDim) 
                            ZDim = 2;
                        elseif ~isempty(ZDim)
                            ZDim = ZDim + 1;
                        end
                    end 
                elseif timestep > MaxNumPEtimesteps || FakeEventCounter_sorted > MaxNumPEComps
                    ZDim = 1;
                end

            end 

            if strcmp(Str_FCLocs(SortedIdx).TypeStr, 'LargeEvent') == 1;
%CG: if the current model event is a Large Event (LE), it must have
%occurred at a repeater location. There can be 2 possibilities: 1) This is
%not the first model event at this location - standard-sized repeater
%events have arrived before it. 2) This is the first model event to have
%occurred at this location, in which case, Ind_PosLELocs will be empty.

                Ind_PosLELocs = find(PosLELocs(:) == Str_FCLocs(SortedIdx).CIndexTS1Str);
%CG: the first standard-sized model event at the current location will have 
%accumulated event numbers for stanard-sized events that appeared at the
%current location before the LE occurred. These events now belong to the
%LE. The first standard-sized model event at the current location should
%also no longer be represented by an ROI.
                Mod_Ind_PosLELocs = Ind_PosLELocs(:);
                if ~isempty(find(Mod_Ind_PosLELocs,1))
%CG: Possibility One: The index value of the current LEs centroid should
%match one, and only one, in the list of indices for all standard-sized
%model events created so far (contained in the array 'PosLELocs').
%'Ind_PosLELocs' represents the index value to find the correct centroid
%index value in the list 'PosLELocs'. 'Ind_PosLELocs' is then used to
%collect the number recording the total number of events previously having
%occurred at this location from the MasterTally array plus 1 (+1 is for the
%current LE). This value should then be entered into the current cell for
%the cell array SLM_CA_EventNum. 

%CG: If a standard-sized event has occurred at this location already, then the 
%script must find the index information which will point to the largest of
%these standard-sized events. Since there cannot be any other LEs at this
%location, this largest standard-sized event must have ROI status and a
%collection of events stored in its cell in the SLM_CA_EventNum cell array,
%which must now be allocated to the current LE. 
                    [rrLocs,~] = ind2sub(size(PosLELocs),Mod_Ind_PosLELocs);
                    
                    FEsToCheck = RepeaterFEs(rrLocs);
                    TSsToCheck = timestepRFE(rrLocs);
%                     [TSsToCheck,ts_idx] = timestepRFE(sort(rrLocs,'ascend')); 
                    ZzzToCheck = zRFE(rrLocs);
                        
                    NumLocsToCheck = numel(Mod_Ind_PosLELocs);
                    EventSizeList = zeros(numel(Mod_Ind_PosLELocs),1);

                    for Event = 1 : NumLocsToCheck
                        EventSizeList(Event) = numel(SLM_CA_EventSize{TSsToCheck(Event), FEsToCheck(Event), ZzzToCheck(Event)});
                    end 
                    LSE_idx = EventSizeList == max(EventSizeList); lim_TSsToCheck = TSsToCheck(LSE_idx);
                    [~, min_idx] = min(lim_TSsToCheck); SpecialStdEvent = find(TSsToCheck == lim_TSsToCheck(min_idx));
%                     SpecialStdEvent = find(EventSizeList(:) == LargestStdEvent);

                    if strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                        && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')
%CG: If the current model event overlaps with a pre-established ROI in the
%cell array ('PE_Cell_OverlapEvent'), then the next available cell along
%the z-axis needs to be used to hold the data for this event. In this case
%the current LargeEvent is also overlapping with a pre-established ROI in
%the actual matrix, so 'PE_Cell_OverlapEvent' status only matters for the
%purposes of filling in the EventSize data. All the other data must be
%entered into the correct cell for the pre-established ROI.
                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                       
                        [rrLocs,~] = ind2sub(size(PosLELocs),Mod_Ind_PosLELocs);
                        GrandMTT = sum(MasterTally(rrLocs)) + 1;
                        SLM_CA_EventNum(DME_row, DME_col, DME_z) = {[]};
                        SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT;

                        AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI(DME_row, DME_col, DME_z) = {[]};
                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

%CG: If the current model event overlaps with a pre-established ROI in the cell
%array ('PE_Cell_OverlapEvent'), then the next available cell along the z-axis needs to be used to
%hold the data for this event.

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_ROIs{timestep, FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                        [rrLocs,~] = ind2sub(size(PosLELocs),Mod_Ind_PosLELocs);
                        GrandMTT = sum(MasterTally(rrLocs)) + 1;
                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = GrandMTT;
%CG: It can happen that there is more than one event with the largest size. 
%In this case, the first event in the list gets priority. Therefore,
%SpecialStdEvent(x) where x = 1. 
                        AllPreFrames = SLM_CA_FrameIEI{TSsToCheck(SpecialStdEvent(end)), FEsToCheck(SpecialStdEvent(end)), ZzzToCheck(SpecialStdEvent(end))};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = AllFrames;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap') 
%CG: It is expected that there will be no practical difference between this condition and 
%the first condition. ZDim is already set, with the SLM_CA_EventSize array
%being built along the z-axis if necessary. 
                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        [rrLocs,~] = ind2sub(size(PosLELocs),Mod_Ind_PosLELocs);
                        GrandMTT = sum(MasterTally(rrLocs)) + 1;
                        SLM_CA_EventNum(DME_row, DME_col, DME_z) = {[]};
                        SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT;

                        AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI(DME_row, DME_col, DME_z) = {[]};
                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '') 

                        SLM_CA_ROIs{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                        [rrLocs,~] = ind2sub(size(PosLELocs),Mod_Ind_PosLELocs);
                        GrandMTT = sum(MasterTally(rrLocs)) + 1;
                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = GrandMTT;

                        AllPreFrames = SLM_CA_FrameIEI{TSsToCheck(SpecialStdEvent(end)), FEsToCheck(SpecialStdEvent(end)), ZzzToCheck(SpecialStdEvent(end))};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = AllFrames;

                        SLM_CA_FrameIEI{TSsToCheck(SpecialStdEvent(end)), FEsToCheck(SpecialStdEvent(end)), ZzzToCheck(SpecialStdEvent(end))} = [];
                        SLM_CA_EventNum{TSsToCheck(SpecialStdEvent(end)), FEsToCheck(SpecialStdEvent(end)), ZzzToCheck(SpecialStdEvent(end))} = [];
                        SLM_CA_ROIs{TSsToCheck(SpecialStdEvent(end)), FEsToCheck(SpecialStdEvent(end)), ZzzToCheck(SpecialStdEvent(end))} = [];

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;
                        SLM_CA_ROIs(TSsToCheck(SpecialStdEvent(end)), FEsToCheck(SpecialStdEvent(end)), ZzzToCheck(SpecialStdEvent(end))) = {[]};

                    end

                elseif isempty(find(Mod_Ind_PosLELocs,1))
%CG: Possibilty Two: This is the first Model Event to occur at this
%location. In this case, we simply assign the current cell in
%SLM_CA_EventNum its first value: 1. Also, the current timestep is entered
%into the appropriate cell of the CA_FrameIEI array.
                    if strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                        && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')
%CG: If the current model event overlaps with a pre-established ROI in the cell
%array ('PE_Cell_OverlapEvent'), then the next available cell along the z-axis needs to be used to
%hold the data for this event. In this case the current LargeEvent is also
%overlapping with a pre-established ROI in the actual matrix, so
%'PE_Cell_OverlapEvent' status only matters for the purposes of filling in
%the EventSize data. All the other data must be entered into the correct
%cell for the pre-established ROI. 
                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                          
                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;

                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = timestep;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');
                       
                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

%CG: If the current model event overlaps with a pre-established ROI in the cell
%array ('PE_Cell_OverlapEvent'), then the next available cell along the z-axis needs to be used to
%hold the data for this event.

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_ROIs{timestep, FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = 1;
%CG: It can happen that there is more than one event with the largest size. 
%In this case, the first event in the list gets priority. Therefore,
%SpecialStdEvent(x) where x = 1. 

                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = timestep;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap') 
%CG: It is expected that there will be no practical difference between this condition and 
%the first condition. ZDim is already set, with the SLM_CA_EventSize array
%being built along the z-axis if necessary. 
                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;

                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = timestep;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '') 

                        SLM_CA_ROIs(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_ROIs{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = 1;

                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = timestep;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                    end
                end

            elseif strcmp(Str_FCLocs(SortedIdx).TypeStr, 'RepeaterEvent') == 1;
%CG: if the current Fake Event is of a standard size and is a repeater event, then
%there are 3 further possibilities: 1) a large event (LE) has already occurred
%at this location 2) No LE has occurred at this location, but a
%standard-sized event has 3) The current event, which is standard-sized, is
%the first event to have occurred at this location.
                Ind_DefLELocs = find(DefLELocs(:) == Str_FCLocs(SortedIdx).CIndexTS1Str);

                if ~isempty(Ind_DefLELocs)
%CG: Possibility 1: If a LE has occurred at the current Model Event's location 
%we have to identify the cell holding the information about the this LE in
%the cell arrays. We extract the total number of events that have been
%already attributed to this LE (this value is represented by
%the variable 'CumulNumEvents')and add 1 (the current Fake Event) to this
%value. The value produced by 'CumulNumEvents + 1' is assigned to the cell
%containing the size information about the LE.
%The current Model Event is represented in SLM_CA_ROIs by the large event,
%so the current cell in SLM_CA_ROIs is left blank - as is the current cell
%in SLM_CA_EventNum. The current timestep is added to the list of frames in
%in the SLM_CA_FrameIEI with the cell coordinates: (MasterTS, MasterFake).

%CG: Possibility 1a: This RepeaterEvent has landed on a
%DistinguishedFixedLoc. In this case, the +1 frequency is not given to the
%established ROI representing the ideal ROI for the LargeEvent in a
%previous model frame, it is give to the pre-established ROI. The same goes
%for FrameIEI and EventType.

%CG: Possibility 1b: This RepeaterEvent's Size data would go to a cell that
%is already occupied. So the EventSize cell array (and possible the others)
%now have to expand along the z-axis. 
                    if strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                         && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')


                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                        if isempty(GrandMTT)
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                        else
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                        end

                        AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        [rrLocs,~] = ind2sub(size(DefLELocs),Ind_DefLELocs);
                        MasterFake = MasterFEs(rrLocs);
                        MasterTS = timestepMFE(rrLocs);
                        MasterZ = zMFE(rrLocs);

                        CumulNumEvents4Mtr = SLM_CA_EventNum{MasterTS, MasterFake, MasterZ} + 1; 
                        SLM_CA_EventNum{MasterTS, MasterFake, MasterZ} = CumulNumEvents4Mtr;
                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                        AllPreFrames = SLM_CA_FrameIEI{MasterTS, MasterFake, MasterZ};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI{MasterTS, MasterFake, MasterZ} = AllFrames;
                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]}; 

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')

                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                        if isempty(GrandMTT)
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                        else
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                        end

                        AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '') 

                        SLM_CA_ROIs(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                        [rrLocs,~] = ind2sub(size(DefLELocs),Ind_DefLELocs);
                        MasterFake = MasterFEs(rrLocs);
                        MasterTS = timestepMFE(rrLocs);
                        MasterZ = zMFE(rrLocs);

                        CumulNumEvents4Mtr = SLM_CA_EventNum{MasterTS, MasterFake, MasterZ} + 1; 
                        SLM_CA_EventNum{MasterTS, MasterFake, MasterZ} = CumulNumEvents4Mtr;
                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                        AllpreFrames = SLM_CA_FrameIEI{MasterTS, MasterFake, MasterZ};
                        AllFrames = cat(2, AllpreFrames, timestep);
                        SLM_CA_FrameIEI{MasterTS, MasterFake, MasterZ} = AllFrames;
                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;
                    end

                elseif isempty(Ind_DefLELocs) && ~isempty(Mod_Ind_PosLELocs)

%CG: Possibility 2: A standard-sized event has already claimed the location 
%of the current event and created a ROI. We need to compare the size of the 
%current event to the size of the standard event which is representing this 
%Repeater location. If the current event is larger, then it will attain ROI
%status and the previous event will lose it. The current event will also
%gain all of the events accumulated by the previous event (in the cell array 
%CA_EventNum. If the current event is smaller or the same size, it will
%count towards the previous event's tally in CA_EventNum. 

%CG: Possibility 2a ....
%CG: Possibility 2b ....
                    Mod_Ind_PosLELocs = Ind_PosLELocs(:);
                    [rrLocs,~] = ind2sub(size(PosLELocs),Mod_Ind_PosLELocs);
                    FEsToCheck = RepeaterFEs(rrLocs);
                    TSsToCheck = timestepRFE(rrLocs); 
                    ZzzToCheck = zRFE(rrLocs);

                    NumLocs = numel(Mod_Ind_PosLELocs);
                    EventSizeList = zeros(numel(Mod_Ind_PosLELocs),1);

                    for Sibling = 1 : NumLocs
                       EventSizeList(Sibling) = numel(SLM_CA_EventSize{TSsToCheck(Sibling), FEsToCheck(Sibling), ZzzToCheck(Sibling)});                   
                    end

                    LargestStdEvent = max(EventSizeList);
                    SizeCurrentEvent = numel(Str_FCLocs(SortedIdx).EventSizeStr);
                    
                    SSE_idx = EventSizeList == max(EventSizeList); lim_TSsToCheck = TSsToCheck(SSE_idx);
                    [~, min_idx] = min(lim_TSsToCheck); SpecialSibling = find(TSsToCheck == lim_TSsToCheck(min_idx));

%                     SpecialSibling = find(EventSizeList(:) == LargestStdEvent);
%CG: SpecialSibling represents the index giving the location in the cell array 
%of the largest standard sized events. Whether this previous event gets the
%number of the current event (i.e. +1 to its cell containing information
%about the number of events at this location) and if it retains its ROI
%status, depends if it is bigger or smaller than the current event.
%Rather than using the extensive swapping of cell contents, it might be
%possible to use the MasterTally array as was done with the condition for
%the occurence of an LE. 
                    if LargestStdEvent >= SizeCurrentEvent

                        if strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')

                            lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                            [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 
                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                            GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                            if isempty(GrandMTT)
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                            else
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                            end

                            AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                            AllFrames = cat(2, AllPreFrames, timestep);
                            SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                        elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                            CumulNumEvents4Mtr = SLM_CA_EventNum{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))}; 
                            SLM_CA_EventNum{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))} = CumulNumEvents4Mtr + 1;
                            SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                            AllPreFrames = SLM_CA_FrameIEI{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))};
                            AllFrames = cat(2, AllPreFrames, timestep);
                            SLM_CA_FrameIEI{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))} = AllFrames;
                            SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]}; 

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                        elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')

                            lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                            [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                            GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                            if isempty(GrandMTT)
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                            else
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                            end

                            AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                            AllFrames = cat(2, AllPreFrames, timestep);
                            SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                        elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '') 

                            SLM_CA_ROIs(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                            CumulNumEvents4Rtr = SLM_CA_EventNum{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))};
                            SLM_CA_EventNum{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))} = CumulNumEvents4Rtr + 1;
                            SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};

                            AllpreFrames = SLM_CA_FrameIEI{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))};
                            AllFrames = cat(2, AllpreFrames, timestep);
                            SLM_CA_FrameIEI{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))} = AllFrames;
                            SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]}; 

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;
                        end

                    elseif LargestStdEvent < SizeCurrentEvent

                        if strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')

                            lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                            [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                            GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                            if isempty(GrandMTT)
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                            else
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                            end

                            AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                            AllFrames = cat(2, AllPreFrames, timestep);
                            SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                            SLM_CA_EventNum(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_FrameIEI(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_ROIs{timestep, FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;
                            SLM_CA_ROIs(TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))) = {[]};

                            CumulNumEvents4Mtr = SLM_CA_EventNum{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))}; 
                            SLM_CA_EventNum(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventNum{timestep, FakeEventCounter_sorted, ZDim} = CumulNumEvents4Mtr + 1;
                            SLM_CA_EventNum(TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))) = {[]};

                            AllPreFrames = SLM_CA_FrameIEI{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))};
                            AllFrames = cat(2, AllPreFrames, timestep);
                            SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = AllFrames;
                            SLM_CA_FrameIEI(TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))) = {[]}; 

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

%                             SLM_CA_ROIs(timestep, FakeEventCounter_sorted, 1) = {[]};
%                             SLM_CA_EventSize(timestep, FakeEventCounter_sorted, 1) = {[]};

                        elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')

                            lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                            [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                            GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                            if isempty(GrandMTT)
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                            else
                                SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                            end

                            AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                            AllFrames = cat(2, AllPreFrames, timestep);
                            SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                        elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                            && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '') 

                            SLM_CA_ROIs(TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))) = {[]};
                            SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_ROIs{timestep, FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                            CumulNumEvents4Rtr = SLM_CA_EventNum{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))};
                            SLM_CA_EventNum(TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))) = {[]};
                            SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = CumulNumEvents4Rtr + 1;

                            AllpreFrames = SLM_CA_FrameIEI{TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))};
                            AllFrames = cat(2, AllpreFrames, timestep);
                            SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = AllFrames;
                            SLM_CA_FrameIEI(TSsToCheck(SpecialSibling(1)), FEsToCheck(SpecialSibling(1)), ZzzToCheck(SpecialSibling(1))) = {[]};

                            SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                            SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                            SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                        end                        
                    end
%CG: In this condition, the first index of Ind_PosLELocs will give the
%index location in the SLM (for a one dimensional SLM). The condition
%statement made earlier in the current function generates the RepeaterFEs
%array and timestepRFE array. RepeaterFake and RepeaterTS store the
%FakeEventCounter and timestep values, respectively, that are congruent
%with the dimensions of the SLM cell arrays. Ind_PosLELocs(1) will always
%point to the first standard-sized Fake Event that occurred in the current
%location. This means that the scalar values of RepeaterFake and RepeaterTS
%should always take the code to the location on the cell array
%SLM_CA_EventNum that represents the first standard-sized event at this
%location.
                elseif isempty(Ind_DefLELocs) && isempty(Mod_Ind_PosLELocs)
%CG: Possibility 3: the simplest case...when the current Fake Event is the
%first event at this location, it establishes the ROI and starts counting
%events for this location. The same thing happens in this condition as if
%the current Fake Event was an LE, but we also assign a ROI for the current
%event (because it is only appropriate to assign an ROI to a standard-sized
%event if it is the first event at this location.

%CG: Possibility 3a ....
%CG: Possibility 3b ....
%                     ZDim = 1;
                    if strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                         && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap')

                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                       
                        GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                        if isempty(GrandMTT)
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                        else
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                        end

                        AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;
                            
                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, 'PE_Cell_OverlapEvent')...
                         && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_ROIs{timestep, FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                        SLM_CA_EventNum(timestep, FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventNum{timestep, FakeEventCounter_sorted, ZDim} = 1;

                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = timestep; 
                            
                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                         && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, 'D_Overlap') 

                        lgc_PreEst = PreEst_cDFL == Str_FCLocs(SortedIdx).FLElement;

                        [DME_row, DME_col, DME_z] = ind2sub(size(SLM_CA_ROIs_PreEst), PreEst_cCell(lgc_PreEst));
%CG: If the current event overlaps with a pre-established ROI in the matrix
%itself ('D_Overlap'), then this LargeEvent does not establish a ROI. 

                        checkmROI = SLM_CA_ROIs{DME_row, DME_col, DME_z};

                        SLM_CA_ROIs(timestep, FakeEventCounter_sorted, ZDim) = {[]};

                        GrandMTT = SLM_CA_EventNum{DME_row, DME_col, DME_z};
                        if isempty(GrandMTT)
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = 1;
                        else
                            SLM_CA_EventNum{DME_row, DME_col, DME_z} = GrandMTT + 1;
                        end

                        AllPreFrames = SLM_CA_FrameIEI{DME_row, DME_col, DME_z};
                        AllFrames = cat(2, AllPreFrames, timestep);
                        SLM_CA_FrameIEI{DME_row, DME_col, DME_z} = AllFrames;
                        
                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = strcat(Str_FCLocs(SortedIdx).TypeStr, 'D');

                    elseif strcmp(Str_FCLocs(SortedIdx).CellOverlapTypeStr, '')...
                         && strcmp(Str_FCLocs(SortedIdx).DOverlapTypeStr, '')

                        SLM_CA_ROIs(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_ROIs{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).NCNhoodStr;

                        SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = 1;

                        SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = timestep;

                        SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;

                        SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                        SLM_CA_EventType{timestep,  FakeEventCounter_sorted} = Str_FCLocs(SortedIdx).TypeStr;
                    end
                end

            elseif strcmp(Str_FCLocs(SortedIdx).TypeStr, 'FreeEvent') == 1;

                SLM_CA_ROIs(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                SLM_CA_ROIs{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).FECNhoodStr;

                SLM_CA_EventNum(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                SLM_CA_EventNum{timestep,  FakeEventCounter_sorted, ZDim} = 1;

                SLM_CA_FrameIEI(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                SLM_CA_FrameIEI{timestep,  FakeEventCounter_sorted, ZDim} = timestep;
                
                SLM_CA_EventSize(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                SLM_CA_EventSize{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).EventSizeStr;
                
                SLM_CA_EventType(timestep,  FakeEventCounter_sorted, ZDim) = {[]};
                SLM_CA_EventType{timestep,  FakeEventCounter_sorted, ZDim} = Str_FCLocs(SortedIdx).TypeStr;

            end

        end
        
        SLM_TouchPrevention(:) = 0;
    end
    
    if ~mod(timestep, wbStep) || timestep == timesteps
        waitbar(timestep/timesteps, wb)
    end
    
end
close(wb)
numel(find(REOrderArray_sc == 0))

if ~strcmp(SaveString, 'NoSave')
    save(strcat('SLM_CellArrays_', DefString, '.mat'), 'SLM_CA_ROIs', 'SLM_CA_EventNum', 'SLM_CA_EventSize',...
        'SLM_CA_EventType', 'SLM_CA_FrameIEI')
end

if turnONimages == 1
    figure('Name', strcat('SLM_', DefString)), imshow3D(SLM); 
end

IdxSz_SLM_CA_ROIs = size(SLM_CA_ROIs);
Sz_SLM_CA_ROIs = IdxSz_SLM_CA_ROIs(1)*IdxSz_SLM_CA_ROIs(2);
TotalNumCells = Sz_SLM_CA_ROIs;
CellCount_ROIs = 0;
CellCount_EventNum = 0;
CellCount_EventSize = 0;

for CurrentCell = 1 : TotalNumCells

    CheckingCell_ROIs = SLM_CA_ROIs{CurrentCell};
    CheckingCell_EventNum = SLM_CA_EventNum{CurrentCell};
    CheckingCell_EventSize = SLM_CA_EventSize{CurrentCell};

    if ~isempty(CheckingCell_ROIs)
        CellCount_ROIs = CellCount_ROIs +1;
    end

    if ~isempty(CheckingCell_EventNum)
        CellCount_EventNum = CellCount_EventNum +1;
    end

    if ~isempty(CheckingCell_EventSize)
        CellCount_EventSize = CellCount_EventSize +1;
    end

end  

if strcmp(LoadOrSaveIt, 'SaveIt')

    save(strcat('FCLocs_sorted_', DefString,'.mat'), 'FCLocs_sorted');
    save(strcat('FEOrderArray_', DefString,'.mat'), 'FEOrderArray_sc','FEOrderArray_fl');
    save(strcat('REOrderArray_', DefString,'.mat'), 'REOrderArray_sc','REOrderArray_fl');
    save(strcat('LEOrderArray_', DefString,'.mat'), 'LEOrderArray_sc','LEOrderArray_fl');
    save(strcat('RESizeList_', DefString,'.mat'), 'RESizeList_sc')
    save(strcat('FESizeList_', DefString,'.mat'), 'FESizeList_sc')
    save(strcat('SLM_', DefString,'.mat'), 'SLM');
end

LEOrderArray_sc_B = LEOrderArray_sc;
REOrderArray_sc_B = REOrderArray_sc;
FESizeList_sc_B = FESizeList_sc;

save('OrderArraysBefore', 'LEOrderArray_sc_B', 'REOrderArray_sc_B', 'FESizeList_sc_B') 

Ind_ROIArrayII = [];
Ind_EventArrayII = [];
for CurrentCell = 1 : TotalNumCells

    Ind_ROIArrayI = SLM_CA_ROIs{CurrentCell};
    Ind_EventArrayI = SLM_CA_EventSize{CurrentCell};

    if isempty(Ind_EventArrayII)
        Ind_EventArrayII = Ind_EventArrayI;
    else
        Ind_EventArrayII = cat(2, Ind_EventArrayII, Ind_EventArrayI);
    end
    
    if isempty(Ind_ROIArrayII)
        Ind_ROIArrayII = Ind_ROIArrayI;
    else
        Ind_ROIArrayII = cat(2, Ind_ROIArrayII, Ind_ROIArrayI);
    end

end

All_ROIs = zeros(xdim, ydim);
All_EventSizes = zeros(xdim, ydim);

All_ROIs(Ind_ROIArrayII) = 1;
All_EventSizes(Ind_EventArrayII) = 1;

if ~strcmp(DefString, 'NoSave')
     save(strcat('SLM_CA2Array_', DefString, '.mat'), 'All_ROIs', 'All_EventSizes');
end

if turnONimages == 1
    figure('Name', strcat('All established ROIs for "', DefString, '"')), imshow(All_ROIs);
    figure('Name', strcat('All model events for "', DefString, '"')), imshow(All_EventSizes);    
end 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    













