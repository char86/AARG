function [CA_id,CA_Distances,subsidROIList] = PickROIs(disbran,skeleton, whiteROI, anotherROI, anotherColor)

datadir = disbran.datadir; Resize = disbran.Resize; fmedianRGB = disbran.fmedianRGB;
DispOnlyProperROIs = 'Off';
condID = '';

if iscell(datadir)
    
    NumDirs = size(datadir,2);
    
    for cDirIdx = 1 : NumDirs
        cd(datadir{cDirIdx})
        dirInfo = dir;                                                                
        NumItems = size(dirInfo,1);
        [~, CellName, ~] = fileparts(datadir{cDirIdx});
        subsidROIList = [];

        for cItemIdx = 1 : NumItems

            cItem = dirInfo(cItemIdx).name;

            if ~isempty(strfind(cItem, '_fmedian_')) && ~isempty(strfind(cItem, '_Chunk_'))
                medaverage = load(cItem);
                fmedian = medaverage.MedianAverage_Section;
            elseif ~isempty(strfind(cItem, '_fmedian_')) && isempty(strfind(cItem, '_Chunk_'))
                medaverage = load(cItem, 'fmedian');
                fmedian = medaverage.fmedian;
            end
            if strcmp(DispOnlyProperROIs,'On') && ~isempty(strfind(cItem,strcat(condID,'_ThresholdData')))
                try
                    susi = load(strcat('CAs_CoreFourKR_', CellName, condID, '.mat'));
                catch
                    susi = load(strcat('CAs_CoreThreeKR_', CellName, condID, '.mat'));
                end
                subsidROIList = susi.subsidROIList;

            end
        end

        if strcmp(Resize,'On')
%             fmedian=myReSize3_bin(fmedian,3);
            disp(strcat('average projection image "fmedian" resized to 170x170')) 
        end
        
        [~,CellName,~] = fileparts(datadir{cDirIdx});
        CA_id = {}; id_idx = [];
        
        screensize = get( groot, 'Screensize' ); 
        fh = figure('Visible', 'On', 'NumberTitle', 'Off');

        left = screensize(3)*0.2;
        bot = screensize(2); wid = screensize(3)*0.66; hei =  screensize(4);
        imagesc(fmedianRGB);
        
        set(gca, 'XTickLabel', '', 'TickLength', [0 0]);
        set(gca, 'YTickLabel', '', 'TickLength', [0 0]);
        set(gca, 'box', 'off'); 
        set(fh, 'OuterPosition', [left, bot, wid, hei]);
        set(fh, 'NumberTitle', 'Off'); axis image; 
        title(strcat(CellName, '; Follow appropriate intructions in MATLAB command line'));
        set(fh, 'Visible', 'On');

        CA_StartPoints = skeleton.CA_StartPoints; Sz_CA_StartPoints = size(CA_StartPoints);
        for cCell = 1 : Sz_CA_StartPoints(1)
            set(CA_StartPoints{cCell}, 'Parent', gca)
            set(CA_StartPoints{cCell}, 'Visible', 'On')
        end
        
        CA_Bones = skeleton.CA_Bones; Sz_CA_Bones = size(CA_Bones);
        for cBoneIdx = 1 : Sz_CA_Bones(1)
            set(CA_Bones{cBoneIdx, 1}, 'Parent', gca)
            set(CA_Bones{cBoneIdx, 1}, 'Visible', 'On')
            set(CA_Bones{cBoneIdx, 1}, 'Color', 'w')
        end

        RightWhiteLines = [];
        CA_ROIc = skeleton.CA_ROIc; Sz_CA_ROIc = size(CA_ROIc);
        
        for cROIidx = 1 : Sz_CA_ROIc(1)
%CG: for-loop resets every connected ROI to red and hides all unconnected
%ROIs. 
            set(CA_ROIc{cROIidx, 2}, 'Parent', gca);
            if ~isempty(CA_ROIc{cROIidx, 6})
                set(CA_ROIc{cROIidx, 6}, 'Parent', gca);
            end
%             if ~isempty(find(whiteROI == CA_ROIc{cROIidx, 3},1))
%                 set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 1 1]);
%                 set(CA_ROIc{cROIidx, 2}, 'EdgeColor', [1 1 1]);
%             end
            if ~isempty(subsidROIList)
%CG: necessary? subsidary ROIs should already be deleted from CA_ROIs and
%should therefore not pass to CA_ROIc. 
                if ~isempty(find(subsidROIList(:,1) == CA_ROIc{cROIidx, 3},1))
                    set(CA_ROIc{cROIidx, 2}, 'Visible', 'off');
                    set(CA_ROIc{cROIidx, 6}, 'Visible', 'off');
                end
            end
            hold on
        end
        
        CA_Distances = skeleton.CA_Distances; k = [];
        
        title(sprintf(strcat('Instructions: (Step 4) Left-click on red ROIs to keep in analysis.\n',...
        'Left-click on blue ROIs to de-select ROIs for analysis.\n',...
        'All ROIs marked in blue will be included.\n',...
        'Press "a" to accept and finish.\n')), 'FontSize', 16)
        
        while isempty(k) 

            k = waitforbuttonpress;

            if k == 1
                KeyPressRes = get(fh,'CurrentCharacter');
                if ~strcmp(KeyPressRes, 'a') && ~strcmp(KeyPressRes, 'z') && ~strcmp(KeyPressRes, 'm')
                    k = [];
                elseif strcmp(KeyPressRes, 'z')
                    if ~isempty(CA_id)
                        Sz_CA_id = size(CA_id);

                        for cCellidx = 1 : Sz_CA_id(1)
                            id_idx = CA_id{cCellidx};
                            Sz_id_idx = size(id_idx);
                            for cidx = 1 : Sz_id_idx(1)
                                cvec_idx = id_idx(cidx,2);
                                if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                                    set(CA_ROIc{cvec_idx,2}, 'FaceColor', [1 1 1])
                                elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                                    set(CA_ROIc{cvec_idx,2}, 'FaceColor', anotherColor)
                                else
                                    set(CA_ROIc{cvec_idx,2}, 'FaceColor', [1 0 0])
                                end
                            end
                        end
                    else
                        disp('nothing to undo')
                    end
                    k = []; id_idx = []; CA_id = {};
                elseif strcmp(KeyPressRes, 'm')
%CG: In this case, all blue ROIs will be merged together as if they were on
%one path. 
                    disp('ROIs merged into single path...')
                    disp('before...')
                    CA_id
                    id_idx = [];
                    if ~isempty(CA_id)
                        Sz_CA_id = size(CA_id);
                        for cellidx = 1 : Sz_CA_id(1)

                            if isempty(id_idx)
                                id_idx = CA_id{cellidx};
                            else
                                id_idx_holder = CA_id{cellidx};
                                id_idx = cat(1, id_idx, id_idx_holder);
                            end

                        end
                    end
                    CA_id = {};
                    disp('after...')
                    CA_id{1,1} = id_idx;
                    CA_id
                    k = [];
                end
            elseif k == 0

                mouseclick = get(fh, 'SelectionType'); PointerStatus = get(fh, 'Pointer');

                if strcmp(mouseclick, 'normal') && strcmp(PointerStatus, 'arrow')

                    cObj = gco; ClassDef = class(cObj);
                    if strcmp(ClassDef, 'matlab.graphics.primitive.Patch')

                        fc = get(cObj, 'FaceColor'); fc_red = fc - [1 0 0];
                        fc_white = fc - [1 1 1]; fc_anotherColor = fc - anotherColor;
                        fc_blue = fc - [0 0 1]; cobjx = get(cObj, 'XData');
                        cobjy = get(cObj, 'YData');

                        if (sum(fc_red) == 0 && isempty(find(fc_red,1))) || ...
                                (sum(fc_white) == 0 && isempty(find(fc_white,1))) || ...
                                    (sum(fc_anotherColor) == 0 && isempty(find(fc_anotherColor,1)))
%CG: if the color of the clicked ROI is red. 
%CG: get all white lines lying between the current ROI and the connecting
%start point.
                            cROIidx = 0;
                            while cROIidx < Sz_CA_ROIc(1)
%CG: find the ROI that has been clicked.
                                cROIidx = cROIidx + 1;
                                cROIx = get(CA_ROIc{cROIidx,2}, 'XData');
                                cROIy = get(CA_ROIc{cROIidx,2}, 'YData');
                                if sum(cobjx - cROIx) == 0 && sum(cobjy - cROIy) == 0
                                    cROI = cROIidx; cROIidx = Sz_CA_ROIc(1);
                                end
                            end

                            id_idx(1,1) = CA_ROIc{cROI,3}; id_idx(1,2) = cROI;
                            id_idx(1,3) = CA_ROIc{cROI,10}; bones = CA_ROIc{cROI,8};
                            id_idx(1,4) = bones(1); set(CA_ROIc{cROI, 2}, 'FaceColor', [0 0 1]);

                            ROIidces = cell2mat(CA_ROIc(:, 3));
                            BoneIdces = CA_ROIc{cROI, 8};
                            BoneIdces = BoneIdces(1);
%CG: Critical DistanceLimits_branches modification: BoneIdces = BoneIdces(1); 
%This means only ROIs connected to the same branch as the selected ROI will
%be included in the same activity map. 
                            iDistances = CA_ROIc{cROI, 9};
                            Sz_BoneIdces = size(BoneIdces);
                            for Idx = 1 : Sz_BoneIdces(2)
                                cBoneROIs = CA_Bones{BoneIdces(Idx), 4};
                                Sz_cBoneROIs = size(cBoneROIs);
                                c_distance = iDistances(Idx);
%CG: all ROIs on the current white line will be set against the value c_distance.
%Any that have a value smaller than this will be included in the heatmap
%display/further analyses. 
%CG: get ROIs closer to the start point and sitting along the same white line
%as the current ROI.                        
                                for cROIidx = 1 : Sz_cBoneROIs(2)
                                    caROIidx = find(ROIidces == CA_ROIc{cBoneROIs(cROIidx),3});
                                    cdis = CA_ROIc{caROIidx, 10};
                                    if cdis < c_distance
                                        set(CA_ROIc{caROIidx, 2}, 'FaceColor', [0 0 1]);
                                        id_idx(end+1,1) = CA_ROIc{caROIidx,3};
                                        id_idx(end,2) = caROIidx;
                                        id_idx(end,3) = CA_ROIc{caROIidx,10};
                                        bones = CA_ROIc{caROIidx,8};
                                        id_idx(end,4) = bones(1);
                                    end
                                end
                            end
                            if isempty(CA_id)
                                CA_id{1,1} = id_idx;
                            else
                                Sz_CA_id = size(CA_id);
                                StartNew = 1; 
%CG: it sometimes happens that ROIs connected to the same white line are
%missed and the user has to select them individually. The for loop checks
%if this selected ROI (which will be in the top row of id_idx) already
%shares a white line with ROIs stored in pre-existing cells in CA_id.
%Logically, if the individually selected ROI connects to the same white
%line as ROIs already stored in CA_id, it should be added to the contents
%of this cell.
                                for cCell = 1 : Sz_CA_id(1)
                                    pre_id_idx = CA_id{cCell};
                                    numberOfROIs = size(id_idx,1);
                                    new_id_idx = pre_id_idx;
                                    for cROIIdx = 1 : numberOfROIs
                                        if isempty(find(pre_id_idx(:,1) == id_idx(cROIIdx,1),1)) &&...
                                                ~isempty(find(pre_id_idx(:,4) == id_idx(cROIIdx,4),1))
                                            new_id_idx = cat(1,new_id_idx,id_idx(cROIIdx,:));
                                            StartNew = 0;
%                                             if pre_id_idx(1,4) == id_idx(1,4)
%                                                id_idx = cat(1,id_idx(1,:),pre_id_idx);
%                                                CA_id{cCell} = id_idx;
% 
% 
%                                                StartNew = 0;
%                                             end
                                        end
                                    end
                                    if StartNew == 0; CA_id{cCell} = new_id_idx; end
                                end
                                if StartNew == 1
                                    CA_id{end+1,1} = id_idx;
                                end
                            end
                            id_idx = [];

                        elseif sum(fc_blue) == 0 && isempty(find(fc_blue,1))

                            Sz_CA_id = size(CA_id);
                            arraysize = Sz_CA_id(1)*Sz_CA_id(2);
                            whitelineNum = []; cROIidx = 0;
                            while isempty(whitelineNum) && cROIidx < Sz_CA_ROIc(1)
                                cROIidx = cROIidx + 1;
                                xes = get(CA_ROIc{cROIidx, 2}, 'XData');
                                ys = get(CA_ROIc{cROIidx, 2}, 'YData');

                                if sum(cobjx - xes) == 0 && sum(cobjy - ys) == 0

                                    whitelineNum = CA_ROIc{cROIidx, 8};
                                    ROIsOnLine = CA_Bones{whitelineNum,4};

                                    NumROIsOnLine = numel(ROIsOnLine);
                                    for cROIon = 1 : NumROIsOnLine
                                        
                                        if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                                            set(CA_ROIc{ROIsOnLine(cROIon),2}, 'FaceColor', [1 1 1])
                                        elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                                            set(CA_ROIc{ROIsOnLine(cROIon),2}, 'FaceColor', anotherColor)
                                        else
                                            set(CA_ROIc{ROIsOnLine(cROIon),2}, 'FaceColor', [1 0 0])
                                        end
                                    end
                                end
                            end
                            whitelinefound = [];linecount = 0;
                            while isempty(whitelinefound) && linecount < Sz_CA_id(1) 
                                linecount = linecount + 1;
                                id_idx = CA_id{linecount, 1};
                                if whitelineNum(1) == id_idx(1,4)
                                    whitelinefound = 1;
                                end
                            end

                            id_idx = CA_id{linecount};
                            Sz_id_idx = size(id_idx);
                            for cidx = 1 : Sz_id_idx(1)
                                cvec_idx = id_idx(cidx,2);
                                if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                                    set(CA_ROIc{cvec_idx,2}, 'FaceColor', [1 1 1])
                                elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                                    set(CA_ROIc{cvec_idx,2}, 'FaceColor', anotherColor)
                                else
                                    set(CA_ROIc{cvec_idx,2}, 'FaceColor', [1 0 0])
                                end
                                    
                            end

                            if arraysize == 1
                                CA_id = {};
                            else
                                CA_id(linecount,:) = [];
                            end
%CG: the previous step might have turned blue labelled ROIs red even though
%they should still be blue labelled because they appear in some of the
%vectors contained in the remaining cells of CA_id.

                            if ~isempty(CA_id)
                                Sz_CA_id = size(CA_id);

                                for cCellidx = 1 : Sz_CA_id(1)
                                    id_idx = CA_id{cCellidx};
                                    Sz_id_idx = size(id_idx);
                                    for cidx = 1 : Sz_id_idx(1)
                                        cvec_idx = id_idx(cidx,2);
                                        set(CA_ROIc{cvec_idx,2}, 'FaceColor', [0 0 1])
                                    end
                                end

                            end
                            id_idx = [];

                        else
                            disp('left-click on ROIs')
                        end
                    else
                        disp('left-click on ROIs')
                    end
                elseif strcmp(PointerStatus, 'arrow')
                    disp('left-click on ROIs')
                end
                k = [];

            end

        end

        for cCell = 1 : Sz_CA_StartPoints(1)
            set(CA_StartPoints{cCell}, 'Visible', 'Off')
        end

        RightWhiteLines = [];
        Sz_CA_id = size(CA_id);
        id_idx = [];
        for cCell = 1 : Sz_CA_id(1)
            if isempty(id_idx)
                id_idx = CA_id{cCell,1};
            else
                id_idx = cat(1,id_idx,CA_id{cCell,1});
            end
        end

        for cROIidx = 1 : Sz_CA_ROIc(1)
            set(CA_ROIc{cROIidx, 2}, 'Parent', gca);
            if ~isempty(id_idx)
                if ~isempty(find(id_idx(:,1) == CA_ROIc{cROIidx, 3},1))
                    
                    if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                        set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 1 1]);
                    elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                        set(CA_ROIc{cROIidx, 2}, 'FaceColor', anotherColor);
                    else
                        set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 0 0]);
                    end
                    
                    if ~isempty(CA_ROIc{cROIidx, 6})
                        set(CA_ROIc{cROIidx, 6}, 'Parent', gca);
                    end
                    wBoneIdces = CA_ROIc{cROIidx, 8};
                    if isempty(RightWhiteLines)
                        RightWhiteLines = wBoneIdces(1);
                    else
                        RightWhiteLines(end+1) = wBoneIdces(1);
                    end
                    set(CA_ROIc{cROIidx, 6}, 'Linewidth', 0.5);
                elseif isempty(find(id_idx(:,1) == CA_ROIc{cROIidx, 3},1))
                    set(CA_ROIc{cROIidx, 2}, 'Visible', 'off');
                    if ~isempty(CA_ROIc{cROIidx, 6})
                        set(CA_ROIc{cROIidx, 6}, 'Visible', 'off');
                    end
                end
            elseif isempty(id_idx)
                %CG: if isempty(id_idx), then this probably means the user
                %has pressed 'a' without clicking any ROIs. This is either
                %an error or meeans the user is doing some other
                %unanticipated activity. In this case, the final image will
                %just show the cell with the dendrites traced and no ROIs
                %or green lines.
                set(CA_ROIc{cROIidx, 2}, 'Visible', 'off');
                if ~isempty(CA_ROIc{cROIidx, 6})
                    set(CA_ROIc{cROIidx, 6}, 'Visible', 'off');
                end
            end
        end

        Sz_CA_Bones = size(CA_Bones);
        for cBoneIdx = 1 : Sz_CA_Bones(1)
            if isempty(find(RightWhiteLines == cBoneIdx,1))
                set(CA_Bones{cBoneIdx, 1}, 'Visible', 'Off')
            else
                set(CA_Bones{cBoneIdx, 1}, 'Linewidth', 0.5)
            end
        end 
    end
end