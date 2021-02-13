function [CA_Bones] = DiscardDeselectedROIs(skeleton, whiteROI, anotherROI, anotherColor)

RightWhiteLines = [];
CA_ROIc = skeleton.CA_ROIc; CA_id = skeleton.CA_id;
Sz_CA_ROIc = size(CA_ROIc); subsidROIList = skeleton.subsidROIList;
CA_StartPoints = skeleton.CA_StartPoints; Sz_CA_StartPoints = size(CA_StartPoints);

if ~isempty(CA_id)
    Sz_CA_id = size(CA_id); NumBranches = Sz_CA_id(1);
else
    NumBranches = 1;
end
RightROIidces = [];

for cBranch = 1 : NumBranches
    
    id_idx = CA_id{cBranch,1};
    
    for cROIidx = 1 : Sz_CA_ROIc(1)

        set(CA_ROIc{cROIidx, 2}, 'Parent', gca);

        if ~isempty(find(id_idx(:,1) == CA_ROIc{cROIidx, 3},1)) 
            if ~isempty(CA_ROIc{cROIidx, 6}) 

                    set(CA_ROIc{cROIidx, 6}, 'Parent', gca);
                    if isempty(RightWhiteLines) 
                        RightWhiteLines = CA_ROIc{cROIidx, 8};
                    elseif ~isempty(RightWhiteLines)
                        BoneIdces = CA_ROIc{cROIidx, 8};
                        RightWhiteLines(end+1) = BoneIdces(1);
                    end
                    set(CA_ROIc{cROIidx, 6}, 'Linewidth', 0.5);
                    set(CA_ROIc{cROIidx, 6}, 'Visible', 'on');
                    set(CA_ROIc{cROIidx, 2}, 'Visible', 'on');
                    if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                        set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 1 1]);
                    elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                        set(CA_ROIc{cROIidx, 2}, 'FaceColor', anotherColor);
                    else
                        set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 0 0]);
                    end
                        
                        
            else
%CG: if there is no connecting green line, then the ROI is excluded from
%analysis and should also be excluded from display.
                set(CA_ROIc{cROIidx, 2}, 'Visible', 'off');                    
            end
            if isempty(RightROIidces)
                RightROIidces = cROIidx;
            else
                RightROIidces(end+1) = cROIidx;
            end

        elseif isempty(find(id_idx(:,1) == CA_ROIc{cROIidx, 3},1)) &&...
                    isempty(find(RightROIidces == cROIidx,1))

                set(CA_ROIc{cROIidx, 2}, 'Visible', 'off');
                if ~isempty(find(whiteROI == CA_ROIc{cROIidx,3}, 1))
                    set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 1 1]);
                elseif ~isempty(find(anotherROI == CA_ROIc{cROIidx,3}, 1))
                    set(CA_ROIc{cROIidx, 2}, 'FaceColor', anotherColor);
                else
                    set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 0 0]);
                end
                if ~isempty(CA_ROIc{cROIidx, 6})
                    set(CA_ROIc{cROIidx, 6}, 'Visible', 'off');

                end
        end
        if ~isempty(subsidROIList)
            if ~isempty(find(subsidROIList(:,1) == CA_ROIc{cROIidx, 3},1))
                set(CA_ROIc{cROIidx, 2}, 'Visible', 'off');
                set(CA_ROIc{cROIidx, 6}, 'Visible', 'off');
            end
        end
        hold on
    end
end

hold off

if ~isempty(RightWhiteLines)
    CA_Bones = skeleton.CA_Bones; Sz_CA_Bones = size(CA_Bones);
    for cBoneIdx = 1 : Sz_CA_Bones(1)
        if isempty(find(RightWhiteLines == cBoneIdx,1))
            set(CA_Bones{cBoneIdx, 1}, 'Visible', 'Off')
        else
            set(CA_Bones{cBoneIdx, 1}, 'Linewidth', 0.5)
        end
        BoneID = ''; BoneLength = []; cROIidx = 1;
%CG: while loop finds dendrite classification. Must assume the user has
%drawn the white lines with this step in mind. Primary => white line is
%directly connected to a start point. Secondary => white line is connected
%to a start point through 1 other white line. >Secondary => white line is
%connected to start point through 2 or more other white lines. 
        while cROIidx < Sz_CA_ROIc(1)
            whitelines = CA_ROIc{cROIidx, 8};
            whitelines = fliplr(whitelines);
            IdxLoc = find(whitelines == cBoneIdx);

            blengths = CA_ROIc{cROIidx, 9};
            blengths = fliplr(blengths);
            if ~isempty(whitelines) && ~isempty(IdxLoc)
                if IdxLoc == 1
                    BoneID = 'Primary';
                    BoneLength = blengths;
                elseif IdxLoc == 2
                    BoneID = 'Secondary';
                    BoneLength = blengths(2) - blengths(1);
                elseif IdxLoc > 2
                    BoneID = '>Secondary';
                    BoneLength = blengths(3) - blengths(2);
                end
                cROIidx = Sz_CA_ROIc(1);
            end
            cROIidx = cROIidx + 1;
        end
        CA_Bones{cBoneIdx, 5} = BoneID;
        CA_Bones{cBoneIdx, 6} = BoneLength;
    end

    for cCell = 1 : Sz_CA_StartPoints(1)
        set(CA_StartPoints{cCell}, 'Visible', 'Off')
    end
end