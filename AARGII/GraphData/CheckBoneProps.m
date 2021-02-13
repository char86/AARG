function CheckBoneProps

%Outline CheckBoneProps: displays rois and connecting lines stored in
%'..._CellBonesCAs.mat' file. This is helpful when the user wants to
%know the distance of start point from the cell body when the cell body is
%not part of the timelapse image. 

% Author: Charlie J Gilbride
% version: 1.0.0

datadir = uigetdir('Select data folder for CellBones');
suffixStrs = {'b_TTAP2'};

if ischar(datadir)
        
    [~,CellName,~] = fileparts(datadir); TDFolders = {};tdfc = 0;
    Sz_suffixStrs = size(suffixStrs); CellBonesFile = strcat(CellName,'_','CellBonesCAs.mat');
    
    for cCell = 1 : Sz_suffixStrs(2)
        cItem = strcat(CellName,suffixStrs{cCell},'_ThresholdData');
        tdfc = tdfc + 1; TDFolders{tdfc,1} = cItem;
        cd(strcat(datadir,'/', cItem))
        fmedian = []; cfile = 0; dirInfo = dir;
        while isempty(fmedian) 
            cfile = cfile + 1; currentfile = dirInfo(cfile).name;
            if ~isempty(strfind(currentfile, 'Chunk_1'))
                data = load(currentfile); fmedian = data.fmedian;
            end
        end
    end
    
    cd(datadir); CellBonesFile = strcat(CellName,'_','CellBonesCAs.mat');
    
    try testing = load(CellBonesFile); CBfilefound = 1;
    catch
        CBfilefound = 0;
    end
end

if ischar(datadir) && CBfilefound == 1
    screensize = get( groot, 'Screensize' );
    
    logFluoImage=log10(fmedian); [xdim, ydim]=size(fmedian);
%     logFluoImage_copy = logFluoImage;
    fh = figure('Visible', 'Off', 'NumberTitle', 'Off', 'Name', 'Connecting ROIs to dendrites GUI');
    
    ax1 = axes; logFluoImage=logFluoImage-min(logFluoImage(:)); 
    logFluoImage=logFluoImage/max(logFluoImage(:));                             
    logFluoImage=uint8(logFluoImage*256);                                       
    cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);
    left = screensize(3)*0.2; bot = screensize(2);
    wid = screensize(3)*0.66; hei =  screensize(4);
    colormap(cyanColorMap); fmedianRGB=ind2rgb(logFluoImage,cyanColorMap);
    xAxisLimit = size(fmedianRGB,1); yAxisLimit = size(fmedianRGB,2);  
    im1 = imagesc(fmedianRGB);
    set(ax1, 'XTickLabel', '', 'TickLength', [0 0]);
    set(ax1, 'YTickLabel', '', 'TickLength', [0 0]);
    set(ax1, 'box', 'off'); set(fh, 'OuterPosition', [left, bot, wid, hei]);
    set(fh, 'NumberTitle', 'Off')
    %axis image; 
    set(fh, 'Visible', 'On');
    
    skeleton = load(CellBonesFile);
    CA_StartPoints = skeleton.CA_StartPoints; 
    Sz_CA_StartPoints = size(CA_StartPoints);
    for cCell = 1 : Sz_CA_StartPoints(1)
        set(CA_StartPoints{cCell}, 'Parent', ax1)
        set(CA_StartPoints{cCell}, 'Visible', 'On')
    end
    CA_Bones = skeleton.CA_Bones; 
    Sz_CA_Bones = size(CA_Bones);
    for cBoneIdx = 1 : Sz_CA_Bones(1)
        disp(cBoneIdx)
        set(CA_Bones{cBoneIdx, 1}, 'Parent', ax1)
        set(CA_Bones{cBoneIdx, 1}, 'Visible', 'on')
        lineCoordinates = CA_Bones{cBoneIdx, 3}; lineLength = size(lineCoordinates,1);
        x_endRaw = lineCoordinates(round(lineLength/2),1); y_endRaw = lineCoordinates(round(lineLength/2),2);
        phandle = patch(x_endRaw,y_endRaw,'y','LineStyle',':','Marker','.', 'MarkerSize', 24);
        phandle.MarkerFaceColor = 'y'; phandle.MarkerEdgeColor = 'y';
        hypotLength = 5*2.083; %CG: should be about 5µm.
        adjacentLength = cos(pi/8)*hypotLength; oppositeLength = sin(pi/8)*hypotLength;
        
        x_beginRaw = x_endRaw + adjacentLength; y_beginRaw = y_endRaw + oppositeLength;
        xMinMax = get(gca,'XLim'); yMinMax = get(gca,'YLim');
        
        hold on 
        xx = [x_beginRaw x_endRaw]; yy = [y_beginRaw y_endRaw];
        plot(xx, yy, 'Color','y','LineWidth',1.5)
        hText = text(xx(1),yy(1),strcat('Dendrite: "', num2str(cBoneIdx),'"'));
        hText.Color = 'y'; hText.FontSize = 12;
        
        
        stophere = 1;
%         set(CA_Bones{cBoneIdx, 1}, 'LineWidth', 2)
%         set(CA_Bones{cBoneIdx, 1}, 'Color', 'w')
    end
    CA_ROIc = skeleton.CA_ROIc; Sz_CA_ROIc = size(CA_ROIc);
    for cROIidx = 1 : Sz_CA_ROIc(1)
        fc = get(CA_ROIc{cROIidx, 2},'FaceColor');
        if sum(fc - [1 0 0]) == 0
            set(CA_ROIc{cROIidx, 2}, 'Parent', ax1); %set(CA_ROIc{cROIidx, 2}, 'FaceColor', [1 0 0])
            set(CA_ROIc{cROIidx, 6}, 'Parent', ax1);
        end
    end
    
    title(CellName)
    
elseif ischar(datadir) && CBfilefound == 0
    disp('No Bones file detected')
end