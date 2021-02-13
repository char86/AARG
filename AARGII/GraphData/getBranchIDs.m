function getBranchIDs(suffixStrs)

%CG: Outline: subroutine for graphAARG. Uses data stored in CA_ROIc to
%show outlines of dendrites drawn by the user (see AARGII). getBranchIDs
%labels each dendrite shown with its ID number. This number should then be
%entered appropriately into the 'Branch ID' text box. For more information
%look at the Help file for AARGII.

% Author: Charlie J Gilbride
% version: 1.0.0
datadir = uipickfiles('prompt','Select data folder to find branch numbers');
if nargin == 0
    suffixStrs = {'b_TTAP2'};
else
    
    suffixStrs = {suffixStrs{end}};
end

if ischar(datadir) || iscell(datadir)
        
    numberOfDirs = size(datadir,2);
    for cDirIdx = 1 : numberOfDirs
        cDir = datadir{cDirIdx}; [~,CellName,~] = fileparts(cDir); 
        TDFolders = {}; tdfc = 0; Sz_suffixStrs = size(suffixStrs); 

        for cCell = 1 : Sz_suffixStrs(2)
            cItem = strcat(CellName,suffixStrs{cCell},'_ThresholdData');
            tdfc = tdfc + 1; TDFolders{tdfc,1} = cItem;
            cd(strcat(cDir,'/', cItem))
            fmedian = []; cfile = 0; dirInfo = dir;
            while isempty(fmedian) 
                cfile = cfile + 1; currentfile = dirInfo(cfile).name;
                if ~isempty(strfind(currentfile, 'Chunk_1'))
                    data = load(currentfile); fmedian = data.fmedian;
                end
            end
        end

        cd(cDir); CellBonesFile = strcat(CellName,'_CellBonesCAs.mat');

        try testing = load(CellBonesFile); CBfilefound = 1;
        catch
            CBfilefound = 0;
        end

        if ischar(cDir) && CBfilefound == 1
            screensize = get( groot, 'Screensize' );

            logFluoImage=log10(fmedian); 
            fh = figure('Visible', 'Off', 'NumberTitle', 'Off', 'Name', 'Connecting ROIs to dendrites GUI');

            ax1 = axes; logFluoImage=logFluoImage-min(logFluoImage(:)); 
            logFluoImage=logFluoImage/max(logFluoImage(:));                             
            logFluoImage=uint8(logFluoImage*256);                                       
            cyanColorMap=([zeros(256,1),linspace(0,1,256)',linspace(0,1,256)']);
            left = screensize(3)*0.2; bot = screensize(2);
            wid = screensize(3)*0.66; hei =  screensize(4);
            colormap(cyanColorMap); fmedianRGB=ind2rgb(logFluoImage,cyanColorMap);
            imagesc(fmedianRGB);
            set(ax1, 'XTickLabel', '', 'TickLength', [0 0]);
            set(ax1, 'YTickLabel', '', 'TickLength', [0 0]);
            set(ax1, 'box', 'off'); set(fh, 'OuterPosition', [left, bot, wid, hei]);
            set(fh, 'NumberTitle', 'Off')
            axis image; 
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
                set(CA_Bones{cBoneIdx, 1}, 'Parent', ax1)
                set(CA_Bones{cBoneIdx, 1}, 'Visible', 'on')
                lineCoordinates = CA_Bones{cBoneIdx, 3}; lineLength = size(lineCoordinates,1);
        %CG: create arrow pointing perpendicularly to the current dendrite.
                x_endRaw = lineCoordinates(round(lineLength/2),1); y_endRaw = lineCoordinates(round(lineLength/2),2);
                phandle = patch(x_endRaw,y_endRaw,'y','LineStyle',':','Marker','.', 'MarkerSize', 24);
                phandle.MarkerFaceColor = 'y'; phandle.MarkerEdgeColor = 'y';
                hypotLength = 5*2.083; %CG: should be about 5µm.
                adjacentLength = cos(pi/8)*hypotLength; oppositeLength = sin(pi/8)*hypotLength;

                x_beginRaw = x_endRaw + adjacentLength; y_beginRaw = y_endRaw + oppositeLength;

                hold on 
                xx = [x_beginRaw x_endRaw]; yy = [y_beginRaw y_endRaw];
                plot(xx, yy, 'Color','y','LineWidth',1.5)
        %CG: label dendrite with its ID number.
                hText = text(xx(1),yy(1),strcat('Dendrite: "', num2str(cBoneIdx),'"'));
                hText.Color = 'y'; hText.FontSize = 12;

            end

            title(strcat('Experiment name: "',CellName,'"'),'FontSize', 14);

        elseif ischar(cDir) && CBfilefound == 0
            disp(strcat('No Bones file detected for experiment: "',CellName,'"'))
        end
    end
end