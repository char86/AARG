function GraphMedianData(CA_BinMedians,suffixStrs,datadirs,SizeBins,configStruct)

% Author: Charlie J Gilbride
% version: 1.0.0

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
    
    yAxisLabel = ''; conditionString = strcat('. Condition: "', suffixStrs{cCondIdx},'"');
    
    if configStruct.displayPeakAmplitude == 1
        fh_median = figure('Name',strcat('Median amplitudes separated according to distance of events from soma',conditionString));
        yAxisLabel = 'Peak amplitude (a.u.)';
    elseif configStruct.displayDecayTimeConstant == 1
        fh_median = figure('Name',strcat('Median event decay time constants separated according to distance from soma', conditionString));
        yAxisLabel = 'tau decay constant (s)';
    end
    xAxisLabel = 'Condition';
    fh_median.NumberTitle = 'off'; fh_median.Position = screenSize.*positionScaled; 
    CA_GreyLineData = {};
    
    set(findall(fh_median, '-property', 'Units'), 'Units', 'Normalized')
    set(fh_median, 'Resize','on') 
    
    for cCondIdx = 1 : NumConditions
        if isempty(configStruct.maximumDistance)
            NumBins = size(CA_BinMedians,2); 
        elseif ~isempty(configStruct.maximumDistance)
            NumBins = size(CA_BinMedians,2); cBinIdx = 1;
            %CG: if maximum distance is non-empty, then the user has
            %specified that s/he only wants to see graphs with distances up
            %to a value of configStruct.maximumDistance
            while cBinIdx <= NumBins
                numberOfStructs = size(CA_BinMedians{cCondIdx,cBinIdx},2);
                dataStruct = CA_BinMedians{cCondIdx,cBinIdx}; structIdx = 1;
                while structIdx <= numberOfStructs 
                    maxMinDistances = dataStruct(structIdx).distance;
                    if ~isempty(maxMinDistances)
                        %CG: maxMinDistances can be empty if the current
                        %directory has no ROIs within the distance range of
                        %the current struct. 
                        if maxMinDistances(2) > configStruct.maximumDistance
                            NumBins = cBinIdx-1; structIdx = numberOfStructs;
                        end
                    end
                    structIdx = structIdx + 1;
                end
                cBinIdx = cBinIdx + 1;
            end
        end
                
        ansx = NumBins/2; 
        if cCondIdx == 1
            baselineData = [];
        end
            
        for cBinIdx = 1 : NumBins
                
            titleStr = strcat(num2str(SizeBins(cBinIdx)),...
                '-',num2str(SizeBins(cBinIdx+1)),'µm');
%CG: CONTINUE FROM HERE....3.3.19    
            binstr = CA_BinMedians{cCondIdx,cBinIdx}; NumDirs = size(binstr,2); 
            amps = []; decays = [];
%             if configStruct.displayPeakAmplitude == 1; amps = []; end
%             if configStruct.displayDecayTimeConstant == 1; decays = []; end
            if configStruct.displayNValues == 1; Nbranches = []; NROIs = []; end

            for cDirIdx = 1 : NumDirs
                if configStruct.displayPeakAmplitude == 1;
                    amps = [amps,binstr(cDirIdx).amps];
                    Nbranches = [Nbranches,binstr(cDirIdx).branchContributionAmps];
                    NROIs = [NROIs, binstr(cDirIdx).ROIContributionAmps];
                end
                if configStruct.displayDecayTimeConstant == 1;
                    decays = [decays,binstr(cDirIdx).decays];
                    Nbranches = [Nbranches,binstr(cDirIdx).branchContributionDecays];
                    NROIs = [NROIs, binstr(cDirIdx).ROIContributionDecays];
                end

            end
            if cCondIdx == 1
                baselineData = [baselineData, amps'];
            end
            Nbranches = numel(Nbranches);
            if ~isempty(amps) && configStruct.displayPeakAmplitude == 1
                Ndirs = numel(find(~isnan(amps)));
            elseif ~isempty(decays) && configStruct.displayDecayTimeConstant == 1
                Ndirs = numel(find(~isnan(decays)));
            end
            
%CG: plot data
            notEnoughData = 0;
            if configStruct.displayPeakAmplitude == 1
                if sum(~isnan(amps)) == 0; notEnoughData = 1; end
            end
            if configStruct.displayDecayTimeConstant == 1;
                if sum(~isnan(decays)) == 0; notEnoughData = 1; end
            end
                
            if notEnoughData == 0
                figure(fh_median); subplot(2,ceil(ansx),cBinIdx);
                if cCondIdx == 1
                    dataStruct = struct;
                    if configStruct.displayPeakAmplitude == 1; 
                        xx = ones(1,numel(amps));
                        p1 = plot(xx, amps, 'o');
                        dataStruct(1).yValues = amps;
                    end
                    if configStruct.displayDecayTimeConstant == 1; 
                        xx = ones(1,numel(decays));
                        p1 = plot(xx, decays, 'o');
                        dataStruct(1).yValues = decays;
                    end
                    dataStruct(1).xValues = xx;
                    
                elseif cCondIdx > 1
                    dataStruct = CA_GreyLineData{cBinIdx};
                    ax1 = gca; hold(ax1,'on');
                    
                    if configStruct.displayPeakAmplitude == 1; 
                        xx = ones(1,numel(amps)); xx = xx.*2;
                        p1 = plot(xx, amps, 'o','Color','r');
                        dataStruct(cCondIdx).yValues = amps;
                    end
                    if configStruct.displayDecayTimeConstant == 1; 
                        xx = ones(1,numel(decays)); xx = xx.*2;
                        p1 = plot(xx, decays, 'o','Color','r');
                        dataStruct(cCondIdx).yValues = decays;
                    end
                    dataStruct(cCondIdx).xValues = xx;

                    title(titleStr)
                    if cBinIdx == 1 || cBinIdx == (NumBins/2)+1; 
                        ylabel(yAxisLabel);
                    else
                        ylabel('')
                        if configStruct.hideUnnecessaryTickLabels == 1
                            ax1.YTickLabel = ''; %ax1.TickLength = [0 0];
                        end
                    end
                    if cBinIdx > (NumBins/2); 
                        ax1.XTickLabel = legendLabels; 
                    else 
                        xlabel(''); 
                        if configStruct.hideUnnecessaryTickLabels == 1
                             ax1.XTickLabel = ''; %ax1.TickLength = [0 0];
                        else
                            ax1.XTickLabel = legendLabels; 
                        end
                    end

                end
                CA_GreyLineData{cBinIdx} = dataStruct;
                ax1.XTick = [1 2];
                ax1.TickDir = 'out'; ax1.LineWidth = 2; 

                ax1 = gca; ax1.Box = 'off';
                ax1.FontSize = fontSize;

                if ~ischar(configStruct.xAxisMaxMedian)
                    ax1.XLim = [configStruct.xAxisMinMedian configStruct.xAxisMaxMedian];
                else
                    if strcmp(configStruct.xAxisMaxMedian, '?') && ~ischar(configStruct.xAxisMinMedian)
                        ax1.XLim = [0 3];
                    end
                end 
                if ~ischar(configStruct.yAxisMaxMedian)
                    ax1.YLim = [configStruct.yAxisMinMedian configStruct.yAxisMaxMedian];
                else
                    if strcmp(configStruct.yAxisMaxMedian, '?') && ~ischar(configStruct.yAxisMinMedian)
                        ax1.YLim = [0 800];

                    end
                end
                
            elseif notEnoughData == 1
                disp('Not enough data error: ned error2')
            end
        end
    end
%CG: now connect the markers with a grey line to indicate which condition belongs to a particular experiment.   

    for cBinIdx = 1 : NumBins
        figure(fh_median); subplot(2,ceil(ansx),cBinIdx);
        ax1 = gca; hold(ax1,'on');
        
        dataStruct = CA_GreyLineData{cBinIdx};
        if size(dataStruct,2) == NumConditions
            xx = []; yy = [];
            for cConditionIdx = 1 : NumConditions
                xx = [xx;dataStruct(cConditionIdx).xValues];
                yy = [yy;dataStruct(cConditionIdx).yValues];
            end

            [~,pp,~,stats] = ttest(yy(1,:), yy(2,:));

            xlims1 = (ax1.XLim); ylims1 = (ax1.YLim);

            labelTStat = strcat('T = "',num2str(stats.tstat),'"');
            t1xy = [xlims1(2)*0.4 ylims1(2)*0.9];
            text(t1xy(1),t1xy(2),labelTStat, 'FontSize', fontSize)

            labelPValue = strcat('P = "',num2str(pp),'"');
            t1xy = [xlims1(2)*0.4 ylims1(2)*0.8];
            text(t1xy(1),t1xy(2),labelPValue, 'FontSize', fontSize)

            NumDataPoints = size(xx,2);
            for cDataPointIdx = 1 : NumDataPoints
                plot(xx(:,cDataPointIdx), yy(:,cDataPointIdx), '-', 'Color', [0.75 0.75 0.75])
            end
        elseif size(dataStruct,2) < NumConditions
            disp('Not enough data error: ned error2')
        elseif size(dataStruct,2) > NumConditions
            disp('Unknown error: unknown error1')
        end
        
    end
   
end