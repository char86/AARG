function ShowEventDurations(CA_dataoutput,strlists,suffixStrs)

% Outline 
% ShowEventDurations plots two panels of frequency distributions
% for events contained within CA_dataoutput. These allow the user to decide
% what Post-Detection Exclusion to use when AARG.m is called.

% Author: Charlie J Gilbride
% version: 1.0.0

PoolingLevel = 2;
%CG: PoolingLevel = 1: this means only split data files will be pooled
%ShowEventDurations will check the 'Condition' field in the 'dataoutput' structures 
%found in each cell. Data belonging to different conditions within the same
%experiment will be kept separate. 
%CG: PoolingLevel = 2: split data files will be pooled and data across
%different experiments (as long as the data is within the same condition)
%will be pooled.
%CG: PoolingLevel = 3: data across all cells in CA_dataoutput will be
%pooled. 
screensize = get( groot, 'Screensize' );

Sz_strlists = size(strlists); numdirs = Sz_strlists(1)*Sz_strlists(2);
for cDir = 1 : numdirs
    if cDir == 1; allNumStrFilesList = strlists(cDir).NumStrFilesList;
    else allNumStrFilesList = cat(2,allNumStrFilesList,strlists(cDir).NumStrFilesList);
    end
end

if PoolingLevel == 2;
    NumConditions = numel(unique(allNumStrFilesList));
    for cCond = 1 : NumConditions
        CondArray = allNumStrFilesList == cCond;
        idxloc = find(CondArray); NumFiles = numel(idxloc);
        pooled_TFD_allEvents = []; pooled_FAM_allEvents = [];
        for cFile = 1 : NumFiles
            dataoutput = CA_dataoutput{idxloc(cFile)};
            TFD_allEvents = dataoutput.TFD_allEvents; FAM_allEvents = dataoutput.FAM_allEvents;
            pooled_TFD_allEvents = [pooled_TFD_allEvents;TFD_allEvents];
            pooled_FAM_allEvents = [pooled_FAM_allEvents;FAM_allEvents];
        end
        figure('NumberTitle', 'Off', 'Name',strcat('Condition: "',suffixStrs{cCond},'"'),...
            'Position',screensize);
        hold on
        subplot(2,1,1);
        histogram(pooled_FAM_allEvents);
        title('Duration of detections between predicted detection peak and detection end',...
            'FontSize', 18);
        ax1 = gca;
        xlabel('Duration (frames)', 'FontSize', 14);
        ylabel('Number of events', 'FontSize', 14);

        subplot(2,1,2);
        histogram(pooled_TFD_allEvents);
        title('Duration of detections between detection start and detection end',...
            'FontSize', 18);
        ax2 = gca;
        xlabel('Duration (frames)', 'FontSize', 14);
        ylabel('Number of events', 'FontSize', 14);

        maxylim = 100;
        ax1xlim = ax1.XLim; ax2xlim = ax2.XLim;
        maxxlim = max([ax1xlim(2), ax2xlim(2)]); 

        ax1.YLim = [0, maxylim];
        ax2.YLim = [0, maxylim];
        ax1.XLim = [0, maxxlim]; 
        ax2.XLim = [0, maxxlim]; 
                        
    end
    
end
