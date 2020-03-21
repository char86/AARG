function [CA_dataoutput] = EventDurations(ThresholdDataFolderList,strlists)

%CG: collects statistics on the duration of events detected by
%DetectEvents. This includes data on the time between the local maxima of
%event spread and the end of an event. 

NumDataFolders = size(ThresholdDataFolderList,1); NumExpFolders = size(ThresholdDataFolderList,2);
CA_dataoutput = {};
%CG: turn off waitbar in spatial filtering function.
wbstr = 0; wb = waitbar(0,'please wait...'); set(wb,'Name','Measuring event durations');
wb_counter = 0;
for cExp = 1 : NumExpFolders
    NumStrFilesList = strlists(cExp).NumStrFilesList;
    for cDF = 1 : NumDataFolders
        dataoutput = struct;
        Condition = NumStrFilesList(cDF,1);
        cd(ThresholdDataFolderList{cDF, cExp});
        dirInfo = dir; numfiles = size(dirInfo,1); 

        dataoutput.ThresholdDataFolder = ThresholdDataFolderList{cDF, cExp};
        dataoutput.Condition = Condition; TFD_allEvents = []; FAM_allEvents = [];
        for cFile = 1 : numfiles
            CurrentItem = dirInfo(cFile).name;
            if numel(CurrentItem)>2
                if ~isempty(strfind(CurrentItem, '_Chunk_'))
                    data = load(CurrentItem);
                    ChunkSize = data.ChunkSize; fmedian = data.fmedian;
                    suprathres_Evals = data.suprathres_Evals;
                    xdim = size(fmedian,1); ydim = size(fmedian,2);
                    blk = zeros(xdim, ydim, ChunkSize);
                    NumEvents = size(suprathres_Evals,2); iix = [];
%CG: put all events into empty matrix 'blk'.
                    for cpp = 1 : NumEvents
                        blk(suprathres_Evals{cpp}) = 1; 
                        ii = suprathres_Evals{cpp}; iix=[iix, ii'];
                    end
%CG: use the same calculation to find local maxima of event spread.
                    ii=2.7; 
                    try 
%CG: inconsistent case across operating systems. myAVG3 for windows and
%myAvg3 for macs
                        dd=myAVG3(blk,ii);
                    catch
                        dd=myAvg3(blk,ii);
                    end
                    dd = mySmooth2D_All3_fft(dd,5.4,5.4,wbstr); 
                    ccx = blk*0;ccx(iix) = 1; dd = dd.*ccx; MaximaMat = myLocalMax3D2(dd,2);
                    blk2 = blk.*0;
                    for cpp = 1 : NumEvents
                        cEvent = suprathres_Evals{cpp};
                        [~, maxind_lmax] = max(MaximaMat(cEvent));
                        [~,~,zz] = ind2sub([xdim,ydim,ChunkSize], cEvent);
                        [~,~,zz_max] = ind2sub([xdim,ydim,ChunkSize], cEvent(maxind_lmax));
                        unique_zz = unique(zz); TotalFrameDuration = numel(unique_zz);
                        FramesAfterMax = unique_zz(end) - zz_max;
                        TFD_allEvents = [TFD_allEvents; TotalFrameDuration];
                        FAM_allEvents = [FAM_allEvents; FramesAfterMax];
                    end
                end
            end
            wb_counter = wb_counter + 1;
            PercentIncr = round((wb_counter/(NumExpFolders*NumDataFolders*numfiles)*100));
            waitbar(PercentIncr/100,wb)
        end
        dataoutput.TFD_allEvents = TFD_allEvents; dataoutput.FAM_allEvents = FAM_allEvents;
        CA_dataoutput{cDF,cExp} = dataoutput;
    end
end
close(wb)
