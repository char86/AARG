function [CA_BinMedians] = AddMedianField(CA_Bins,suffixStrs,NumDirs)

%CG: CA_BinMedians will be the main input argument for GraphMedianData. 
NumConditions = size(suffixStrs,2); NumBins = size(CA_Bins,2); CA_BinMedians = {}; 
fh_hist = figure('Visible','Off');
for cConditionIdx = 1 : NumConditions
    medianData = struct;
    for cBinIdx = 1 : NumBins
        cBin = CA_Bins{cConditionIdx,cBinIdx}; 
        if ~isempty(cBin)
            NumROIs = size(cBin,2); contributions = [];
            for cROIIdx = 1 : NumROIs
                contributions = [contributions; cBin(cROIIdx).contribution];
            end

            for cDirIdx = 1 : NumDirs
                
                
                
                directorySelection = contributions(:,1);
                targetIndces = find(directorySelection == cDirIdx);
                ampsForCurrentDir = []; decaysForCurrentDir = []; distances = []; 
                branchContributionAmps = []; branchContributionDecays = []; 
                ROIContributionAmps = 0; ROIContributionDecays = 0;
                for cTargetIdx = 1 : numel(targetIndces)
                    ampsForCurrentDir = [ampsForCurrentDir, cBin(targetIndces(cTargetIdx)).amps];
                    decaysForCurrentDir = [decaysForCurrentDir, cBin(targetIndces(cTargetIdx)).decays];
                    distances = [distances, cBin(targetIndces(cTargetIdx)).distance];

                    if cDirIdx == 5
    %                     aa = contributions(targetIndces(cTargetIdx),2); disp(num2str(aa))
                        stophere = 1;
                    end

                    if sum(~isnan(cBin(targetIndces(cTargetIdx)).amps)) > 0
                        branchContributionAmps = [branchContributionAmps, contributions(targetIndces(cTargetIdx),2)];
                        ROIContributionAmps = ROIContributionAmps + 1;
                    end


                    if sum(~isnan(cBin(targetIndces(cTargetIdx)).decays)) > 0 && sum(cBin(targetIndces(cTargetIdx)).decays,'omitnan') < 5000000000000000
                        branchContributionDecays = [branchContributionDecays, contributions(targetIndces(cTargetIdx),2)];
                        ROIContributionDecays = ROIContributionDecays + 1;
                    end
                end
                
                decaysLogicFilter = decaysForCurrentDir == 5000000000000000;
                decaysForCurrentDir(decaysLogicFilter) = NaN;
                %CG: Some events with very slow decay times may have curve fits
                %with exponential decays that are extremely large. In such cases, the
                %most appropriate action would be to discard the event or
                %accept only its amplitude measurement. The decay of such
                %events may be distorted by a second, overlapping event, but
                %this may not have been apparent when classifying the events.
                %It seems any such events have decay constants of
                %"5000000000000000". 
                
%                 [medianAmp, medianDecay] = medianFromProbabilityDistribution(ampsForCurrentDir, decaysForCurrentDir, fh_hist);
                
                medianData(cDirIdx).amps = median(ampsForCurrentDir, 'omitnan');
%                 medianData(cDirIdx).amps = medianAmp;
                
                medianData(cDirIdx).decays = median(decaysForCurrentDir, 'omitnan');
%                 medianData(cDirIdx).decays = medianDecay;
                
                medianData(cDirIdx).distance = [min(distances) max(distances)];

                medianData(cDirIdx).branchContributionAmps = unique(branchContributionAmps);
                medianData(cDirIdx).branchContributionDecays = unique(branchContributionDecays);
                medianData(cDirIdx).ROIContributionAmps = ROIContributionAmps;
                medianData(cDirIdx).ROIContributionDecays = ROIContributionDecays;
            end
            CA_BinMedians{cConditionIdx, cBinIdx} = medianData;
        end
    end
end
        
        
