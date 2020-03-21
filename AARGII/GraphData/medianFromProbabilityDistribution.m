function [medianAmp, medianDecay] = ...
    medianFromProbabilityDistribution(ampsForCurrentDir, decaysForCurrentDir, fh_hist)

if ~isempty(ampsForCurrentDir)
    %CG: collect the median amplitude...
    figure(fh_hist); fh_hist.Visible = 'Off';
    ampsNaNLogic = isnan(ampsForCurrentDir); ampsForCurrentDir(ampsNaNLogic) = [];
    h1 = histogram(ampsForCurrentDir, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
    h1.BinWidth = 5;

    amps_sorted = sort(ampsForCurrentDir, 'ascend');
    probabilityValues = h1.Values; 
    if ~isnan(probabilityValues)
        cumulativeProb = 0; valueCount = 0;
        zeroValues = probabilityValues == 0; probabilityValues_ZerosOut = probabilityValues; 
        probabilityValues_ZerosOut(zeroValues) = [];
        minProbValue = min(probabilityValues_ZerosOut);

        if numel(find(probabilityValues_ZerosOut == minProbValue)) <= 1
            errordlg('frequency of minProbValue unexpectedly low in the list of probability values. Contact AARG developers for assistance')
        end

        medianAmp = []; cIdx = 1;
        while cIdx <= numel(probabilityValues) && isempty(medianAmp)
            if probabilityValues(cIdx) ~= 0
                cumulativeProb = cumulativeProb + probabilityValues(cIdx);
                valueCount = valueCount + probabilityValues(cIdx)/minProbValue;
                if cumulativeProb > 0.5;
                    valueCount = valueCount - ((cumulativeProb-0.5)/minProbValue);
                    medianAmp = amps_sorted(round(valueCount));
                end
            end
            cIdx = cIdx + 1;
        end
    else
        medianAmp = NaN;
    end
elseif isempty(ampsForCurrentDir)
    medianAmp = NaN;
end

if ~isempty(decaysForCurrentDir)
    %CG: collect the median decay...

    figure(fh_hist); fh_hist.Visible = 'Off';
    decaysNaNLogic = isnan(decaysForCurrentDir); decaysForCurrentDir(decaysNaNLogic) = [];
    h1 = histogram(decaysForCurrentDir, 'DisplayStyle', 'stairs', 'Normalization', 'probability');
    h1.BinWidth = 0.01;

    decays_sorted = sort(decaysForCurrentDir, 'ascend');
    probabilityValues = h1.Values; 
    if ~isnan(probabilityValues)
        cumulativeProb = 0; valueCount = 0;
        zeroValues = probabilityValues == 0; probabilityValues_ZerosOut = probabilityValues; 
        probabilityValues_ZerosOut(zeroValues) = [];
        minProbValue = min(probabilityValues_ZerosOut);

        if numel(find(probabilityValues_ZerosOut == minProbValue)) <= 1
            errordlg('frequency of minProbValue unexpectedly low in the list of probability values. Contact AARG developers for assistance')
        end

        medianDecay = []; cIdx = 1;
        while cIdx <= numel(probabilityValues) && isempty(medianDecay)
            if probabilityValues(cIdx) ~= 0
                cumulativeProb = cumulativeProb + probabilityValues(cIdx);
                valueCount = valueCount + probabilityValues(cIdx)/minProbValue;
                if cumulativeProb > 0.5;
                    valueCount = valueCount - ((cumulativeProb-0.5)/minProbValue);
                    medianDecay = decays_sorted(round(valueCount));
                end
            end
            cIdx = cIdx + 1;
        end
    else
        medianDecay = NaN;
    end
elseif isempty(ampsForCurrentDir)
    medianDecay = NaN;
end