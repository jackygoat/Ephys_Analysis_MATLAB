clear;

mainDir = uigetdir();
cd(mainDir)
subj = 'XXXX';
load([subj '__FilterCluster_SUMU.mat'], 'SUclus');
load(['F344AD_' subj '_TimeStamps.mat'], 'Trial_pulse', 'lastseq')
load(['F344AD_' subj '_SU_Waveform_Output_Extraction'], 'unitIDs', 'duration')
validUnitID = unitIDs(duration>0) + 1;

% Parameters
samplingRate = 30000;
binSize = 0.0005; 
duration = 0.06;
gaussSigma = 0.005;
gaussFactor = 8;
percentPoisson = 99;
rangeStart = 0.0015;
rangeEnd = 0.004;

timePeriods = {
    'baseline', [0 Trial_pulse(1,1)];
    'stim1', [Trial_pulse(1,1) Trial_pulse(end,1)];
    'poststim1', [Trial_pulse(end,1) Trial_pulse(1,2)];
    'stim2', [Trial_pulse(1,2) Trial_pulse(end,2)];
    'poststim2', [Trial_pulse(end,2) Trial_pulse(1,3)];
    'stim3', [Trial_pulse(1,3) Trial_pulse(end,3)];
    'poststim3', [Trial_pulse(end,3) Trial_pulse(1,4)];
    'stim4', [Trial_pulse(1,4) Trial_pulse(end,4)];
    'poststim4', [Trial_pulse(end,4) Trial_pulse(1,5)];
    'stim5', [Trial_pulse(1,5) Trial_pulse(end,5)];
    'poststim5', [Trial_pulse(end,5) Trial_pulse(1,6)];
    'stim6', [Trial_pulse(1,6) Trial_pulse(end,6)];
    'poststim6', [Trial_pulse(end,6) Trial_pulse(1,7)];
    'stim7', [Trial_pulse(1,7) Trial_pulse(end,7)];
    'poststim7', [Trial_pulse(end,7) Trial_pulse(1,8)];
    'stim8', [Trial_pulse(1,8) Trial_pulse(end,8)];
    'poststim8', [Trial_pulse(end,8) Trial_pulse(1,9)];
    'stim9', [Trial_pulse(1,9) Trial_pulse(end,9)];
    'poststim9', [Trial_pulse(end,9) Trial_pulse(1,10)];
    'stim10', [Trial_pulse(1,10) Trial_pulse(end,10)];
    'poststim10', [Trial_pulse(end,10) lastseq(1,1)];
    'stim11', [lastseq(1,1) lastseq(end,1)];
    'poststim', [lastseq(end,1) SUclus.spike_sample(end,1)];
};

ccgResults = struct();
SUclusPeriod = struct(); 

for p = 1:size(timePeriods, 1)
    periodName = timePeriods{p, 1};
    periodRange = timePeriods{p, 2};

    disp(['Processing period: ' periodName])

    periodSpikes = SUclus.spike_sample >= periodRange(1) & SUclus.spike_sample <= periodRange(2);
    
    SUclusPeriod.(periodName).spike_sample = SUclus.spike_sample(periodSpikes);
    SUclusPeriod.(periodName).clusterN = SUclus.clusterN(periodSpikes);

    validCCGUnitID = [];
    
    parfor u1 = 1:size(validUnitID, 2)
        validCCGUnitIDLocal = [];
        for u2 = u1+1:size(validUnitID, 2)
            unit1 = validUnitID(1, u1);
            unit2 = validUnitID(1, u2);
    
            % disp(['Testing Units: ' num2str(unit1) ' and ' num2str(unit2)])
    
            validCCGUnitIDLocal = ccgEval(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod.(periodName), rangeStart, rangeEnd, validCCGUnitIDLocal);
            validCCGUnitIDLocal = ccgEval(unit2, unit1, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod.(periodName), rangeStart, rangeEnd, validCCGUnitIDLocal);

        end
        validCCGUnitID = [validCCGUnitID; validCCGUnitIDLocal];
    end
    
    ccgResults.(periodName) = validCCGUnitID;
end

save([subj '_CCGResults.mat'], 'ccgResults', 'SUclusPeriod')

function [ccg, t] = CCG(spikeTimes1, spikeTimes2, binSize, duration)
    halfBins = round(duration / binSize / 2);
    nBins = 2 * halfBins + 1;
    t = (-halfBins:halfBins)' * binSize;

    ccg = zeros(nBins, 1);

    for i = 1:length(spikeTimes1)
        diffs = spikeTimes2 - spikeTimes1(i);
        binIndices = round(diffs / binSize) + halfBins + 1;
        validIndices = binIndices > 0 & binIndices <= nBins;
        for j = find(validIndices)'
            ccg(binIndices(j)) = ccg(binIndices(j)) + 1;
        end
    end
    ccg = ccg / (length(spikeTimes1) * binSize);
end

function validCCGUnitID = ccgEval(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod, rangeStart, rangeEnd, validCCGUnitID)
    spikeTimes1 = SUclusPeriod.spike_sample(SUclusPeriod.clusterN == unit1) ./ samplingRate;
    spikeTimes2 = SUclusPeriod.spike_sample(SUclusPeriod.clusterN == unit2) ./ samplingRate;
    
    [ccg, t] = CCG(spikeTimes1, spikeTimes2, binSize, duration);

    % Gaussian filtering
    windowSize = round(gaussFactor * gaussSigma / binSize);
    gaussWindow = gausswin(windowSize, gaussSigma / binSize);
    gaussWindow = gaussWindow / sum(gaussWindow);
    baselinePred = conv(ccg, gaussWindow, 'same');

    alpha = 1 - (percentPoisson / 100);

    % Thresholds
    thresholdUpper = poissinv(1 - alpha, baselinePred);
    thresholdLower = poissinv(alpha, baselinePred);

    % Range of Interest
    rangeIndices = find(t >= rangeStart & t <= rangeEnd);
    beforeRangeIndices = find(t < rangeStart);

    % Significant bins within the region of interest
    significantBinsAbove = find(ccg(rangeIndices) > thresholdUpper(rangeIndices));
    significantBinsBelow = find(ccg(rangeIndices) < thresholdLower(rangeIndices));

    % Check for consecutive bins criteria
    isValidComparisonAbove = any(diff(significantBinsAbove) == 1);
    isValidComparisonBelow = any(diff(significantBinsBelow) == 1);

    % Get the first cluster of two or more consecutive bins within the region of interest
    firstConsecutiveAbove = getFirstConsecutiveCluster(significantBinsAbove);
    firstConsecutiveBelow = getFirstConsecutiveCluster(significantBinsBelow);

    consecutiveAboveBefore = getConsecutiveClusters(ccg(beforeRangeIndices) > thresholdUpper(beforeRangeIndices));
    consecutiveBelowBefore = getConsecutiveClusters(ccg(beforeRangeIndices) < thresholdLower(beforeRangeIndices));

    if ~isempty(consecutiveAboveBefore)
        sumAboveBefore = max(cellfun(@(x) sum(ccg(beforeRangeIndices(x))), consecutiveAboveBefore, 'UniformOutput', true));
    else
        sumAboveBefore = 0;
    end
    
    if ~isempty(consecutiveBelowBefore)
        sumBelowBeforeidx = find(cellfun(@(x) numel(x) >= length(significantBinsBelow), consecutiveBelowBefore));
        sumBelowBefore = min(cellfun(@(x) sum(ccg(beforeRangeIndices(x))), consecutiveBelowBefore(sumBelowBeforeidx), 'UniformOutput', true));
    else
        sumBelowBefore = 0;
    end

    sumAboveRegion = sum(ccg(rangeIndices(firstConsecutiveAbove)));
    sumBelowRegion = sum(ccg(rangeIndices(firstConsecutiveAbove)));

    % If only isValidComparisonAbove
    if isValidComparisonAbove && ~isValidComparisonBelow
        if sumAboveRegion > sumAboveBefore
            validCCGUnitID = [validCCGUnitID; unit1, unit2, 1];
            disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' significant excitatory'])
        end
    % If only isValidComparisonBelow
    elseif isValidComparisonBelow && ~isValidComparisonAbove && all(thresholdLower(rangeIndices) > 0) % HELLO
        if sumBelowRegion < sumBelowBefore
            validCCGUnitID = [validCCGUnitID; unit1, unit2, 0];
            disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' significant inhibitory'])
        end
    % If both criteria are satisfied
    elseif isValidComparisonAbove && isValidComparisonBelow && all(thresholdLower(rangeIndices) > 0) % HELLO
        % Compare based on the first cluster of consecutive significant bins
        if ~isempty(firstConsecutiveAbove) && ~isempty(firstConsecutiveBelow) 
            if firstConsecutiveAbove(1) < firstConsecutiveBelow(1)
                % Excitatory first
                if sumAboveRegion > sumAboveBefore
                    validCCGUnitID = [validCCGUnitID; unit1, unit2, 1];
                    disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' significant excitatory'])
                elseif sumBelowRegion < sumBelowBefore
                    validCCGUnitID = [validCCGUnitID; unit1, unit2, 0];
                    disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' significant inhibitory'])
                end
            elseif firstConsecutiveBelow(1) < firstConsecutiveAbove(1)
                % Inhibitory first
                if sumBelowRegion < sumBelowBefore
                    validCCGUnitID = [validCCGUnitID; unit1, unit2, 0];
                    disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' significant inhibitory'])
                elseif sumAboveRegion > sumAboveBefore
                    validCCGUnitID = [validCCGUnitID; unit1, unit2, 1];
                    disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' significant excitatory'])
                end
            end
        end
    end
end

function firstCluster = getFirstConsecutiveCluster(significantBins)
    if length(significantBins) >= 2
        diffs = diff(significantBins);
        idx = find(diffs == 1, 1); 
        if ~isempty(idx)
            firstCluster = significantBins(idx:idx+1);  
        else
            firstCluster = [];
        end
    else
        firstCluster = [];
    end
end

function clusters = getConsecutiveClusters(logicalArray)
    clusters = {};  
    
    if any(logicalArray)
        logicalArray = [false; logicalArray(:); false]; 
        diffArray = diff(logicalArray);
        startIndices = find(diffArray == 1);  
        endIndices = find(diffArray == -1) - 1; 
        
        for i = 1:length(startIndices)
            clusters{end+1} = startIndices(i):endIndices(i);  
        end
    end
end
