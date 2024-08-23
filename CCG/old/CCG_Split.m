clear;

mainDir = uigetdir();
cd(mainDir)
subj = '4760_HPC';
load([subj '__FilterCluster_SUMU.mat'], 'SUclus');
load(['F344AD_' subj '_TimeStamps.mat'], 'Trial_pulse', 'lastseq')
load(['F344AD_' subj '_SU_Waveform_Output_Extraction'], 'unitIDs', 'duration')
validUnitID = unitIDs(duration>0) + 1;

% parameter
samplingRate = 30000;
binSize = 0.0005; 
duration = 0.06;
gaussSigma = 0.005;
gaussFactor = 8;
percentPoisson = 99;
rangeStart = 0.0015;
rangeEnd = 0.004;

% time periods
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

validCCGUnitID = [];
ccgResults = struct();

for p = 1:size(timePeriods, 1)
    periodName = timePeriods{p, 1};
    periodRange = timePeriods{p, 2};

    disp(['Processing period: ' periodName])

    periodSpikes = SUclus.spike_sample >= periodRange(1) & SUclus.spike_sample <= periodRange(2);
    SUclusPeriod = struct('spike_sample', SUclus.spike_sample(periodSpikes), 'clusterN', SUclus.clusterN(periodSpikes));

    for u1 = 1:size(validUnitID, 2)
        for u2 = u1+1:size(validUnitID, 2)
            unit1 = validUnitID(1, u1);
            unit2 = validUnitID(1, u2);

            disp(['Testing Units: ' num2str(unit1) ' and ' num2str(unit2)])

            validCCGUnitID = ccgEval(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod, rangeStart, rangeEnd, validCCGUnitID);
            validCCGUnitID = ccgEval(unit2, unit1, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod, rangeStart, rangeEnd, validCCGUnitID);
        end
    end

    ccgResults.(periodName) = validCCGUnitID;
    validCCGUnitID = []; 
end

save([subj '_CCGResults.mat'], 'ccgResults')

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

function validCCGUnitID = ccgEval(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclus, rangeStart, rangeEnd, validCCGUnitID)
    spikeTimes1 = SUclus.spike_sample(SUclus.clusterN == unit1) ./ samplingRate;
    spikeTimes2 = SUclus.spike_sample(SUclus.clusterN == unit2) ./ samplingRate;
    
    [ccg, t] = CCG(spikeTimes1, spikeTimes2, binSize, duration);

    windowSize = round(gaussFactor * gaussSigma / binSize);
    gaussWindow = gausswin(windowSize, gaussSigma / binSize);
    gaussWindow = gaussWindow / sum(gaussWindow);
    baselinePred = conv(ccg, gaussWindow, 'same');

    alpha = 1 - (percentPoisson / 100);

    thresholdUpper = poissinv(1 - alpha, baselinePred);
    thresholdLower = poissinv(alpha, baselinePred);

    rangeIndices = find(t >= rangeStart & t <= rangeEnd);
    significantBinsAbove = find(ccg(rangeIndices) > thresholdUpper(rangeIndices));
    significantBinsBelow = find(ccg(rangeIndices) < thresholdLower(rangeIndices));

    beforeRangeIndices = find(t < rangeStart);
    afterRangeIndices = find(t > rangeEnd);
    zeroBinsBefore = find(ccg(beforeRangeIndices) == 0);
    zeroBinsAfter = find(ccg(afterRangeIndices) == 0);

    hasConsecutiveZerosBefore = any(diff(zeroBinsBefore) == 1);
    hasConsecutiveZerosAfter = any(diff(zeroBinsAfter) == 1);

    isValidComparisonAbove = any(diff(significantBinsAbove) == 1);
    isValidComparisonBelow = any(diff(significantBinsBelow) == 1) && all(thresholdLower(rangeIndices) > 0) && ~(hasConsecutiveZerosBefore || hasConsecutiveZerosAfter);
    
    if isValidComparisonAbove && all(thresholdLower(rangeIndices) > 0)
        validCCGUnitID = [validCCGUnitID; unit1, unit2, 1];
        disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' above threshold ----------------'])
    end
    if isValidComparisonBelow && all(thresholdLower(rangeIndices) > 0)
        validCCGUnitID = [validCCGUnitID; unit1, unit2, 0];   
        disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' below baseline ----------------'])
    end
end
