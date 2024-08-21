mainDir = 'INSERT_DIRECTORY'

cd(mainDir)
subj = '####';
load([subj '__FilterCluster_SUMU.mat'], 'SUclus');
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

% init
validCCGUnitID = [];

parfor u1 = 1:size(validUnitID,2)
    localValidCCGUnitID = [];  % Local variable to avoid conflict in parfor
    for u2 = u1+1:size(validUnitID,2)
        unit1 = validUnitID(1,u1);
        unit2 = validUnitID(1,u2);

        disp(['Testing Units: ' num2str(unit1) ' and ' num2str(unit2)])

        localValidCCGUnitID = ccgEval(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclus, rangeStart, rangeEnd, localValidCCGUnitID);
        localValidCCGUnitID = ccgEval(unit2, unit1, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclus, rangeStart, rangeEnd, localValidCCGUnitID);
    end
    % Append local results to the global variable
    validCCGUnitID = [validCCGUnitID; localValidCCGUnitID];
end

save([subj '_validCCGUnitID.mat'], 'validCCGUnitID')

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

    % gaussian
    windowSize = round(gaussFactor * gaussSigma / binSize);
    gaussWindow = gausswin(windowSize, gaussSigma / binSize);
    gaussWindow = gaussWindow / sum(gaussWindow);
    baselinePred = conv(ccg, gaussWindow, 'same');

    % poisson
    alpha = 1 - (percentPoisson./100);

    % thresholds
    thresholdUpper = poissinv(1 - alpha, baselinePred);
    thresholdLower = poissinv(alpha, baselinePred);

    % significance
    rangeIndices = find(t >= rangeStart & t <= rangeEnd);
    significantBinsAbove = find(ccg(rangeIndices) > thresholdUpper(rangeIndices));
    significantBinsBelow = find(ccg(rangeIndices) < thresholdLower(rangeIndices));

    % check for consecutive bins outside the range that are 0
    beforeRangeIndices = find(t < rangeStart);
    afterRangeIndices = find(t > rangeEnd);
    zeroBinsBefore = find(ccg(beforeRangeIndices) == 0);
    zeroBinsAfter = find(ccg(afterRangeIndices) == 0);

    hasConsecutiveZerosBefore = any(diff(zeroBinsBefore) == 1);
    hasConsecutiveZerosAfter = any(diff(zeroBinsAfter) == 1);

    % if significant and thresholdLower is greater than 0 and no consecutive zeros before/after the range
    isValidComparisonAbove = any(diff(significantBinsAbove) == 1);
    isValidComparisonBelow = any(diff(significantBinsBelow) == 1) && ...
                             all(thresholdLower(rangeIndices) > 0) && ...
                             ~(hasConsecutiveZerosBefore || hasConsecutiveZerosAfter);

    if isValidComparisonAbove && all(thresholdLower(rangeIndices) > 0)
        validCCGUnitID = [validCCGUnitID; unit1, unit2, 1];
        disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' above threshold ----------------'])
    end
    if isValidComparisonBelow && all(thresholdLower(rangeIndices) > 0)
        validCCGUnitID = [validCCGUnitID; unit1, unit2, 0];
        disp(['Unit: ' num2str(unit1) ' -> ' num2str(unit2) ' below baseline ----------------'])
    end
end