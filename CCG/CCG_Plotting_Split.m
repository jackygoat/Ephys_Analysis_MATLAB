%% CCG Plotting

mainDir = uigetdir();
cd(mainDir)
subj = 'XXXX';

% Load CCG Results and other necessary data
load([subj '_CCGResults.mat'], 'ccgResults', 'SUclusPeriod');
load([subj '__FilterCluster_SUMU.mat'], 'SUclus', 'clusY');
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

% Iterate through the different time periods and plot CCG
timePeriods = fieldnames(ccgResults);

periodName = ['XXXX'];
field = ccgResults.(periodName);


for k = 1:size(field, 1)
    unit1 = field(k, 1);
    unit2 = field(k, 2);
    ccgPlotting(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod.(periodName), rangeStart, rangeEnd, clusY, periodName);
    pause;
end

% ccgPlotting(284, 181, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod.(periodName), rangeStart, rangeEnd, clusY, periodName);


function ccgPlotting(unit1, unit2, samplingRate, binSize, duration, gaussFactor, gaussSigma, percentPoisson, SUclusPeriod, rangeStart, rangeEnd, clusY, periodName)
    spikeTimes1 = SUclusPeriod.spike_sample(SUclusPeriod.clusterN == unit1) ./ samplingRate;
    spikeTimes2 = SUclusPeriod.spike_sample(SUclusPeriod.clusterN == unit2) ./ samplingRate;
    
    [ccg, t] = CCG(spikeTimes1, spikeTimes2, binSize, duration);

    % Gaussian 
    windowSize = round(gaussFactor * gaussSigma / binSize);
    gaussWindow = gausswin(windowSize, gaussSigma / binSize);
    gaussWindow = gaussWindow / sum(gaussWindow);
    baselinePred = conv(ccg, gaussWindow, 'same');
    
    % Poisson thresholds
    alpha = 1 - (percentPoisson / 100);
    thresholdUpper = poissinv(1 - alpha, baselinePred);
    thresholdLower = poissinv(alpha, baselinePred);

    depth1 = clusY(clusY(:,1) == unit1, 2);
    depth2 = clusY(clusY(:,1) == unit2, 2);

    % Plot the CCG
    % figure;
    bar(t, ccg, 'hist');
    hold on;
    xline(rangeStart, 'k--');
    xline(rangeEnd, 'k--');
    plot(t, thresholdLower, 'r--', 'LineWidth', 1);
    plot(t, thresholdUpper, 'r--', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Counts');
    xlim([-duration/2 duration/2]);
    title(['CCG ' periodName ': Unit ' num2str(unit1) ' (' num2str(depth1) 'um) -> Unit ' num2str(unit2) ' (' num2str(depth2) 'um)']);
    hold off;
    set(gcf, "Theme", "Light")
end

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

%%

stimNames = {'stim1', 'stim2', 'stim3', 'stim4', 'stim5', 'stim6', 'stim7', 'stim8','stim9', 'stim10', 'stim11'};
commonResults = struct(); 

for i = 1:length(stimNames)-1
    for j = i+1:length(stimNames)
        stimA = stimNames{i};
        stimB = stimNames{j};
        commonRows = intersect(ccgResults.(stimA)(:,1:2), ccgResults.(stimB)(:,1:2), 'rows');
        fieldName = ['common_', stimA, '_', stimB];
        compareNames = ['percent_' stimA '_' stimB];
        commonResults.(fieldName) = commonRows;
        commonResults.(compareNames) = [size(commonRows,1) ./ size(ccgResults.(stimA),1), size(commonRows,1) ./ size(ccgResults.(stimB),1)];
    end
end

combinedMatrix = cell(length(stimNames));

for i = 1:length(stimNames)
    for j = 1:length(stimNames)
        if isnan(percentageMatrixA(i, j)) && isnan(percentageMatrixB(i, j))
            combinedMatrix{i, j} = ''; 
        else
            combinedMatrix{i, j} = sprintf('%.3f / %.3f', percentageMatrixA(i, j), percentageMatrixB(i, j));
        end
    end
end
combinedTable = cell2table(combinedMatrix, 'VariableNames', stimNames, 'RowNames', stimNames)