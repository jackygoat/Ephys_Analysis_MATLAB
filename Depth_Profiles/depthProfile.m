clear;

mainDir = uigetdir();
cd(mainDir)

subj = 'XXXX';
% subj = extractAfter(mainDir,'Cohort'); subj = subj(2:5);
% subj = extractAfter(mainDir,'AD_'); subj = subj(1:4);
% subj = extractAfter(mainDir,'Emma\'); subj = subj(10:13);
load ([subj '__FilterCluster_SUMU.mat'])
load (['F344AD_' subj '_SU_Waveform_Output_Extraction'])

SUduration = duration;
SUpeaktotrough = peaktroughratio;
SUdepth = SUdis(:, 2);

% clear duration peaktroughratio;
% load (['F344AD_' subj '_SU_Waveform_Output_Extraction'])
% MUduration = duration;
% MUpeaktotrough = peaktroughratio;
% MUdepth = MUdis(:, 2);

combinedDuration = SUduration;
combinedPeaktrough = SUpeaktotrough;
combinedDepth = SUdepth';

% ----------peak to trough 
figure;
scatter(combinedPeaktrough, combinedDepth, 'filled');
title([subj '-Cluster- Peak to Trough Ratio vs. Depth']);
xlabel('Peak to Trough Ratio');
ylabel('Depth');
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'YDir','reverse');
grid on;
set(gcf, 'Color', 'w', 'Theme', 'light');
hold on

% Fitting peak to trough
binSize = 200;
stepSize =10;
depthRange = min(combinedDepth):stepSize:max(combinedDepth);

unitSums  = zeros(1, length(depthRange));

for i = 1:length(depthRange)
    binStart = depthRange(i) - binSize/2;
    binEnd = depthRange(i) + binSize/2;
    unitSums(i) = median(combinedPeaktrough(combinedDepth >= binStart & combinedDepth < binEnd));
end

% Smoothing peak to trough
w = gausswin(13,6);
w = w/sum(w);
smoothedSums = filter(w, 1, unitSums);

% Plot line peak to trough
resizeFactor = 1.5;
plot(smoothedSums*resizeFactor, depthRange, 'r', 'LineWidth', 2); 

saveas (gcf, [subj '-Cluster- Peak to Trough Ratio vs. Depth.fig'])

% ---------Duration threshhold
rmvdindx = find(combinedDuration > -20);
combinedDuration = combinedDuration(rmvdindx);
combinedDepth = combinedDepth(rmvdindx);

% duration
figure;
scatter(combinedDuration, combinedDepth, 'filled');
title([subj '-Cluster- Duration vs. Depth']);
xlabel('Duration');
ylabel('Depth');
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'YDir','reverse');
grid on;
set(gcf, 'Color', 'w', 'Theme', 'light');
hold on

% Fitting duration
binSize = 200;
stepSize =20;
depthRange = min(combinedDepth):stepSize:max(combinedDepth);

unitSums  = zeros(1, length(depthRange));
    
SigncombinedDuration = sign(combinedDuration);

for i = 1:length(depthRange)
    binStart = depthRange(i) - binSize/2;
    binEnd = depthRange(i) + binSize/2;
    unitSums(i) = sum(SigncombinedDuration(combinedDepth >= binStart & combinedDepth < binEnd));
end

% Smoothing duration
w = gausswin(11,6);
w = w/sum(w);
smoothedSums = filter(w, 1, unitSums);

% Plot line duration
resizeFactor = 1;
plot(smoothedSums*resizeFactor, depthRange, 'b', 'LineWidth', 2); 

saveas (gcf, [subj '-Cluster- Duration vs. Depth.fig'])

% --------- Duration vs Depth (below 0 in red, above 0 in blue)
figure;
scatter(combinedDuration(combinedDuration < 0), combinedDepth(combinedDuration < 0), 'r', 'filled');
hold on;
scatter(combinedDuration(combinedDuration >= 0), combinedDepth(combinedDuration >= 0), 'b', 'filled');
title([subj '-Cluster- Duration vs. Depth (Separated)']);
xlabel('Duration');
ylabel('Depth');
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'YDir','reverse');
grid on;
set(gcf, 'Color', 'w', 'Theme', 'light');

% Fitting duration for below 0
binSize = 500;
stepSize =20;
unitSumsBelow  = zeros(1, length(depthRange));
for i = 1:length(depthRange)
    binStart = depthRange(i) - binSize/2;
    binEnd = depthRange(i) + binSize/2;
    unitSumsBelow(i) = mean(combinedDuration(combinedDepth >= binStart & combinedDepth < binEnd & combinedDuration < 0));
end

% Smoothing duration below 0
smoothedSumsBelow = filter(w, 1, unitSumsBelow);

% Plot line duration below 0
% plot(smoothedSumsBelow*resizeFactor, depthRange, 'r', 'LineWidth', 2); 

% Fitting duration for above 0
unitSumsAbove  = zeros(1, length(depthRange));
for i = 1:length(depthRange)
    binStart = depthRange(i) - binSize/2;
    binEnd = depthRange(i) + binSize/2;
    unitSumsAbove(i) = mean(combinedDuration(combinedDepth >= binStart & combinedDepth < binEnd & combinedDuration >= 0));
end

% Smoothing duration above 0
smoothedSumsAbove = filter(w, 1, unitSumsAbove);

% Plot line duration above 0
% plot(smoothedSumsAbove*resizeFactor, depthRange, 'b', 'LineWidth', 2); 
hold off;

saveas (gcf, [subj '-Cluster- Duration vs. Depth (Separated).fig'])

% --------- Histograms of Duration vs Depth ---------
figure;
hold on;
binCount = 50;
[counts1, centers1] = hist(combinedDepth(combinedDuration < 0), binCount);
barh(centers1, -counts1, 'FaceColor', 'r');
[counts2, centers2] = hist(combinedDepth(combinedDuration >= 0), binCount);
barh(centers2, counts2, 'FaceColor', 'b');
title([subj '-Cluster- Duration vs. Depth Histogram']);
xlabel('Count');
ylabel('Depth');
set(gca, 'FontSize', 12, 'FontName', 'Arial', 'YDir','reverse');
grid on;
set(gcf, 'Color', 'w', 'Theme', 'light');
hold off;

saveas (gcf, [subj '-Cluster- Duration vs. Depth Histogram.fig'])
