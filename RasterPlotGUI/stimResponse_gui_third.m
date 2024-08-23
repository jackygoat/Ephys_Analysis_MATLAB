function create_gui()
    % Create the main GUI window with a grid layout
    fig = uifigure('Name', 'Stimulus Response Analysis', 'Position', [100 100 1400 800]);
    gl = uigridlayout(fig, [3, 2]); 
    gl.RowHeight = {'fit', '1.4x', 'fit'};
    gl.ColumnWidth = {'1x', '1x'};

    % Create a panel for scatter plots
    scatterPanel = uipanel(gl, 'Title', 'Scatter Plots');
    scatterPanel.Layout.Row = 2;
    scatterPanel.Layout.Column = 1;

    % Create a panel for the raster plot
    rasterPanel = uipanel(gl, 'Title', 'Raster Plot');
    rasterPanel.Layout.Row = 2;
    rasterPanel.Layout.Column = 2;

    % Create a panel for input fields
    inputPanel = uipanel(gl, 'Title', 'Parameters');
    inputPanel.Layout.Row = 1;
    inputPanel.Layout.Column = 1;
    inputGrid = uigridlayout(inputPanel, [1, 8]);
    inputGrid.ColumnWidth = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};

    % Create input fields for parameters
    trialLabel = uilabel(inputGrid, 'Text', 'Trial:');
    trialInput = uieditfield(inputGrid, 'numeric', 'Value', 0);

    secBeforeLabel = uilabel(inputGrid, 'Text', 'Sec Before:');
    secBeforeInput = uieditfield(inputGrid, 'numeric', 'Value', 0.1);

    secAfterLabel = uilabel(inputGrid, 'Text', 'Sec After:');
    secAfterInput = uieditfield(inputGrid, 'numeric', 'Value', 0.3);
    
    binWidthLabel = uilabel(inputGrid, 'Text', 'Bin Width (s):');
    binWidthInput = uieditfield(inputGrid, 'numeric', 'Value', 0.004); % Default bin width value

    % Create a panel to display clicked point details
    detailsPanel = uipanel(gl, 'Title', 'Clicked Point Details');
    detailsPanel.Layout.Row = 1;
    detailsPanel.Layout.Column = 2;
    detailsGrid = uigridlayout(detailsPanel, [1, 8]);
    uilabel(detailsGrid, 'Text', 'Field:');
    fieldLabel = uilabel(detailsGrid, 'Text', 'N/A');
    uilabel(detailsGrid, 'Text', 'Unit Index:');
    xValueLabel = uilabel(detailsGrid, 'Text', 'N/A');
    uilabel(detailsGrid, 'Text', 'Unit ID:');
    yValueLabel = uilabel(detailsGrid, 'Text', 'N/A');
    uilabel(detailsGrid, 'Text', '% Stimuli:');
    zValueLabel = uilabel(detailsGrid, 'Text', 'N/A');

    % Create a panel for line plot
    linePanel = uipanel(gl, 'Title', '% Stimulation Over Trials');
    linePanel.Layout.Row = 3; 
    linePanel.Layout.Column = 1;
    
    % Create a panel for histogram
    histogramPanel = uipanel(gl, 'Title', 'Spike Timing Histogram');
    histogramPanel.Layout.Row = 3;
    histogramPanel.Layout.Column = 2;

    % Load data and initialize variables
    mainDir = cd;
    subj = extractAfter(mainDir,'Cohort'); subj = subj(2:5);
    load([subj '_stoch_response_FR.mat'], 'percentageSTIM_excite_cortex', 'percentageSTIM_excite_hipp', 'percentageSTIM_inhib_cortex', 'percentageSTIM_inhib_hipp');
    load(['F344AD_' subj '_TimeStamps.mat'], 'Trial_pulse');
    load([subj '__FilterCluster_SUMU.mat'], 'clusTab');
    load([subj '_layered_PExcitePInhib.mat'], 'excite_cortex_ID', 'excite_hipp_ID', 'inhib_cortex_ID', 'inhib_hipp_ID');

    % Ensure data is squeezed
    percentageSTIM_excite_cortex = squeeze(percentageSTIM_excite_cortex);
    percentageSTIM_excite_hipp = squeeze(percentageSTIM_excite_hipp);
    percentageSTIM_inhib_cortex = squeeze(percentageSTIM_inhib_cortex);
    percentageSTIM_inhib_hipp = squeeze(percentageSTIM_inhib_hipp);

    % Create scatter plots
    % Here we just create the layout, the actual data plotting will happen in the callback
    t = tiledlayout(scatterPanel, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    % Plot % Stim Response Excite Cortex
    ax1 = nexttile(t);
    title(ax1, [subj ' % Stim Response Excite Cortex']);

    % Plot % Stim Response Excite Hipp
    ax2 = nexttile(t);
    title(ax2, [subj ' % Stim Response Excite Hipp']);

    % Plot % Stim Response Inhib Cortex
    ax3 = nexttile(t);
    title(ax3, [subj ' % Stim Response Inhib Cortex']);

    % Plot % Stim Response Inhib Hipp
    ax4 = nexttile(t);
    title(ax4, [subj ' % Stim Response Inhib Hipp']);

    % Add click callback
    scatterPanel.UserData = struct('ax1', ax1, 'ax2', ax2, 'ax3', ax3, 'ax4', ax4, ...
                                   'percentageSTIM_excite_cortex', percentageSTIM_excite_cortex, ...
                                   'percentageSTIM_excite_hipp', percentageSTIM_excite_hipp, ...
                                   'percentageSTIM_inhib_cortex', percentageSTIM_inhib_cortex, ...
                                   'percentageSTIM_inhib_hipp', percentageSTIM_inhib_hipp, ...
                                   'trialInput', trialInput, 'subj', subj, ...
                                   'rasterPanel', rasterPanel, ...
                                   'linePanel', linePanel, ...
                                   'histogramPanel', histogramPanel, ...
                                   'fieldLabel', fieldLabel, ...
                                   'xValueLabel', xValueLabel, ...
                                   'yValueLabel', yValueLabel, ...
                                   'zValueLabel', zValueLabel,...
                                   'Trial_pulse', Trial_pulse, ...
                                   'clusTab', clusTab, ...
                                   'excite_cortex_ID', excite_cortex_ID, ...
                                   'excite_hipp_ID', excite_hipp_ID, ...
                                   'inhib_cortex_ID', inhib_cortex_ID, ...
                                   'inhib_hipp_ID', inhib_hipp_ID, ...
                                   'secBeforeInput', secBeforeInput, ...
                                   'secAfterInput', secAfterInput, ...
                                   'binWidthInput', binWidthInput);

    % Trigger the initial plot
    updateScatterPlots(scatterPanel);

    % Update scatter plots on trialInput change
    trialInput.ValueChangedFcn = @(src, event)updateScatterPlots(scatterPanel);
end

function updateScatterPlots(scatterPanel)
    userData = scatterPanel.UserData;
    trialInput = userData.trialInput.Value;

    if trialInput == 0
        % Use the mean across all trials
        data_excite_cortex = mean(userData.percentageSTIM_excite_cortex, 2);
        data_excite_hipp = mean(userData.percentageSTIM_excite_hipp, 2);
        data_inhib_cortex = mean(userData.percentageSTIM_inhib_cortex, 2);
        data_inhib_hipp = mean(userData.percentageSTIM_inhib_hipp, 2);
    else
        % Use the specific trial
        data_excite_cortex = userData.percentageSTIM_excite_cortex(:, trialInput);
        data_excite_hipp = userData.percentageSTIM_excite_hipp(:, trialInput);
        data_inhib_cortex = userData.percentageSTIM_inhib_cortex(:, trialInput);
        data_inhib_hipp = userData.percentageSTIM_inhib_hipp(:, trialInput);
    end

    scatter(userData.ax1, 1:size(data_excite_cortex, 1), data_excite_cortex, 'filled', 'b');
    scatter(userData.ax2, 1:size(data_excite_hipp, 1), data_excite_hipp, 'filled', 'b');
    scatter(userData.ax3, 1:size(data_inhib_cortex, 1), data_inhib_cortex, 'filled', 'r');
    scatter(userData.ax4, 1:size(data_inhib_hipp, 1), data_inhib_hipp, 'filled', 'r');

    % Update titles to reflect current trial
    trialLabel = 'Mean of All Trials';
    if trialInput ~= 0
        trialLabel = ['Trial: ' num2str(trialInput)];
    end

    title(userData.ax1, [userData.subj ' % Stim Response Excite Cortex ' trialLabel]);
    title(userData.ax2, [userData.subj ' % Stim Response Excite Hipp ' trialLabel]);
    title(userData.ax3, [userData.subj ' % Stim Response Inhib Cortex ' trialLabel]);
    title(userData.ax4, [userData.subj ' % Stim Response Inhib Hipp ' trialLabel]);

    % Axis labels
    xlabel(userData.ax1, 'Unit Index');
    xlabel(userData.ax2, 'Unit Index');
    xlabel(userData.ax3, 'Unit Index');
    xlabel(userData.ax4, 'Unit Index');
    ylabel(userData.ax1, '% Stimulation');
    ylabel(userData.ax2, '% Stimulation');
    ylabel(userData.ax3, '% Stimulation');
    ylabel(userData.ax4, '% Stimulation');

    % Add click callback for each scatter plot
    addClickCallback(userData.ax1, userData, 'excite_cortex');
    addClickCallback(userData.ax2, userData, 'excite_hipp');
    addClickCallback(userData.ax3, userData, 'inhib_cortex');
    addClickCallback(userData.ax4, userData, 'inhib_hipp');
end

function addClickCallback(ax, userData, dataType)
    scatterPlots = findobj(ax, 'Type', 'Scatter');
    for i = 1:length(scatterPlots)
        scatterPlots(i).ButtonDownFcn = @(src, event)scatterClickCallback(src, event, ...
            userData.rasterPanel, userData.linePanel, userData.histogramPanel, userData.subj, dataType, ...
            userData.Trial_pulse, userData.clusTab, userData.excite_cortex_ID, ...
            userData.excite_hipp_ID, userData.inhib_cortex_ID, userData.inhib_hipp_ID, ...
            userData.trialInput, userData.secBeforeInput, userData.secAfterInput, ...
            userData.fieldLabel, userData.xValueLabel, userData.yValueLabel, userData.zValueLabel,...
            userData.percentageSTIM_excite_cortex, userData.percentageSTIM_excite_hipp, ...
            userData.percentageSTIM_inhib_cortex, userData.percentageSTIM_inhib_hipp, ...
            userData.binWidthInput);
    end
end

function scatterClickCallback(src, event, rasterPanel, linePanel, histogramPanel, subj, dataType, Trial_pulse, clusTab, excite_cortex_ID, excite_hipp_ID, inhib_cortex_ID, inhib_hipp_ID, trialInput, secBeforeInput, secAfterInput, fieldLabel, xValueLabel, yValueLabel, zValueLabel, percentageSTIM_excite_cortex, percentageSTIM_excite_hipp, percentageSTIM_inhib_cortex, percentageSTIM_inhib_hipp, binWidthInput)
    % Get the clicked point's X value
    clickedPoint = round(event.IntersectionPoint(1));
    yValue = event.IntersectionPoint(2);

    % Update display fields
    fieldLabel.Text = dataType;
    xValueLabel.Text = num2str(clickedPoint);
    zValueLabel.Text = num2str(yValue);

    % Get input values for parameters
    trial = trialInput.Value;
    secBefore = secBeforeInput.Value;
    secAfter = secAfterInput.Value;
    binWidth = binWidthInput.Value;

    % Select the appropriate unit IDs
    switch dataType
        case 'excite_cortex'
            unitIDs = excite_cortex_ID(:, 1);
            color = 'c';
            percentageSTIM = percentageSTIM_excite_cortex;
        case 'excite_hipp'
            unitIDs = excite_hipp_ID(:, 1);
            color = 'c';
            percentageSTIM = percentageSTIM_excite_hipp;
        case 'inhib_cortex'
            unitIDs = inhib_cortex_ID(:, 1);
            color = 'm';
            percentageSTIM = percentageSTIM_inhib_cortex;
        case 'inhib_hipp'
            unitIDs = inhib_hipp_ID(:, 1);
            color = 'm';
            percentageSTIM = percentageSTIM_inhib_hipp;
        otherwise
            error('Unknown data type');
    end

    % Get the unit ID
    unitID = unitIDs(clickedPoint);
    yValueLabel.Text = num2str(unitID);

    % Define parameters
    samplingRate = 30000;

    % Filter cluster ID
    filtered_clusTab = clusTab(clusTab.clusterN == unitID, :);

    % Initialize variables
    clusTabVar = [];
    stimNum = 0;

    if trial == 0
        % Define offset for each trial to separate them on the y-axis
        trialOffset = 100; % Adjust this value as needed for separation

        % Loop through all trials and pulses
        for tr = 1:size(Trial_pulse, 2)
            for i = 1:size(Trial_pulse, 1)
                % Define timeframe
                time_point = Trial_pulse(i, tr);
                start_time = time_point - secBefore * samplingRate;
                end_time = time_point + secAfter * samplingRate;

                % Filter data
                stimNum = stimNum + 1;
                temp = filtered_clusTab(filtered_clusTab.spike_sample >= start_time & filtered_clusTab.spike_sample <= end_time, :);
                normTimestamp = temp.spike_sample - time_point;
                normTimestampSec = normTimestamp / samplingRate;
                clusTabVar = [clusTabVar; temp.spike_sample, normTimestamp, normTimestampSec, temp.clusterN, repmat(stimNum + (tr - 1) * trialOffset, size(temp.clusterN, 1), 1)];
            end
        end
    else
        for i = 1:size(Trial_pulse, 1)
            % Define timeframe
            time_point = Trial_pulse(i, trial);
            start_time = time_point - secBefore * samplingRate;
            end_time = time_point + secAfter * samplingRate;

            % Filter data
            stimNum = stimNum + 1;
            temp = filtered_clusTab(filtered_clusTab.spike_sample >= start_time & filtered_clusTab.spike_sample <= end_time, :);
            normTimestamp = temp.spike_sample - time_point;
            normTimestampSec = normTimestamp / samplingRate;
            clusTabVar = [clusTabVar; temp.spike_sample, normTimestamp, normTimestampSec, temp.clusterN, repmat(stimNum, size(temp.clusterN, 1), 1)];
        end
    end

    % Clear the previous raster plot
    delete(rasterPanel.Children);

    % Create a grid layout to fill the rasterPanel
    rasterGrid = uigridlayout(rasterPanel, [1, 1]);

    % Create raster plot
    ax = uiaxes(rasterGrid);
    scatter(ax, clusTabVar(:, 3), clusTabVar(:, 5), 4, 'filled', color);
    xline(ax, 0, 'r', 'LineWidth', 1);
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Stimulus #');
    title(ax, 'Raster Plot');

    % Adjust y-axis to separate trials when plotting all trials
    if trial == 0
        yLimits = [0, max(clusTabVar(:, 5))];
        if numel(yLimits) == 2 && yLimits(1) < yLimits(2)
            ax.YLim = yLimits;
        end
    end

    % Make the plot fill the entire rasterPanel
    ax.Layout.Row = 1;
    ax.Layout.Column = 1;

    % Plot % Stimulation over trials for the specific data point
    plotStimulationOverTrials(linePanel, percentageSTIM, clickedPoint, color);

    % Plot histogram of spike timing
    plotHistogram(histogramPanel, clusTabVar(:, 3), binWidth, secBefore, secAfter, color);
end

function plotStimulationOverTrials(linePanel, percentageSTIM, clickedPoint, color)
    % Clear the previous line plot
    delete(linePanel.Children);

    % Create a grid layout to fill the linePanel
    lineGrid = uigridlayout(linePanel, [1, 1]);

    % Create line plot
    ax = uiaxes(lineGrid);
    plot(ax, 1:size(percentageSTIM, 2), percentageSTIM(clickedPoint, :), 'Color', color, 'LineWidth', 2);
    xlabel(ax, 'Trial');
    ylabel(ax, '% Stimulation');
    ylim(ax, [0 100]);
    title(ax, '% Stimulation Over Trials');
    ax.Layout.Row = 1;
    ax.Layout.Column = 1;
end

function plotHistogram(histogramPanel, data, binWidth, secBefore, secAfter, color)
    % Clear the previous histogram plot
    delete(histogramPanel.Children);

    % Calculate the number of bins
    nbins = ceil((max(data) - min(data)) / binWidth);

    % Create a grid layout to fill the histogramPanel
    histogramGrid = uigridlayout(histogramPanel, [1, 1]);

    % Create histogram plot
    ax = uiaxes(histogramGrid);
    hold(ax, 'on');
    h = histfit(ax, data, nbins, 'kernel');
    h(1).FaceColor = color;
    h(1).FaceAlpha = 0.5;
    h(1).EdgeColor = 'none';
    h(2).Color = color;

    % Find the max point of the fitted line
    [~, idx] = max(h(2).YData);
    maxPoint = [h(2).XData(idx), h(2).YData(idx)];

    % Plot the max point
    plot(ax, maxPoint(1), maxPoint(2), 'o', 'MarkerEdgeColor', color, 'MarkerFaceColor', color);

    hold(ax, 'off');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Count');
    title(ax, 'Spike Timing Histogram');
    xlim(ax, [-secBefore secAfter]);
    ax.Layout.Row = 1;
    ax.Layout.Column = 1;
end
