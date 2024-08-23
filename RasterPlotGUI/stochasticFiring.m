clear;
mainDir = uigetdir();
cd(mainDir)

subj = 'XXXX';
% subj = extractAfter(mainDir,'Cohort'); subj = subj(2:5);
% subj = extractAfter(mainDir, 'AD_'); subj = subj(1:4);
% subj = extractAftsqer(mainDir,'Emma\'); subj = subj(10:13);

load([subj '__FilterCluster_SUMU.mat'], 'clusY')
load([subj '_stochastic responsiveness .mat'], 'SUResTrial', 'MUResTrial')
load(['F344AD_SU_' subj '_Below_Threshold.mat'], 'unitIDs_below')
load(['F344AD_SU_' subj '_Over_Threshold.mat'], 'unitIDs_over')

cortexThresh = XXXX; 
hippThreshLower = XXXX; 
hippThreshUpper = XXXX; 

inhib_cortex = [];
inhib_cortex_ID = [];
excite_cortex = [];
excite_cortex_ID = [];
inhib_hipp = [];
inhib_hipp_ID = [];
excite_hipp = [];
excite_hipp_ID = [];

% SU
for i = 1:size(SUResTrial, 1)
    unitID = clusY(i, 1);
    depth = clusY(i, 2);
    
    SUTrialFR = SUResTrial(i, :, :);
    
    % inhibitory or excitatory
    if ismember(unitID, unitIDs_below+1)
        % Inhibitory 
        if depth < cortexThresh
            inhib_cortex = [inhib_cortex; SUTrialFR];
            inhib_cortex_ID = [inhib_cortex_ID; unitID, depth];
        elseif depth >= hippThreshLower && depth <= hippThreshUpper
            inhib_hipp = [inhib_hipp; SUTrialFR];
            inhib_hipp_ID = [inhib_hipp_ID; unitID, depth];
        else
            continue
        end
    elseif ismember(unitID, unitIDs_over+1)
        % Excitatory 
        if depth < cortexThresh
            excite_cortex = [excite_cortex; SUTrialFR];
            excite_cortex_ID = [excite_cortex_ID; unitID, depth];
        elseif depth >= hippThreshLower && depth <= hippThreshUpper
            excite_hipp = [excite_hipp; SUTrialFR];
            excite_hipp_ID = [excite_hipp_ID; unitID, depth];
        else
            continue
        end
    end
end

[x_cortex, y_stim, z_trials] = size(excite_cortex);
[x_hipp, ~, ~] = size(excite_hipp);
[x_inhib_cortex, ~, ~] = size(inhib_cortex);
[x_inhib_hipp, ~, ~] = size(inhib_hipp);

% calculate percentage of responding neurons
calculate_percentage = @(data, x_dim, y_dim) sum(data, 1) / x_dim * 100;

percentage_excite_cortex = calculate_percentage(excite_cortex, x_cortex, y_stim);
percentage_excite_hipp = calculate_percentage(excite_hipp, x_hipp, y_stim);
percentage_inhib_cortex = calculate_percentage(inhib_cortex, x_inhib_cortex, y_stim);
percentage_inhib_hipp = calculate_percentage(inhib_hipp, x_inhib_hipp, y_stim);

% calculate percentage of responding stimulus
calculate_percentage_clusters = @(data, y_dim) sum(data, 2) / y_dim * 100;
percentageSTIM_excite_cortex = calculate_percentage_clusters(excite_cortex, y_stim);
percentageSTIM_excite_hipp = calculate_percentage_clusters(excite_hipp, y_stim);
percentageSTIM_inhib_cortex = calculate_percentage_clusters(inhib_cortex, y_stim);
percentageSTIM_inhib_hipp = calculate_percentage_clusters(inhib_hipp, y_stim);

percentage_cortex_MU = [];
percentage_hipp_MU = [];

% MU
% load(['F344AD_MU_' subj '_Below_Threshold.mat'], 'unitIDs_below')
% load(['F344AD_MU_' subj '_Over_Threshold.mat'], 'unitIDs_over')
% unitIDs_MU = [unitIDs_below, unitIDs_over];
% 
% cortex = [];
% cortex_ID = [];
% hipp = [];
% hipp_ID = [];
% 
% for i = 1:size(MUResTrial, 1)
%     unitID = clusY(i, 1);
%     depth = clusY(i, 2);
%     
%     MUTrialFR = MUResTrial(i, :, :);
%     if ismember(unitID, unitIDs_MU+1)
%         if depth < cortexThresh
%             cortex = [cortex; MUTrialFR];
%             cortex_ID = [cortex_ID; unitID, depth];
%         elseif depth >= hippThreshLower && depth <= hippThreshUpper
%             hipp = [hipp; MUTrialFR];
%             hipp_ID = [hipp_ID; unitID, depth];
%         else
%             continue
%         end
%     end
% end
% 
% [x_cortex, y_stim, z_trials] = size(cortex);
% [x_hipp, ~, ~] = size(hipp);
% 
% percentage_cortex_MU = calculate_percentage(cortex, x_cortex, y_stim);
% percentage_hipp_MU = calculate_percentage(hipp, x_hipp, y_stim);

% cd G:\getWave\Waveform_Extraction
save ([subj '_stoch_response_FR.mat'], 'percentage_inhib_hipp', 'percentage_inhib_cortex', 'percentage_excite_hipp', 'percentage_excite_cortex', 'percentage_cortex_MU', 'percentage_hipp_MU', 'percentageSTIM_excite_cortex', 'percentageSTIM_excite_hipp', 'percentageSTIM_inhib_cortex', 'percentageSTIM_inhib_hipp')
save ([subj '_layered_PExcitePInhib.mat'], 'inhib_cortex_ID', 'inhib_hipp_ID', 'excite_cortex_ID', 'excite_hipp_ID')
cd(mainDir)
clear
