%% data extraction _ stochastic responsiveness - Hipp only

clear;
mainDir = uigetdir();
cd(mainDir)

subj = 'XXXX';
% subj = extractAfter(mainDir,'Cohort'); subj = subj(2:5);
% subj = extractAfter(mainDir, 'AD_'); subj = subj(1:4);
% subj = extractAfter(mainDir,'Emma\'); subj = subj(10:13);

load([subj '_HPC__FilterCluster_SUMU.mat'], 'clusY')
load([subj '_HPC_stochastic responsiveness .mat'], 'SUResTrial', 'MUResTrial')
load(['F344AD_SU_' subj '_HPC_Below_Threshold.mat'], 'unitIDs_below')
load(['F344AD_SU_' subj '_HPC_Over_Threshold.mat'], 'unitIDs_over')

hippThresh_lower = XXXX;
hippThresh_upper = XXXX;

inhib_hipp = [];
inhib_hipp_ID = [];
excite_hipp = [];
excite_hipp_ID = [];

% SU
for i = 1:size(SUResTrial, 1)
    unitID = clusY(i, 1);
    depth = clusY(i, 2);
    
    SUTrialFR = SUResTrial(i, :, :);
    
    if ismember(unitID, unitIDs_below+1)
        % Inhibitory 
        if depth > hippThresh_lower && depth < hippThresh_upper
            inhib_hipp = [inhib_hipp; SUTrialFR];
            inhib_hipp_ID = [inhib_hipp_ID; unitID, depth];
        else
            continue
        end
    elseif ismember(unitID, unitIDs_over+1)
        % Excitatory 
        if depth > hippThresh_lower && depth < hippThresh_upper
            excite_hipp = [excite_hipp; SUTrialFR];
            excite_hipp_ID = [excite_hipp_ID; unitID, depth];
        else
            continue
        end
    end
end

[x_hipp, y_stim, z_trials] = size(excite_hipp);
[x_inhib_hipp, ~, ~] = size(inhib_hipp);

% calculate percentage of responding neurons
calculate_percentage = @(data, x_dim, y_dim) sum(data, 1) / x_dim * 100;

percentage_excite_hipp = calculate_percentage(excite_hipp, x_hipp, y_stim);
percentage_inhib_hipp = calculate_percentage(inhib_hipp, x_inhib_hipp, y_stim);

% calculate percentage of responding stimulus
calculate_percentage_clusters = @(data, y_dim) sum(data, 2) / y_dim * 100;
percentageSTIM_excite_hipp = calculate_percentage_clusters(excite_hipp, y_stim);
percentageSTIM_inhib_hipp = calculate_percentage_clusters(inhib_hipp, y_stim);

percentage_hipp_MU = [];
percentage_inhib_cortex = [];
percentage_excite_cortex = [];
percentage_cortex_MU = [];
percentageSTIM_excite_cortex = [];
percentageSTIM_inhib_cortex = [];
inhib_cortex_ID = [];
excite_cortex_ID = [];

% cd G:\getWave\Waveform_Extraction
save ([subj '_HPC_stoch_response_FR.mat'], 'percentage_inhib_hipp', 'percentage_inhib_cortex', 'percentage_excite_hipp', 'percentage_excite_cortex', 'percentage_cortex_MU', 'percentage_hipp_MU', 'percentageSTIM_excite_cortex', 'percentageSTIM_excite_hipp', 'percentageSTIM_inhib_cortex', 'percentageSTIM_inhib_hipp')
save ([subj '_HPC_layered_PExcitePInhib.mat'], 'inhib_cortex_ID', 'inhib_hipp_ID', 'excite_cortex_ID', 'excite_hipp_ID')
cd(mainDir)
clear
