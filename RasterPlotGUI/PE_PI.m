clear;
animal_ids = [XXXX];

threshold = XXXX;

for i = 1:length(animal_ids)
    filename = sprintf('F344AD_%d_SU_Waveform_Output_Extraction.mat', animal_ids(i)); % SUMU
    
    data = load(filename);
    
    % variables
    duration = data.duration;
    idx = data.idx;
    meanwaveForms = data.meanwaveForms;
    peaktroughratio = data.peaktroughratio;
    unitIDs = data.unitIDs;

    % Remove any duration that is less than 0
    duration = duration / 30;
    valid_idx = duration >= 0;
    duration = duration(valid_idx);
    idx = idx(valid_idx);
    meanwaveForms = meanwaveForms(valid_idx, :);
    peaktroughratio = peaktroughratio(valid_idx);
    unitIDs = unitIDs(valid_idx);

    % Split 
    over_threshold_idx = duration > threshold;
    below_threshold_idx = duration <= threshold;
    
    % over
    duration_over = duration(over_threshold_idx);
    idx_over = idx(over_threshold_idx);
    meanwaveForms_over = meanwaveForms(over_threshold_idx, :);
    peaktroughratio_over = peaktroughratio(over_threshold_idx);
    unitIDs_over = unitIDs(over_threshold_idx);
    
    % below 
    duration_below = duration(below_threshold_idx);
    idx_below = idx(below_threshold_idx);
    meanwaveForms_below = meanwaveForms(below_threshold_idx, :);
    peaktroughratio_below = peaktroughratio(below_threshold_idx);
    unitIDs_below = unitIDs(below_threshold_idx);
    
  
    save(sprintf('F344AD_SU_%d_Over_Threshold.mat', animal_ids(i)), 'duration_over', 'idx_over', 'meanwaveForms_over', 'peaktroughratio_over', 'unitIDs_over');
    save(sprintf('F344AD_SU_%d_Below_Threshold.mat', animal_ids(i)), 'duration_below', 'idx_below', 'meanwaveForms_below', 'peaktroughratio_below', 'unitIDs_below');
end
