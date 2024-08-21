addpath 'INSERT_DIRECTORY'
load ('INSERT_DIRECTORY/4759__FilterCluster_SUMU.mat')
dataDir = 'INSERT_DIRECTORY';
kilosubj = '####';

% SU
sp = loadKSdir(dataDir);
gwfparams.sr = sp.sample_rate;
gwfparams.dataDir = dataDir;
gwfparams.fileName = 'temp_wh.dat';
gwfparams.dataType = sp.dtype;
gwfparams.nCh = sp.n_channels_dat;
gwfparams.wfWin = [-(0.001*gwfparams.sr) 0.002*gwfparams.sr];
gwfparams.nWf = 5000;
gwfparams.spikeTimes = ceil(sp.st * sp.sample_rate);
gwfparams.spikeClusters = sp.clu;
gwfparams.newClus = uint32(SUdis(:,1)-1);

wf = getWaveForms(gwfparams);

numSU = size(wf.waveFormsMean,1);
WinSize = size(wf.waveFormsMean,3);
meanwaveForms = nan(numSU,WinSize);
for SUid = 1:size(wf.waveFormsMean,1)
    pp_Ch = max(wf.waveFormsMean(SUid,:,:),[],3) -min(wf.waveFormsMean(SUid,:,:),[],3);
    [~,idx(SUid)] = max(pp_Ch);
    pp_event = max(wf.waveForms(SUid,:,idx(SUid),:),[],4) -min(wf.waveForms(SUid,:,idx(SUid),:),[],4);
    threshold = nanmean(pp_event)+2*nanstd(pp_event);
    pp_event_idx = pp_event < threshold;
    meanwaveForms (SUid,:) = squeeze(mean(wf.waveForms(SUid,pp_event_idx,idx(SUid),:),2));
    [val_min I_min] = min(wf.waveFormsMean(SUid,idx(SUid),:),[],3);
    [val_max I_max] = max(wf.waveFormsMean(SUid,idx(SUid),:),[],3);
    duration(SUid) = I_max - I_min;
    baseline = mean(wf.waveFormsMean(SUid,idx(SUid),1:20),3);
    peaktroughratio(SUid) = (val_max - baseline) /(val_min - baseline);
end
filen1 = ['F344AD_' kilosubj '_SU_Waveform_Output_Extraction.mat'];
disp(['saving SU ' kilosubj])

unitIDs = wf.unitIDs;

save(filen1,'unitIDs','meanwaveForms','idx','duration','peaktroughratio')
disp(['saved SU' kilosubj])

% MU
clear wf; clear SUid; clear unitIDs; clear meanwaveForms; clear idx, clear duration; clear peaktroughratio;
gwfparams.newClus = uint32(MUdis(:,1)-1);

wf = getWaveForms(gwfparams);

numSU = size(wf.waveFormsMean,1);
WinSize = size(wf.waveFormsMean,3);
meanwaveForms = nan(numSU,WinSize);
for SUid = 1:size(wf.waveFormsMean,1)
    pp_Ch = max(wf.waveFormsMean(SUid,:,:),[],3) -min(wf.waveFormsMean(SUid,:,:),[],3);
    [~,idx(SUid)] = max(pp_Ch);
    pp_event = max(wf.waveForms(SUid,:,idx(SUid),:),[],4) -min(wf.waveForms(SUid,:,idx(SUid),:),[],4);
    threshold = nanmean(pp_event)+2*nanstd(pp_event);
    pp_event_idx = pp_event < threshold;
    meanwaveForms (SUid,:) = squeeze(mean(wf.waveForms(SUid,pp_event_idx,idx(SUid),:),2));
    [val_min I_min] = min(wf.waveFormsMean(SUid,idx(SUid),:),[],3);
    [val_max I_max] = max(wf.waveFormsMean(SUid,idx(SUid),:),[],3);
    duration(SUid) = I_max - I_min;
    baseline = mean(wf.waveFormsMean(SUid,idx(SUid),1:20),3);
    peaktroughratio(SUid) = (val_max - baseline) /(val_min - baseline);
end
filen1 = ['F344AD_' kilosubj '_MU_Waveform_Output_Extraction.mat'];
disp(['saving MU ' kilosubj])

unitIDs = wf.unitIDs;

save(filen1,'unitIDs','meanwaveForms','idx','duration','peaktroughratio')
disp(['saved MU ' kilosubj])

% FUNC
function wf = getWaveForms(gwfparams)
fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});
chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);
unitIDs = unique(gwfparams.newClus)';
numUnits = size(unitIDs,2);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(numUnits,nChInMap,wfNSamples);
for curUnitInd=1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    
    parfor curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = mmf.Data.x(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
    end
    waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
    disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
end
% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;
end
