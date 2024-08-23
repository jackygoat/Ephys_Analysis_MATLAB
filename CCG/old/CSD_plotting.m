%% LFP 
clear all
clc
% load raw LFP data
kilopath = uigetdir();
cd(kilopath)
kilosubj = extractAfter(kilopath,'Cohort'); kilosubj = kilosubj(2:5);
% kilosubj = extractAfter(kilopath,'AD_'); kilosubj = kilosubj(1:4);
% kilosubj = extractAfter(kilopath,'Emma\');kilosubj = kilosubj(10:13);
load(['F344AD_' kilosubj '_HPC_TimeStamps.mat']) % load time stamp
% load('chanMap_Site1_ChanMap_v2_08062024')  
load('neuropixPhase3B2_kilosortChanMap.mat')

% load('chanMap_100sitesAway_NP1') 
% ycoords = sort(ycoords,'descend');
% input1 = input('channel map');
% if input1 == 0;
% ycoords = sort(ycoords,'descend');
% end
    %load channel map

[filen,path] = uigetfile([kilopath(1:93) '\structure.oebin']);
LFP         = load_open_ephys_binary([path filen],'continuous',2,'mmap');
OEfac       = LFP.bitVolts;
freqLFP     = LFP.Header.sample_rate;
% check and remove malfunc channel, merge 4 sites
pickRandom  = randi([1, Trial_pulse(1,1)-10*fs]);
baselineTtmp = pickRandom:pickRandom+10*fs; 
range1      = round(baselineTtmp/fs *freqLFP);
LFPraw_baselinetmp      = double(LFP.Data.Data.mapped(:,range1(1):range1(end)));
LFPraw_baselinetmp      = LFPraw_baselinetmp.*OEfac;
rmssig1     = rms(LFPraw_baselinetmp,2);
rngNoise    = prctile(rmssig1,[0.5,99.5]);
indxChn     = rmssig1>rngNoise(1)&rmssig1<rngNoise(2);
[ycoords_sorted,sortIdx] = sort(ycoords);
chnKeep     = indxChn(sortIdx);
nChanAvg    = 4;chanOverlap = 0;
totNPchans  = size(LFP.Data.Data.mapped,1);
extratime = 5000;
for loop = 1:size(Trial_pulse,2)
    Pulsetmp = Trial_pulse(:,loop);
    tmp = Pulsetmp(~isnan(Pulsetmp));
    time_stim = ceil([tmp(1) 2*tmp(end)-tmp(end-1)]/fs*freqLFP);
    % BaselineT_total = [Trial_pulse(1,1)-6*60*fs Trial_pulse(1,1)];
    % time_baseline = round(BaselineT_total/fs *freqLFP);
        newChanNum  = floor(totNPchans/(nChanAvg-chanOverlap));
    %     tVec        = 0:(1/freqLFP):(1/freqLFP)*(time_stim(2)-time_stim(1));
    %     LFPd        = zeros(newChanNum,length(tVec));
        LFPraw_baseline = double(LFP.Data.Data.mapped(:,time_stim(1)-extratime:time_stim(2)+extratime));
    for nc = 1:newChanNum
            index1     = (nc-1)*(nChanAvg-chanOverlap)+1:(nc-1)*(nChanAvg-chanOverlap)+nChanAvg;
            if max(index1) > totNPchans
                break
            end
            logicKeep  = chnKeep(index1);                                      % discard corrupted chans
            chns2avg   = sortIdx(index1);
            depth(nc) = mean(ycoords(chns2avg));
            LFPd(nc,:) = mean(LFPraw_baseline(chns2avg(logicKeep),:).*OEfac,1);
    end
        % 60 Hz Notch Filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',freqLFP);
    for chan = 1:size(LFPd,1)
        LFPd_notch(chan,:) = filtfilt(d,LFPd(chan,:));
    end
    
    % CAR filter 
    carsig= median(LFPd_notch);
    LFPbaseline_carsigfilter = LFPd_notch-carsig;
    
    % CSD Bandpass filter
    csdlowfreq=2;csdhighfreq= 120;
    [bcsd, acsd] = butter(1, [csdlowfreq csdhighfreq]/(freqLFP*0.5));
    for ch = 1:size(LFPd,1)
        LFP_Filtered{loop}(ch,:) = filter(bcsd,acsd,double(LFPbaseline_carsigfilter(ch,:)), [],2);
    end
    
    
    lfpmapt = size(LFP_Filtered{loop},2);
    lintime=linspace(1,lfpmapt,lfpmapt)/freqLFP;
    csd_t = lintime;
    
    [tmpgx,tmpgy]=gradient(LFP_Filtered{loop}',1:96,csd_t);
    [tmpgx2,tmpgy]=gradient(tmpgx,1:96,csd_t);
    tmpgx3=tmpgx2';
    
    startTimepoint = time_stim(1)-extratime;
    stimtrigcsd=round(tmp/fs*freqLFP)-startTimepoint;
    cb = 1.5;% seconds after or before trig
    ca = 1.5;
    csdon=tmpgx3(:,round(stimtrigcsd(1)-cb*freqLFP):round(stimtrigcsd(1)+ca*freqLFP));
    firststimulicsd(:,:,loop) = csdon;
    lfpon=LFP_Filtered{loop}(:,round(stimtrigcsd(1)-cb*freqLFP):round(stimtrigcsd(1)+ca*freqLFP));
    
%     for csdcnt = 2:length(stimtrigcsd)
%         csdtmp=tmpgx3(:,round(stimtrigcsd(csdcnt)-cb*freqLFP):round(stimtrigcsd(csdcnt)+ca*freqLFP));
%         lfptmp=LFP_Filtered{loop}(:,round(stimtrigcsd(csdcnt)-cb*freqLFP):round(stimtrigcsd(csdcnt)+ca*freqLFP));
%         csdon=csdon+csdtmp;
%         lfpon=lfpon+lfptmp;
%     end
%     csdon_tot(:,:,loop) = csdon;
%     lfpon_tot(:,:,loop) = lfpon;
% 
    clear LFPd LFPd_notch
end
% First Trial

f1 = figure('Name',[kilosubj ' CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz First Stimuli']);
csdon1_avg = mean(firststimulicsd(:,:,:),3);
cb = 1.5;% seconds after or before trig
ca = 1.5;
lintime=linspace(-cb*1000,ca*1000,size(csdon1_avg,2));
% h = pcolor(lintime,1:95,firststimulicsd(:,:,ii));
h = pcolor(lintime,flip(depth),csdon1_avg);
set(h,'EdgeColor','none');shading interp;
set(gca, 'YDir','reverse')
colormap jet
xlim([-1500 1500]);
xlabel('Time (ms)')
ylabel('Depth(µm)')
title([kilosubj ' CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz First Stimuli Trial'])
caxis([-20 20])
saveas(f1,[f1.Name '.fig'])