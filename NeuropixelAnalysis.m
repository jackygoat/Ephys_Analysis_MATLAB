
clc
close all
clear all
% %% load the LFP data
% [filen,path] = uigetfile("F:\Neuropixel\20240318_F344AD_4432\2024-03-14_12-41-40\Record Node 101\experiment1\recording1\structure.oebin");
% LFP         = load_open_ephys_binary([path filen],'continuous',2,'mmap');
% OEfac       = LFP.bitVolts;
% freqLFP     = LFP.Header.sample_rate;
% LFPraw_baseline      = double(LFP.Data.Data.mapped(60,:));
% LFPraw_baseline      = LFPraw_baseline.*OEfac;

%% load AP data
[filen,path] = uigetfile("E:\Neuropixel\20240229_F344AD_4414\2024-02-29_14-26-42\Record Node 101\experiment1\recording1\structure.oebin");
AP         = load_open_ephys_binary([path filen],'continuous',1,'mmap');
OEfac       = AP.bitVolts;
freqAP     = AP.Header.sample_rate; %30000 Hz

i = 100;
APraw      = double(AP.Data.Data.mapped(i,:));
APraw      = APraw.*OEfac;
%% Time Stamps
Tsamp = 1:length(APraw);

h=figure;
plot(APraw)
avg = mean(APraw);StdEV = std(APraw);Maxval = max(APraw); Minval = min(APraw);
title(['Average: ',num2str(avg),'  Std: ',num2str(StdEV),'  Max: ',num2str(Maxval),'  Min: ',num2str(Minval)]);
MaxN = abs((Maxval-avg))/StdEV; MinN = abs((Minval-avg))/StdEV;
disp (['Potential N for max: ',num2str(MaxN),'  ; N for min: ',num2str(MinN)])
drawnow;
ThresmaxNum = input('What is the n for max threshold? ');
hold on
yline(avg+ThresmaxNum*StdEV,'-.r')
drawnow;
change = input('want to change threshold? '); 
if change ~= 0
    ThresmaxNum = change;
    yline(avg+ThresmaxNum*StdEV,'k','LineWidth',1)
end
BaselineT = input('Identify Baseline Time start: ');
% ThresminNum = input('What is the n for min threshold? ');
% hold on
% yline(avg-ThresminNum*StdEV,'-.k')

posthres = Tsamp(APraw > (avg+ThresmaxNum*StdEV))';
% tmp_negthres = Tsamp(APraw < (avg-ThresminNum*StdEV))';
testpos = diff(posthres);
% test2 = diff(tmp_negthres);
idx = find(testpos>(freqAP/8-100));


if testpos(idx(1)) > 20000 & testpos(idx(1)) < 8900000
    Stamp(1) = posthres(idx(1)+1);
else
    Stamp(1) = posthres(1);
end

ct =2;
for i = 2:length(testpos(idx))
    if testpos(idx(i)) > 590000
        Stamp(ct) = posthres(idx(i)+1);
        ct = ct+1;
    else testpos(idx(i)) > 3600 & testpos(idx(i)) < 8900000;
        if i+1 < length(testpos(idx)) ;
            if testpos(idx(i)) + testpos(idx(i+1)) > 8900000 & testpos(idx(i+1)) < (8900000-1000) ;
            Stamp(ct) = posthres(idx(i+1)+1);
            ct = ct+1;
            continue
            elseif  testpos(idx(i)) + testpos(idx(i+1)) + testpos(idx(i+2))> 8900000 & all([testpos(idx(i)) testpos(idx(i+1)) testpos(idx(i+2))] < 8900000);
                            Stamp(ct) = posthres(idx(i+2)+1);
                            ct = ct+1;
                            continue
                            continue
            end
        end
    end
    
    if i > length(testpos(idx))
        break
    end
end


TrialStamp = unique(Stamp);
idx_pulse = find(testpos < (freqAP/3+1000) & testpos > 3000);
for loop = 1:length(TrialStamp)
    firstnumber = find(idx_pulse>=(find(posthres == TrialStamp(loop))),1,'first');
    if loop <= length(TrialStamp)-1;
    Trial_pulse(1,loop) = TrialStamp(loop);
    lastnumber = find(idx_pulse<(find(posthres == TrialStamp(loop+1))),1,'last');
    Trial_pulse(2:(lastnumber-firstnumber+2),loop) = posthres(idx_pulse(firstnumber:lastnumber)+1);
    elseif loop == length(TrialStamp);
%     Trial_pulse(2:end, loop) = posthres(idx_pulse(firstnumber:end)+1);
    lastseq = posthres(idx_pulse(firstnumber:end)+1);
    end
end
Trial_pulse(Trial_pulse == 0) = nan;
BaselineT(2) = BaselineT(1) + Trial_pulse(end,1)+Trial_pulse(2,1)-2*Trial_pulse(1,1);
cd(path)
StrIx0 = strfind(path,'AD_');subj = extractAfter(path,'AD_');subj = subj(1:4);
filen1 = ['F344AD_' subj '_TimeStamps.mat'];
save(filen1,'Trial_pulse','BaselineT','lastseq')
% 'lastseq'

%% CSD plot (frequencies: 1-120Hz;  3.5-7.5Hz; 30-100Hz)
clear all
clc
close all


kilopath = uigetdir("/Volumes/Jack's Gate/Sunnybrook/Data/20240314_F344AD_4402/2024-03-14_12-41-40/Record Node 101/experiment1/recording1/continuous/Neuropix-PXI-108.ProbeA-AP/");
cd(kilopath)
% kilosubj = extractAfter(kilopath,'AD_');kilosubj = kilosubj(1:4);
kilosubj = '3928';
load(['F344AD_' kilosubj '_TimeStamps.mat']) % load time stamp
load('chanMap_ChannelMaP_KC') %load channel map

[filen,path] = uigetfile([kilopath(1:93) '\structure.oebin']);
LFP         = load_open_ephys_binary([path filen],'continuous',2,'mmap');
OEfac       = LFP.bitVolts;
freqLFP     = LFP.Header.sample_rate;
time_baseline = round(BaselineT/fs *freqLFP);
extratime = 2500; %1s
LFPraw_baseline      = double(LFP.Data.Data.mapped(:,time_baseline(1)-extratime:time_baseline(2)+extratime));
LFPraw_baseline      = LFPraw_baseline.*OEfac;

for loop = 1:size(Trial_pulse,2)
    Pulsetmp = Trial_pulse(:,loop);
    tmp = Pulsetmp(~isnan(Pulsetmp));
    time_stim = ceil([tmp(1) 2*tmp(end)-tmp(end-1)]/fs*freqLFP);
    LFPraw_stim      = double(LFP.Data.Data.mapped(:,time_stim(1)-extratime:time_stim(2)+extratime));
    LFPraw_stim      = LFPraw_stim.*OEfac;

    a = sort(ycoords,'ascend');
    [ycoords_sorted,sortIdx] = sort(ycoords,'descend');
    LFPbaseline_chansorted = LFPraw_baseline(sortIdx,:);
    LFPstim_chansorted = LFPraw_stim(sortIdx,:);

% average 4 sites 
LFPbaseline_avg=[];LFPstim_avg = [];csd_d = [];
for k = 1:(size(LFPraw_baseline,1)/4)
    LFPbaseline_avg = [LFPbaseline_avg;mean(LFPbaseline_chansorted(4*(k-1)+1:4*k,:))];
    LFPstim_avg = [LFPstim_avg;mean(LFPstim_chansorted(4*(k-1)+1:4*k,:))];
    csd_d = [csd_d;mean(ycoords_sorted(4*(k-1)+1:4*k))];
end
csd_d_cor = csd_d-min(csd_d);

carsig_baseline= max(LFPbaseline_avg);
LFPbaseline_avg= LFPbaseline_avg -carsig_baseline;
carsig_stim = max(LFPstim_avg);
LFPstim_avg = LFPstim_avg-carsig_stim;


% filter the LFP
csdlowfreq=1;csdhighfreq= 120;
[bcsd, acsd] = butter(1, [csdlowfreq csdhighfreq]/(freqLFP*0.5));
for channel = 1:size(LFPstim_avg,1)
    csdLFPbaseline_Filtered{loop}(channel,:) = filter(bcsd,acsd,double(LFPbaseline_avg(channel,:)), [],2);
    csdLFPstim_Filtered{loop}(channel,:) = filter(bcsd,acsd,double(LFPstim_avg(channel,:)), [],2);
end


lfpmapt = size(csdLFPstim_Filtered{loop},2);
lintime=linspace(1,lfpmapt,lfpmapt)/freqLFP;
csd_t = lintime;
% for channel=2:size(LFPstim_Filtered,1)
%     csd_t=cat(1,csd_t,lintime);
% end

[tmpgx,tmpgy]=gradient(csdLFPstim_Filtered{loop}',csd_d_cor,csd_t);
[tmpgx2,tmpgy]=gradient(tmpgx,csd_d_cor,csd_t);
tmpgx3=tmpgx2';

startTimepoint = time_stim(1)-extratime;
stimtrigcsd=round(tmp/fs*freqLFP)-startTimepoint;
cb = 0.12;% seconds after or before trig
ca = 0.3;
csdon=tmpgx3(:,round(stimtrigcsd(1)-cb*freqLFP):round(stimtrigcsd(1)+ca*freqLFP));
firststimulicsd(:,:,loop) = csdon;
lfpon=csdLFPstim_Filtered{loop}(:,round(stimtrigcsd(1)-cb*freqLFP):round(stimtrigcsd(1)+ca*freqLFP));

for csdcnt = 2:length(stimtrigcsd)
    csdtmp=tmpgx3(:,round(stimtrigcsd(csdcnt)-cb*freqLFP):round(stimtrigcsd(csdcnt)+ca*freqLFP));
    lfptmp=csdLFPstim_Filtered{loop}(:,round(stimtrigcsd(csdcnt)-cb*freqLFP):round(stimtrigcsd(csdcnt)+ca*freqLFP));
    csdon=csdon+csdtmp;
    lfpon=lfpon+lfptmp;
end
csdon_tot(:,:,loop) = csdon;
lfpon_tot(:,:,loop) = lfpon;


% csdon_avg = mean(csdon_tot,3);
% StrIx0 = strfind(path,'AD_');subj = extractAfter(path,'AD_');subj = subj(1:4);
% cd('C:\Users\skytechpc1\OneDrive\PostDoc\AD_NEI\Prelim Analysis\04162024_meeting')
% f = figure('Name',[subj '_CSD Repetition ' num2str(loop)]);
% lintime=linspace(-cb*1000,ca*1000,size(csdon,2));
% h = pcolor(lintime,csd_d_cor,csdon);
% set(h,'EdgeColor','none');shading interp;
% set(gca, 'YDir','reverse')
% colormap jet
% xlim([-100 300]);
% caxis([-0.15 0.15])
% xlabel('Time (ms)')
% ylabel('Depth (\mum)')
% % saveas(f,[f.Name '.fig'])

clear csdon carsig_baseline carsig_stim
end

cd(kilopath)
filen0 = [kilosubj '_CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz band plotting'];
save(filen0,'firststimulicsd','csdon_tot','lfpon_tot','csd_t','csd_d_cor')
%
csdon_avg = mean(csdon_tot,3);
f1 = figure('Name',[kilosubj '_CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz theta band repetition avged']);
lintime=linspace(-cb*1000,ca*1000,size(csdon_avg,2));
h1 = pcolor(lintime,csd_d_cor,csdon_avg);
set(h1,'EdgeColor','none');shading interp;
set(gca, 'YDir','reverse')
colormap jet
xlim([-50 250]);
caxis([-0.06 0.06])
xlabel('Time (ms)')
ylabel('Depth (\mum)')
title([kilosubj '_CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz Repetition-avg trial-avg CSD'])
saveas(f1,[f1.Name '.fig'])

%
% for ii = 1:10
f2 = figure('Name',[kilosubj '_CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz band First Stimuli']);
csdon1_avg = mean(firststimulicsd(:,:,:),3);
lintime=linspace(-cb*1000,ca*1000,size(csdon1_avg,2));
% h = pcolor(lintime,csd_d_cor,firststimulicsd(:,:,ii));
h = pcolor(lintime,csd_d_cor,csdon1_avg);
set(h,'EdgeColor','none');shading interp;
set(gca, 'YDir','reverse')
colormap jet
xlim([-100 300]);
xlabel('Time (ms)')
ylabel('Depth (\mum)')
title([kilosubj '_CSD ' num2str(csdlowfreq) '-' num2str(csdhighfreq) 'Hz band First Stimuli Trial CSD'])
caxis([-0.012 0.012])
xlim([-50 250])
saveas(f2,[f2.Name '.fig'])
% end
%% Power Spectrum Analysis
% % % % % figure
% % % % % tiledlayout('flow')
% % % % % nexttile
% % % % % pspectrum(mean(LFPstim_avg_raw(:,:),1),freqLFP,"power","FrequencyLimits",[1 300])
% % % % % nexttile
% % % % % pspectrum(mean(LFPstim_avg{loop}(:,:),1),freqLFP,"power","FrequencyLimits",[1 300])
% % % % % 
% % % % % 
% % % % % figure
% % % % % tiledlayout('flow')
% % % % % nexttile
% % % % % pspectrum(mean(LFPstim_avg_raw(:,:),1),freqLFP,"spectrogram","FrequencyLimits",[1 120])
% % % % % %  caxis([35 75])
% % % % % nexttile
% % % % % pspectrum(mean(LFPstim_avg{loop}(:,:),1),freqLFP,"spectrogram","FrequencyLimits",[1 120])
% % % % % %  caxis([35 75])

%% Power Spectrum Analysis - 1min 

clear all
clc
close all


kilopath = uigetdir("E:\Neuropixel\20240229_F344AD_4414\2024-02-29_14-26-42\Record Node 101\experiment1\recording1\continuous\Neuropix-PXI-108.ProbeA-AP\");
cd(kilopath)
kilosubj = extractAfter(kilopath,'AD_');kilosubj = kilosubj(1:4);
load(['F344AD_' kilosubj '_TimeStamps.mat']) % load time stamp
load('chanMap_ChannelMaP_KC') %load channel map

[filen,path] = uigetfile([kilopath(1:93) '\structure.oebin']);
LFP         = load_open_ephys_binary([path filen],'continuous',2,'mmap');
OEfac       = LFP.bitVolts;
freqLFP     = LFP.Header.sample_rate;
% time_baseline = round(BaselineT/fs *freqLFP);
time_baseline = [2833,156515];
extratime = 2500; %1s
LFPraw_baseline      = double(LFP.Data.Data.mapped(:,time_baseline(1)-extratime:time_baseline(2)+extratime));
LFPraw_baseline      = LFPraw_baseline.*OEfac;

for loop = 1:size(Trial_pulse,2)
    Pulsetmp = Trial_pulse(:,loop);
    tmp = Pulsetmp(~isnan(Pulsetmp));
    time_stim = ceil([tmp(1) 2*tmp(end)-tmp(end-1)]/fs*freqLFP);
    LFPraw_stim      = double(LFP.Data.Data.mapped(:,time_stim(1)-extratime:time_stim(2)+extratime));
    LFPraw_stim      = LFPraw_stim.*OEfac;

%     a = sort(ycoords,'ascend');
    [ycoords_sorted,sortIdx] = sort(ycoords,'descend');
    LFPbaseline_chansorted = LFPraw_baseline(sortIdx,:);
    LFPstim_chansorted = LFPraw_stim(sortIdx,:);

% average 4 sites 
LFPbaseline_avg_raw=[];LFPstim_avg_raw = [];csd_d = [];
for k = 1:(size(LFPraw_baseline,1)/4)
    LFPbaseline_avg_raw= [LFPbaseline_avg_raw;mean(LFPbaseline_chansorted(4*(k-1)+1:4*k,:))];
    LFPstim_avg_raw = [LFPstim_avg_raw;mean(LFPstim_chansorted(4*(k-1)+1:4*k,:))];
    csd_d = [csd_d;mean(ycoords_sorted(4*(k-1)+1:4*k))];
end
csd_d_cor = csd_d-min(csd_d);

carsig_baseline= max(LFPbaseline_avg_raw);
LFPbaseline_avg{loop} = LFPbaseline_avg_raw-carsig_baseline;
carsig_stim = max(LFPstim_avg_raw);
LFPstim_avg{loop} = LFPstim_avg_raw-carsig_stim;

for channel = 1:(size(LFPraw_baseline,1)/4)
[p_stim(channel,:,loop),f1,t1] = pspectrum(LFPstim_avg{loop}(channel,:),freqLFP,"Power","FrequencyLimits",[1 120]);
[p_baseline(channel,:,loop),f2,t2] = pspectrum(LFPbaseline_avg{loop}(channel,:),freqLFP,"Power","FrequencyLimits",[1 120]);
nor_Power(channel,:,loop) = 10*log10(p_stim(channel,:,loop)) - 10*log10(p_baseline(channel,:,loop));
end

% [Spectrogram_stim(:,:,channel,loop),f11,t11] = pspectrum(LFPstim_avg{loop}(channel,:),freqLFP,"spectrogram","FrequencyLimits",[1 120]);
% [Spectrogram_baseline(:,:,channel,loop),f22,t22] = pspectrum(LFPbaseline_avg{loop}(channel,:),freqLFP,"spectrogram","FrequencyLimits",[1 120]);
% nor_PowerSpectrogram(:,:,channel,loop) = 10*log10(Spectrogram_stim(:,:,channel,loop) ) - 10*log10(Spectrogram_baseline(:,:,channel,loop));
end

%%% Power Spectrum _ plotting check
% % % % % % figure
% % % % % % tiledlayout('flow')
% % % % % % nexttile
% % % % % % pspectrum(mean(LFPstim_avg_raw(:,:),1),freqLFP,"power","FrequencyLimits",[1 300])
% % % % % % nexttile
% % % % % % pspectrum(mean(LFPstim_avg{loop}(:,:),1),freqLFP,"power","FrequencyLimits",[1 300])
% % % % % % 
% % % % % % 
% % % % % % figure
% % % % % % tiledlayout('flow')
% % % % % % nexttile
% % % % % % pspectrum(mean(LFPstim_avg_raw(:,:),1),freqLFP,"spectrogram","FrequencyLimits",[1 120])
% % % % % % %  caxis([35 75])
% % % % % % nexttile
% % % % % % pspectrum(mean(LFPstim_avg{loop}(:,:),1),freqLFP,"spectrogram","FrequencyLimits",[1 120])
% % % % % % %  caxis([35 75])

%%% power spectrum_ data save
StrIx0 = strfind(path,'AD_');subj = extractAfter(path,'AD_');subj = subj(1:4);
fig1 = figure('Name',[subj '_Normalized Power Magnitude']);
pcolor(f1,csd_d_cor,mean(nor_Power,3));set(gca, 'YDir','reverse');shading interp;colorbar;cb = colorbar();yl = ylabel(cb,'Normalized Power (dB)');
xlabel('Frequency (Hz)')
ylabel('Depth (\mum)')
title([subj ' Normalized Power Magnitude'])
ax = gca;ax.FontSize = 12;

cd(kilopath)
filen0 = [subj '_Power Spectrum analysis'];
save(filen0,'LFPbaseline_avg','LFPstim_avg','p_stim','p_baseline','nor_Power')
saveas(fig1,[fig1.Name '.fig'])
% Power spectrogram does not look very interesting
% fig2 = figure('Name',[subj '_Normalized Power Magnitude Spectrogram']);
% imagesc(t11,f11,mean(mean(nor_PowerSpectrogram,3),4));set(gca, 'YDir','reverse');shading interp;colorbar;cb = colorbar();yl = ylabel(cb,'Power (dB)');
% xlabel('Time (ms)')
% ylabel('Frequency (Hz)')
% title([subj ' Normalized Power Magnitude Spectrogram'])
% ax = gca;ax.FontSize = 12;