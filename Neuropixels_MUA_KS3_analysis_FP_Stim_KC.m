clc; clear all; close all

mainDir = uigetdir();
cd(mainDir)
StrIx0      = strfind(mainDir,'AD_');
% subj        = extractAfter(mainDir,'AD_');
% subj        = subj(1:4);
subj        = '4402_HPC';
oldDir      = cd(mainDir);
load('rez.mat')                                                         % load KS3 results
% load stim times from LFP-analysis folder
pathDirs    = strsplit(mainDir,'\');
Indx        = find(contains(pathDirs,'recording'));
oebinDir    = fullfile(pathDirs{1:Indx});
alldirs     = dir(oebinDir);
subdirs0    = alldirs([alldirs(:).isdir]);
subdirs0    = subdirs0(~ismember({subdirs0(:).name},{'.','..'}));
subdirs0    = struct2cell(subdirs0);
subdirs0    = subdirs0(1,:);

% load(['F344AD_' subj '_TimeStamps.mat'])

% cluster analysis
clustN      = unique(rez.st3(:,2));
vNames      = {'spike_sample','clusterN','Amplitude','X','Y'};
xypos1      = flip(rez.xy,2);
xypos1(:,2) = abs(xypos1(:,2)-max(xypos1(:,2)));                            % KS has this weird thing of calc depth as 'distance from tip'
KSout       = array2table([rez.st3(:,1:3) xypos1],'VariableNames',vNames);
fs          = rez.ops.fs;
tEnd0       = rez.ops.tend;
Good        = rez.good;
tvec0       = 0:1/fs:(1/fs)*(tEnd0-1);

% set some metrics and discard some units
minFR       = 0.1;                                                        % min firing rate for each cluster
minCovr     = 0.75;                                                     % min coverage of spikes over the tot time of recording 
[rngAmp] = prctile(KSout.(vNames{3}), [0.1, 92]);

clusTab     = table();                                                  % clusters to keep
exclClu     = table();                                                  % discarded clusters

SUclus      = table();
MUclus      = table();

for cl = 1:max(clustN)
    clInd   = KSout.(vNames{2})== cl;
    tmpTcl  = KSout(clInd,:); 
    amplit  = tmpTcl.(vNames{3});
    avgAmp  = mean(amplit);
    if avgAmp > rngAmp(1) && avgAmp < rngAmp(2)
        
        % [N,edges] = histcounts(amplit, 'Normalization', 'probability','BinWidth',1);
        [f,xi]  = ksdensity(amplit,'Bandwidth',2);                          % calculate the probab density function
        [idx0,Prom] = islocalmax(f,'SamplePoints',xi,'MinProminence',0.02); % find main peaks in the data
        if sum(idx0)>1
            Prom2 = sort(Prom,'descend');
            idx0 = Prom == Prom2(1);                                        % look for highest peak if more than 1 peak is found
        end
        pdfDer1 = abs(diff(f));
        RminInd = find(idx0) + find(islocalmin(pdfDer1(find(idx0):end)), 1, 'first') +1;
        if isempty(RminInd)
            RminInd = length(f);
        end
        LminInd = find(islocalmin(pdfDer1(1:find(idx0))), 1, 'first');
        if isempty(LminInd)
            LminInd = 1;
        end
        ampIdx  = amplit>xi(LminInd) & amplit<xi(RminInd);
        tmpTcl  = tmpTcl(ampIdx,:);                                         % keep only spike amplitudes that fall within distribution
        tmpGood = table(repmat(logical(Good(cl)),height(tmpTcl),1),'VariableNames',{'Good'});
        tmpTcl  = cat(2,tmpTcl, tmpGood);
        % [h,p,ad_stat] = adtest(tmpTcl.(vNames{3}),'Distribution','norm','Alpha',0.01);                         % check if amplitude is normally distributed, if it is, h =0 and p>0.05
        
        spikeT  = tmpTcl.(vNames{1}).*fs^-1;
        avgFR1  = 1/mean(diff(spikeT));                                     % recalc avg firing and amplitude
        avgAmp  = median(tmpTcl.(vNames{3})); 
        avgFR   = repmat(avgFR1,height(tmpTcl),1);
        avgAmpl = repmat(avgAmp,height(tmpTcl),1);
        tmpTcl  = addvars(tmpTcl,avgFR,avgAmpl,'After',vNames{3},'NewVariableNames',{'avgFR','avgAmplitude'});
        spCover = (max(spikeT)-min(spikeT))/tvec0(end);                     % want to have >75% coverage during the recording for the unit!
        if avgFR1 > minFR && spCover>= minCovr % && ~h
            clusTab = cat(1,clusTab,tmpTcl);
                if Good(cl) == 1
                    SUclus = cat(1,SUclus,tmpTcl);
                elseif Good(cl) == 0
                    MUclus = cat(1,MUclus,tmpTcl);
                end
        else
            tTable  = table(cl,avgFR1,spCover);
            exclClu = [exclClu; tTable];
        end
    end

end
% end

% save filtered clusters
cd(mainDir)
filen0 = [subj '_FilteredClusters.mat'];
% save(filen0,'clusTab','exclClu')

f1  = figure('Name',[subj '_clusterPos'],'NumberTitle','off','color','w');
vNames2 = clusTab.Properties.VariableNames;     % same as SUclus and MUclus
maxFR   = max(clusTab.(vNames2{4}));
minFR   = min(clusTab.(vNames2{4}));
newClus = unique(clusTab.(vNames{2}))';newSUClus = unique(SUclus.(vNames{2}))';newMUClus = unique(MUclus.(vNames{2}))';

clusY   = zeros(numel(newClus),2);
avgClFR = zeros(numel(newClus),3);
SUdis = zeros(numel(newSUClus),3);MUdis = zeros(numel(newMUClus),3);

n = 1;nsu =1; nmu =1;
for i = newClus
    indx0   = clusTab.(vNames{2}) == i;
    tmpTcl1 = clusTab(indx0,:);
    avgFR2  = unique(tmpTcl1.(vNames2{4}));
    clPosX  = tmpTcl1.(vNames{4});
    clPosX1 = mean(clPosX,'omitnan'); 
    clPosY  = tmpTcl1.(vNames{5});
    clPosY1 = mean(clPosY,'omitnan');
    if ~any(tmpTcl1.(vNames2{8}))
        MUdis(nmu,:) = [i, clPosY1,avgFR2];
        nmu = nmu +1;
    else
        SUdis(nsu,:) = [i, clPosY1,avgFR2];
        nsu = nsu +1;
    end
    figure(f1); hold on
    sz      = avgFR2/maxFR*200; % cluster size refecting firing rate of cluster? 
    scatter(clPosX1/1000,clPosY1/1000,sz,'filled') % y and x position
    clusY(n,:) = [i, clPosY1];                                              % save cluster depth
    avgClFR(n,:) = [i,avgFR2,avgFR2/maxFR];                                 % and average firing rate
    n = n +1;
    % axis('equal')
end
xlim([0.01 0.06])
ylabel('depth (mm)')
xlabel('mm')
set(gca, 'YDir','reverse')
saveas(f1,[f1.Name '.fig'])

binwidth = 150; edges = 50:binwidth:4050;
SUcounts = histcounts(SUdis(:,2),edges);MUcounts = histcounts(MUdis(:,2),edges);
Barcenter = edges(1:end-1) + binwidth/2;

cd(mainDir)
filen1 = [subj '__FilterCluster_SUMU'];
save(filen1,'clusTab','exclClu','SUclus','MUclus','newClus','newSUClus','newMUClus','SUdis','MUdis','clusY','clPosX1')
%
%%
baseT   = BaselineT./fs;                                                               % time in sec to check FR of the cluster prior to stim
stblkn  = 1:size(Trial_pulse,2);
n = 1;nmu = 1; nsu = 1;
for l = newClus
        indx1   = clusTab.(vNames2{2}) == l;
        tmpTcl2 = clusTab(indx1,:);
        spkT    = tmpTcl2.(vNames2{1}); 
        spkT    = sort(spkT./fs);
        clusterISI  = diff(sort(spkT));
                                             
                % null distribution
                nper = 1e4;
                nullISIdelta = zeros(nper,1);
                parfor nperm = 1:nper
                    % shuffle Firing Rates
                    shuffindex2     = randperm(numel(clusterISI));                  % shuffle firing rates
                    ISIShuff2       = clusterISI(shuffindex2);
                    indxShuf2       = floor(numel(ISIShuff2) /2);                   % divide the shuffled ISI in 2 consecutive blocks
                    tmpBaseShISI    = ISIShuff2(1:indxShuf2);
                    tmpStimShISI    = ISIShuff2(indxShuf2+1:indxShuf2*2);
                    TmpISImed       = [median(tmpBaseShISI),median(tmpStimShISI)];  % calc median FR for all
                    nullISIdelta(nperm,:) = TmpISImed(2)-TmpISImed(1);              % pseudo stim - baseline delta for ISI
                end
                baseISI         = diff(spkT(spkT>baseT(1) & spkT<=baseT(2)));                 % inter-spike intervals at baseline
                baseFR          = 1./mean(baseISI);
                baseStd         = 1./std(baseISI);

        for j = stblkn
            zapBlok = Trial_pulse(:,j);
            nStims  = sum(~isnan(zapBlok));                
            stimT = zapBlok./fs;
            dt      = median(diff(stimT),'omitnan');
            upperR  = floor(dt*10)/10;
            lowerR  = floor((dt-upperR)*10)/100;
            spkT1           = spkT(spkT>stimT(1) & spkT<=stimT(end));                             % spike times following the stim
            stimISI         = diff(spkT1);
            ISIstimDelt     = median(stimISI)-median(baseISI);                      % delta ISI stim vs baseline
            pValISIstim     = sum(abs(nullISIdelta(:,1)) >= abs(ISIstimDelt)) / size(nullISIdelta,1);

                if pValISIstim < 0.05
                    SigResRep (n,1) = l;
                    SigResRep (n,j+1) = 1;
                end


        for k = 1:nStims                                                    % loop through the blocks
            t0      = zapBlok(k)/fs;
                           % lower bound
            rng0    = [t0-lowerR, t0+upperR];                               % range to look for stim in the spike times vec
            spkT2   = spkT(spkT>t0 & spkT<rng0(2));                         % spike times following the stim

            if length(spkT2) >=2
                stimISIindTrial         = diff(spkT2);
                ISIstimindTrialDelt     = median(stimISIindTrial)-median(baseISI); 
                pValTrial     = sum(abs(nullISIdelta(:,1)) >= abs(ISIstimindTrialDelt)) / size(nullISIdelta,1);
                if pValTrial < 0.05
                    SigResTrial_ID(n) = l;
                    SigResTrial (n,k,j) = 1;
                      if ~any(tmpTcl2.(vNames2{8}))
                          MUResTrial_ID(n) = 1;
                          MUResTrial (n,k,j) = 1;
                      else
                          SUResTrial_ID(n) = 1;
                          SUResTrial (n,k,j) = 1;
                      end
                end
            end
        end
        end
        n = n+1;
end


cd(mainDir)
filen0 = [subj '_stochastic responsiveness '];
save(filen0,'SigResTrial_ID','SigResTrial','MUResTrial_ID','MUResTrial','SUResTrial_ID','SUResTrial')

%     figure;pcolor(standratio);set(gca, 'YDir','reverse');shading interp;colorbar; xlabel('Trial Repetition'); ylabel('Cluster #')