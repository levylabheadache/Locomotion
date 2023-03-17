%% Fit a master HMM using all available data
trainMouse = 'DL89';
xMouse = strcmpi(expt.mouse, trainMouse)
[transEst, emitEst, speedDiscreteExpt] = TrainLocoHMM(loco, 'dir',sprintf('%s%s', mainDir, );

%% Find locomotive and CSD bouts
sigThresh = 0.01;
boutSumm = repmat( struct('totBout',NaN), 1, Nexpt ); onStruct = struct(); 
effectPool = cell2struct( repmat({struct('short',[], 'long',[])}, NallVars, 1), allVars, 1);
periBout = cell(1,Nexpt); periStat = cell(1,Nexpt); periParam = cell(1,Nexpt); 
csdBout = cell(1,Nexpt); csdStat = cell(1,Nexpt); %periParam = cell(1,Nexpt); 
%clearvars csdBout csdStat; %csdBout = []; % repmat(struct(), 1, Nexpt);
xCSDbout = []; minBout = 10;
tailVal = 'right';
tic
[loco, transEst, emitEst, ~] = GetLocoState(loco, 'show',false);
%xCSD = find([expt.csd]>0);
for x = xPresent
    % Get peri-locomotion and peri-CSD bout data
    for run = 1:expt(x).Nruns % flip()
        [periBout{x}(run), periStat{x}(run), periParam{x}(run), loco{x}(run)] = ...
            PeriLoco3D(runInfo{x}(run), Tscan{x}{run}, loco{x}(run), deform{x}(run), fluor{x}(run).z.ROI, defVars, 'base',5, 'run',4, 'on',4, 'Nshuff',10); %
        if ismember(x, xCSD) && run == expt(x).csd
            [csdBout{x}(run), csdStat{x}(run)] = ...
                PeriCSD(runInfo{x}(run), Tscan{x}{run}, loco{x}(run), deform{x}(run), fluor{x}(run), defVars, periParam{x}(run), 'show',false); 
        end
    end
    boutSumm(x).totBout = sum([periBout{x}.Nbout]);
    boutSumm(x).Tstart = [periBout{x}.Tstart]';
    boutSumm(x).dur = [periBout{x}.dur]';
    
    for v = 1:NallVars
        if boutSumm(x).totBout > 0
            % How many bouts drove deformation/activity in at least one ROI/plane?
            tempVarStat = [periStat{x}.(allVars{v})];
            tempEffect = vertcat( tempVarStat.effect );
            boutSumm(x).(allVars{v}).effect = tempEffect(:,:,3); % 1 = onset effect, 2 = offset effect, 3 = overall effect (see BoutEffect)
            tempPvar = vertcat(tempVarStat.pEffect);
            boutSumm(x).(allVars{v}).effectSig = tempPvar(:,:,3) < sigThresh & ~isnan(boutSumm(x).(allVars{v}).effect) & abs(boutSumm(x).(allVars{v}).effect) > 0;
            boutSumm(x).(allVars{v}).boutSig = find( any( tempPvar(:,:,3) < sigThresh, 2 ) );
            boutSumm(x).(allVars{v}).fracBout = numel(boutSumm(x).(allVars{v}).boutSig)/boutSumm(x).totBout;
            % How many planes/ROI were consistently affected by locomotion bouts?
            boutSumm(x).(allVars{v}).sigUnit = find( median(tempPvar(:,:,3), 1, 'omitnan') < sigThresh );
            boutSumm(x).(allVars{v}).fracUnit = numel(boutSumm(x).(allVars{v}).sigUnit)/size(tempPvar,2);
        else
            boutSumm(x).(allVars{v}).effect = [];
            boutSumm(x).(allVars{v}).effectSig = [];
            boutSumm(x).(allVars{v}).sigBout = [];
            boutSumm(x).(allVars{v}).fracBout = NaN;
            boutSumm(x).(allVars{v}).sigUnit = [];
            boutSumm(x).(allVars{v}).fracUnit = NaN;
        end
    end
    % Read off effects from locomotion onset
    onStruct(x).T = periBout{x}(1).on.T;
    onStruct(x).ind.base = find(onStruct(x).T < -1);
    onStruct(x).ind.short = find(onStruct(x).T > 0 & onStruct(x).T < 2);
    onStruct(x).ind.long = find(onStruct(x).T >= 2 ); % & onStruct(x).T < 4
    onTemp = [periBout{x}.on];
    onStruct(x).velocity = cat(2, onTemp.velocity );
    if ismember(x, xCSD)
        onStruct(x).csd.pre.bout = 1:sum( [periBout{x}(1:expt(x).csd-1).Nbout] );
        onStruct(x).csd.pre.Nbout = numel(onStruct(x).csd.pre.bout);
        onStruct(x).csd.post.bout = onStruct(x).csd.pre.bout(end)+1:boutSumm(x).totBout;
        onStruct(x).csd.post.Nbout = numel(onStruct(x).csd.post.bout);
    else
        onStruct(x).csd.pre.bout = 1:boutSumm(x).totBout;
        onStruct(x).csd.pre.Nbout = boutSumm(x).totBout;
        onStruct(x).csd.post.bout = [];
        onStruct(x).csd.post.Nbout = 0;
    end
    for v = 1:NallVars
        if strcmpi(allVars{v}, 'speed')
            onStruct(x).(allVars{v}).data = permute( cat(2, onTemp.(allVars{v})), [1,3,2]); % data format: time x units x bouts 
        else
            onStruct(x).(allVars{v}).data = cat(3, onTemp.(allVars{v}));
        end
        onStruct(x).(allVars{v}).base = mean(onStruct(x).(allVars{v}).data(onStruct(x).ind.base,:,:), 1,'omitnan');
        onStruct(x).(allVars{v}).short = mean(onStruct(x).(allVars{v}).data(onStruct(x).ind.short,:,:), 1,'omitnan');
        onStruct(x).(allVars{v}).long = mean(onStruct(x).(allVars{v}).data(onStruct(x).ind.long,:,:), 1,'omitnan');
        onStruct(x).(allVars{v}).shortEffect = permute(onStruct(x).(allVars{v}).short - onStruct(x).(allVars{v}).base, [3,2,1]); 
        onStruct(x).(allVars{v}).longEffect = permute(onStruct(x).(allVars{v}).long - onStruct(x).(allVars{v}).base, [3,2,1]);
        if strcmpi(allVars{v}, 'fluor') 
            [~, onStruct(x).(allVars{v}).pShort] = ttest( onStruct(x).(allVars{v}).shortEffect );
            [~, onStruct(x).(allVars{v}).pLong] = ttest( onStruct(x).(allVars{v}).longEffect );
            onStruct(x).(allVars{v}).responder = find(onStruct(x).(allVars{v}).pShort < sigThresh | onStruct(x).(allVars{v}).pLong < sigThresh );
        end
        
        if ismember(x, xCSD)
            % Get data and effects from CSD itself
            onStruct(x).csd.(allVars{v}).data = csdBout{x}(expt(x).csd).on.(allVars{v});
            onStruct(x).csd.(allVars{v}).base = mean(onStruct(x).csd.(allVars{v}).data(onStruct(x).ind.base,:,:), 1,'omitnan');
            onStruct(x).csd.(allVars{v}).short = mean(onStruct(x).csd.(allVars{v}).data(onStruct(x).ind.short,:,:), 1,'omitnan');
            onStruct(x).csd.(allVars{v}).long = mean(onStruct(x).csd.(allVars{v}).data(onStruct(x).ind.long,:,:), 1,'omitnan');
            onStruct(x).csd.(allVars{v}).shortEffect = permute(onStruct(x).(allVars{v}).short - onStruct(x).(allVars{v}).base, [3,2,1]); 
            onStruct(x).csd.(allVars{v}).longEffect = permute(onStruct(x).(allVars{v}).long - onStruct(x).(allVars{v}).base, [3,2,1]);
            % Split up bout onset effects between pre and post CSD, averaging over planes/ROI
            % Pre/Post-CSD short-latency effects
            onStruct(x).csd.pre.(allVars{v}).short = mean( onStruct(x).(allVars{v}).shortEffect(onStruct(x).csd.pre.bout,:), 2, 'omitnan');
            onStruct(x).csd.post.(allVars{v}).short = mean( onStruct(x).(allVars{v}).shortEffect(onStruct(x).csd.post.bout,:), 2, 'omitnan');
            [~, onStruct(x).csd.post.(allVars{v}).pShort] =  ttest2(onStruct(x).csd.pre.(allVars{v}).short, onStruct(x).csd.post.(allVars{v}).short);
            % Pre/Post-CSD long-latency effects
            onStruct(x).csd.pre.(allVars{v}).long = mean( onStruct(x).(allVars{v}).longEffect(onStruct(x).csd.pre.bout,:), 2, 'omitnan');
            onStruct(x).csd.post.(allVars{v}).long = mean( onStruct(x).(allVars{v}).longEffect(onStruct(x).csd.post.bout,:), 2, 'omitnan');
            [~, onStruct(x).csd.post.(allVars{v}).pLong] =  ttest2(onStruct(x).csd.pre.(allVars{v}).long, onStruct(x).csd.post.(allVars{v}).long);
        end

        effectPool.(allVars{v}).short = [effectPool.(allVars{v}).short; onStruct(x).(allVars{v}).shortEffect(:)];
        effectPool.(allVars{v}).long = [effectPool.(allVars{v}).long; onStruct(x).(allVars{v}).longEffect(:)];
    end
    % Identify ROI driven by locomotion
    [~, pPre] = ttest( onStruct(x).fluor.shortEffect(onStruct(x).csd.pre.bout,:), 0, 'tail',tailVal );
    rLocoPre{x} = find(pPre < sigThresh);
    % Are there enough loco bouts pre and post to be useful?
    if onStruct(x).csd.pre.Nbout >= minBout && onStruct(x).csd.post.Nbout >= minBout
        xCSDbout = [xCSDbout, x]; 
        [~, pPost] = ttest( onStruct(x).fluor.shortEffect(onStruct(x).csd.post.bout,:), 0, 'tail',tailVal );
        rLocoPost{x} = find(pPost < sigThresh);
    end
end