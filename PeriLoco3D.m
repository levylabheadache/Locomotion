function [periBout, periParam, loco] = PeriLoco3D( expt, T, loco, deform, fluor, defVars, varargin ) % , periStat
% PeriBoutData gather relevant data from a defined period of time around the onset of each bout of running. 
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'loco', @isstruct )
addRequired( IP, 'deform', @isstruct )
addRequired( IP, 'fluor', @isnumeric )
addRequired( IP, 'defVars', @iscell )
addParameter( IP, 'iso', [], @isnumeric ) % bout must be separated by at least this much time [before, after] from the next bout
addParameter( IP, 'base', 15, @isnumeric )
addParameter( IP, 'run', 5, @isnumeric ) % bout must last at least this long
addParameter( IP, 'on', 5, @isnumeric )
addParameter( IP, 'min_vel_on', 0, @isnumeric)
addParameter( IP, 'Nshuff', 0, @isnumeric )
addParameter( IP, 'sig', 0.05, @isnumeric )
parse( IP, expt, T, loco, deform, fluor, defVars, varargin{:} ); 
periParam.base = IP.Results.base;
periParam.run = IP.Results.run;
periParam.on = IP.Results.on;
periParam.iso = IP.Results.iso;
periParam.min_vel_on = IP.Results.min_vel_on;
if isempty(periParam.iso), periParam.iso = periParam.base*[1,0]; end
if numel(periParam.iso) == 1, periParam.iso = periParam.iso*[1,1]; end
periParam.Nshuff = IP.Results.Nshuff;
periParam.sigThresh = IP.Results.sig;
if periParam.on > periParam.run
    warning('Onset time exceeds minimum run time. Setting on to run');
    periParam.on = periParam.run;
end
%periParam.base = 3; periParam.run = 2; binDur = 0.5; periParam.iso = 2*periParam.base;
[Nscan, Nstate] = size( loco.stateBinary );
NdefVars = numel(defVars);
periParam.dT = 1/expt.scanRate;
periParam.NbaseScan = round(expt.scanRate*periParam.base); % ceil
periParam.NrunScan = round(expt.scanRate*periParam.run); % 
periParam.T = periParam.dT*(-periParam.NbaseScan:periParam.NrunScan); %linspace( -periParam.base, periParam.run, NsubScan); 
%{
Tfull = periParam.dT*(0:NmergedScan-1)';
plot( T, loco.state ); box off; ylim([-0.01, 1.01]); xlim([0,Inf]);
%}
runState = loco.stateBinary(:,Nstate); %logical(loco.locoState(:,Nstate))'; % plot( runState ); ylim([0,1.3]);
runBout = regionprops( runState, 'Area', 'PixelIdxList' );
Nputative = numel(runBout);
putativeBout = repmat( struct('dur',NaN, 'iso',nan(1,2), 'scan',[], 'Nscan',NaN), 1, Nputative );
for p = 1:numel(runBout)
    % how isolated is the run bout from neighboring bouts?
    preIsoScan = find( runState(1:runBout(p).PixelIdxList(1)-1), 1, 'last' ); % where was the last scan from the previous putatative bout?
    if isempty(preIsoScan)
        preIsoScan = 0; % 1
    end
    postIsoScan = find( runState(runBout(p).PixelIdxList(end)+1:end), 1, 'first' ); % where's the first scan of the next putatative bout?
    if ~isempty(postIsoScan) 
        postIsoScan = postIsoScan + runBout(p).PixelIdxList(end); 
    else 
        postIsoScan = Nscan+1;  % Nscan
    end
    tempIsoScan = [runBout(p).PixelIdxList(1)-(preIsoScan+1),  postIsoScan-(runBout(p).PixelIdxList(end)+1)]; % add ones since isolation is number of still scans between putative bouts
    % get its basic info
    putativeBout(p).dur = periParam.dT*runBout(p).Area;
    putativeBout(p).iso = periParam.dT*tempIsoScan;
    putativeBout(p).scan = runBout(p).PixelIdxList';
    putativeBout(p).Nscan = runBout(p).Area;
end
%{
figure('Units','normalized','OuterPosition',[0.16,0.13,0.7,0.8], 'Color','w');
subplot(1,2,1); ecdf( [putativeBout.dur] ); hold on; line(periParam.run*[1,1],[0,1],'color','r','lineStyle','--'); xlabel('Putative Duration (s)'); ylabel('CDF'); axis square;
subplot(1,2,2); ecdf( [putativeBout.iso] ); hold on; line(periParam.iso,[0,1],'color','r','lineStyle','--'); xlabel('Isolation (s)'); set(gca,'Xscale','log'); axis square; % ylabel('CDF'); 
%}
putIso = reshape( [putativeBout.iso], 2, Nputative )';
pGood = find( [putativeBout.dur]' >= periParam.run & all([putIso(:,1) > periParam.iso(1),putIso(:,2) > periParam.iso(2)],2) )';  % enforce minimum bout duration, isolation

%pGood = find( [putativeBout.dur]' >= periParam.run & all(reshape( [putativeBout.iso], 2, Nputative )' >= periParam.iso, 2) )';  % enforce minimum bout duration, isolation
boutScans = {putativeBout(pGood).scan};
% {
figure('Units','normalized','OuterPosition',[0.16,0.13,0.7,0.8], 'Color','w');
sp(1) = subplot(3,1,1); plot( loco.speedDown ); xlim([-Inf,Inf]); %ylim([-0.01, 1.01]);
sp(2) = subplot(3,1,2); plot( runState ); xlim([-Inf,Inf]); ylim([-0.01, 1.01]);
goodBouts = zeros(size(runState) ); goodBouts( [putativeBout(pGood).scan] ) = 1;
sp(3) = subplot(3,1,3); plot(goodBouts); xlim([-Inf,Inf]); ylim([-0.01, 1.01]);
linkaxes(sp,'xy');
%}

% Construct the final bouts
periBout = GetBoutData( boutScans, T, loco, deform, fluor, defVars, periParam );
periBout.iso = vertcat(putativeBout(pGood).iso);

%periBout.dur = [putativeBout(pGood).dur];


% Summary statistics for each bout
subStruct = struct('pre',[], 'run',[], 'post',[], 'effect',[], 'effectCI',[], 'effectSig',[]);
if periBout.Nbout > 0 
    loco.bout( [periBout.scan{:}] ) = 1;
    % Summarize effect on readout variables before/during/after bout
    periBout.stat.speed = BoutEffect( periBout.speed, periBout.preScan, periBout.boutScan, periBout.postScan ); % effect = [run - prerun, run-post, overall]
    periBout.stat.fluor = BoutEffect( periBout.fluor, periBout.preScan, periBout.boutScan, periBout.postScan );
    for v = 1:NdefVars
        periBout.stat.(defVars{v}) = BoutEffect( periBout.(defVars{v}), periBout.preScan, periBout.boutScan, periBout.postScan );
    end
    % Compare effect sizes to those obtained from shuffled data
    if periParam.Nshuff > 0 
        %{
        periBout.stat.speed.pEffect = BoutEffectPval(periBout, periBout.stat, periParam, 'speed');
        periBout.stat.fluor.pEffect = BoutEffectPval(periBout, periBout.stat, periParam, 'fluor');
        for v = 1:NdefVars
            periBout.stat.(defVars{v}).pEffect = BoutEffectPval(periBout, periBout.stat, periParam, defVars{v});
        end
        %}
        [periBout.stat.speed.effectCI, periBout.stat.speed.effectSig] = BoutEffectConfInt(periBout, periBout.stat, periParam, 'speed');
        [periBout.stat.fluor.effectCI, periBout.stat.fluor.effectSig] = BoutEffectConfInt(periBout, periBout.stat, periParam, 'fluor');
        %{
        for v = 1:NdefVars [periBout.stat.(defVars{v}).effectCI,
        periBout.stat.(defVars{v}).effectSig] = BoutEffectConfInt(periBout,
        periBout.stat, periParam, defVars{v}); end
        %}
    else
        periBout.stat.speed.effectSig = [];
        periBout.stat.fluor.effectSig = [];
        for v = 1:NdefVars
            periBout.stat.(defVars{v}).effectSig = [];
        end
    end
else
    periBout.stat.speed = subStruct;
    periBout.stat.fluor = subStruct;
    for v = 1:NdefVars
        periBout.stat.(defVars{v}) = subStruct;
    end
end
end