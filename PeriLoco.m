function [periBout, periParam, loco] = PeriLoco( expt, T, loco, varargin ) % , periStat , defVars , deform, fluor
% PeriBoutData gather relevant data from a defined period of time around the onset of each bout of running.
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'loco', @isstruct )
addOptional( IP, 'deform', [], @isstruct )
addOptional( IP, 'fluor', [], @isnumeric )
addParameter( IP, 'merge', false, @islogical)
addParameter( IP, 'iso', [], @isnumeric ) % bout must be separated by at least this much time [before, after] from the next bout
addParameter( IP, 'base', 15, @isnumeric ) % get this much time prior to bout onset, to use as a baseline
addParameter( IP, 'run', 5, @isnumeric ) % bout must last at least this long
addParameter( IP, 'on', 5, @isnumeric )
addParameter( IP, 'min_vel_on', 0, @isnumeric)
addParameter( IP, 'Nshuff', 0, @isnumeric )
addParameter( IP, 'show', false, @islogical )
parse( IP, expt, T, loco, varargin{:} ); %defVars, deform, fluor, 
deform = IP.Results.deform;
fluor = IP.Results.fluor;
periParam.base = IP.Results.base;
periParam.run = IP.Results.run;
periParam.on = IP.Results.on;
periParam.merge = IP.Results.merge;
periParam.iso = IP.Results.iso;
periParam.min_vel_on = IP.Results.min_vel_on;
if isempty(periParam.iso), periParam.iso = periParam.base*[1,0]; end
if numel(periParam.iso) == 1, periParam.iso = periParam.iso*[1,1]; end
show = IP.Results.show;
if periParam.on > periParam.run
    warning('Onset time exceeds minimum run time. Setting on to run');
    periParam.on = periParam.run;
end
%periParam.base = 3; periParam.run = 2; binDur = 0.5; periParam.iso = 2*periParam.base;
[Nscan, Nstate] = size( loco.stateBinary );
%NdefVars = numel(defVars);
periParam.dT = 1/expt.scanRate;
periParam.NbaseScan = round(expt.scanRate*periParam.base); % ceil
periParam.NrunScan = round(expt.scanRate*periParam.run); %
periParam.T = periParam.dT*(-periParam.NbaseScan:periParam.NrunScan); %linspace( -periParam.base, periParam.run, NsubScan);
%{
Tfull = periParam.dT*(0:NmergedScan-1)';
plot( T, loco.state ); box off; ylim([-0.01, 1.01]); xlim([0,Inf]);
%}
runState = loco.stateBinary(:,Nstate); %logical(loco.locoState(:,Nstate))'; % plot( runState ); ylim([0,1.3]);
runBout = regionprops( runState, loco.Vdown, 'Area', 'PixelIdxList', 'MaxIntensity' );
filtState = false(size(runState));
filtState(vertcat(runBout([runBout.Area] >= periParam.NrunScan & [runBout.MaxIntensity] >= periParam.min_vel_on).PixelIdxList)) = true; % filter the runState to exclude weak/brief running
if periParam.merge
    mergeState = imdilate(filtState, strel("square",round(2*expt.scanRate*periParam.base))); % pad the filtState bouts 
else
    mergeState = filtState;
end
prelimBout = regionprops( mergeState, loco.Vdown, 'Area', 'PixelIdxList', 'MaxIntensity' );
boutScans = {prelimBout.PixelIdxList};

% Construct the final bouts
periBout = GetMergedBoutData( periParam, T, boutScans, filtState, loco, deform, fluor ); % boutScans, filtState, T, loco, periParam


periBout.iso = nan(periBout.Nbout,2);
periState = false(size(runState));
if periBout.Nbout > 0
    periState([periBout.scan{:}]) = true;
    for p = 1:periBout.Nbout
        % how isolated is each run bout from neighboring bouts?
        if p == 1
            periBout.iso(p,1) = periBout.Tstart(1)-T(1);
            if periBout.Nbout > 1
                periBout.iso(p,2) = periBout.Tstart(2) - periBout.Tstop(1);
            else
                periBout.iso(p,2) = T(end) - periBout.Tstop(1);
            end
        elseif p == periBout.Nbout
            periBout.iso(p,1) = periBout.Tstart(p) - periBout.Tstop(p-1);
            periBout.iso(p,2) = T(end) - periBout.Tstop(p);
        else
            periBout.iso(p,1) = periBout.Tstart(p) - periBout.Tstop(p-1);
            periBout.iso(p,2) = periBout.Tstart(p+1) - periBout.Tstop(p);
        end
    end
end

if show
    figure('Units','normalized','OuterPosition',[0.16,0.13,0.7,0.8], 'Color','w');
    sp(1) = subplot(2,1,1);
    plot( T, loco.speedDown );
    hold on;
    for p = 1:periBout.Nbout
        plot(periBout.Tstart(p), loco.speedDown(periBout.scan{p}(periBout.boutScan{p}(1))), 'o' );
        plot(periBout.Tstop(p), loco.speedDown(periBout.scan{p}(periBout.boutScan{p}(end))), 'x' );
        %line(periBout.scan{p}([1,end]), -1*[1,1], 'color','k','linewidth',2)
        %plot(periBout.scan)
    end
    xlim([-Inf,Inf]); %ylim([-0.01, 1.01]);

    sp(2) = subplot(2,1,2);
    plot( T, runState, 'r' ); hold on;
    plot( T, mergeState, 'b' );
    plot( T, periState, 'k' );
    legend(["run","merged","peri"])
    linkaxes(sp,'x');
    xlim([-Inf,Inf]); ylim([-0.01, 1.01]);
end
end

%{
    plot( filtState, 'r' ); hold on;
    plot( mergeState, 'b' );
    %plot( periState, 'k' );
xlim([-Inf,Inf]); ylim([-0.01, 1.01]);
%}