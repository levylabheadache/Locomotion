function [stillEpoch, stillBin, stillSumm] = BinStillEpochs(expt, T, loco, fluor, deform, csdBout, varargin)
% Parse inputs
IP = inputParser;
checkData = @(x)(isempty(x) || isstruct(x));
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'T', @iscell )
addRequired( IP, 'loco', @isstruct )
addRequired( IP, 'fluor', checkData )
addRequired( IP, 'deform', checkData )
addRequired( IP, 'csdBout', checkData )
addParameter( IP, 'criterion', 'bout', @ischar )
addParameter( IP, 'minStillDur', 90, @isnumeric ) % initial epochs must be at least this long
addParameter( IP, 'trim', 15, @isnumeric ) % trim the first/last seconds of each initial epoch
addParameter( IP, 'bin', 60, @isnumeric ) % divide each epoch into bin seconds intervals
addParameter( IP, 'maxMissing', 0.1, @isnumeric ) % exclude epochs where the median ROI is missing at least this fraction of data
addParameter( IP, 'pre', -3, @isnumeric ) % epochs/bins that started at least this many minutes before CSD bout are considered "pre-CSD" 
addParameter( IP, 'acute', 15, @isnumeric )  % epochs/bins that started after pre, and less than this many minutes after CSD bout are considered "acute" 
%addParameter( IP, 'post', 60, @isnumeric ) % post-CSD phase ends post minutes after beginning of CSD run
addParameter( IP, 'defVars', ["transMag","scaleMag","shearMag","shiftZ"], @isstring);
addParameter( IP, 'show', false, @islogical )
parse( IP, expt, T, loco, fluor, deform, csdBout, varargin{:} );
criterion = IP.Results.criterion;
minStillDur = IP.Results.minStillDur;
trim = IP.Results.trim;
stillChunkDur = IP.Results.bin;
if stillChunkDur > minStillDur
    warning('Chunk size > min still dur')
    stillChunkDur = minStillDur - 2*trim; 
end
if stillChunkDur <= 0, error('Trim size not compatible with bin size'); end
maxMissingFrac = IP.Results.maxMissing; %0.1;
preLimit = IP.Results.pre;
acuteLimit = IP.Results.acute;
%postLimit = IP.Results.post;
defVars = IP.Results.defVars;
NdefVars = numel(defVars);
show = IP.Results.show;
ROItype = 'ROI'; 
if ~isempty(fluor) && isfield(fluor(1).F, 'axon'), ROItype = 'axon';  end
trimScan = round(trim*expt.scanRate);
stillChunkNscan = round(stillChunkDur*expt.scanRate);
stillEpoch = repmat(struct('scan_full',[], 'Nscan_full',NaN, 'T_full',[], 'Tstart_full',NaN, 'Tmid_full',NaN, 'Tstop_full',NaN, 'dur_full',[], 'dur_trim',[], 'raster_full',[], 'Froi',[], 'Fnp',[], 'Fo',[], 'dFF',[], 'scale',[],...
    'Nmissing',[], 'missingFrac',[], 'missingFracMed',[], 'roiDur',[], 'roiDurMed',[] ), 1, expt.Nruns); % 'run',NaN, 
stillBin = repmat(struct('run',NaN, 'epoch',NaN, 'scan',[], 'Nscan',NaN, 'T',[], 'Tstart',NaN, 'Tmid',NaN, 'Tstop',NaN, 'dur',NaN, 'raster',[], 'Froi',[], 'Fnp',[], 'Fo',[], 'dFF',[], 'scale',[],...
    'Nmissing',[], 'missingFrac',[], 'missingFracMed',[], 'roiDur',[], 'roiDurMed',[] ), 0, 1);
stillSumm = struct('Nepoch',NaN, 'epoch_run',[], 'dur_full',[], 'dur_full_tot',[], 'dur_trim',[], 'dur_trim_tot',[]);
S = 0; E = 0;
for runs = 1:expt.Nruns
    if show
        stateVec = nan(length(loco(runs).stateDown), 1); %loco(runs).stateDown > 1;
        stateVec(loco(runs).stateDown == 2) = 1;
        figure('WindowState','maximized');
        plot( loco(runs).Vdown); hold on; % T{runs}/60,
        plot(-1.5*stateVec, 'color','r', 'LineWidth',2 );
        plot(-1*loco(runs).bout, 'color','b', 'LineWidth',2 ); % T{runs}/60,
        xlim([-Inf,Inf])
    end
    % Find sufficiently long periods of stillness (epochs)
    if strcmpi(criterion, 'bout')
        tempStillCC = bwconncomp(loco(runs).bout ~= 1); %bwconncomp(loco(run).stateDown == 1);
    elseif strcmpi(criterion, 'state')
        tempStillCC = bwconncomp(loco(runs).stateDown == 1);
    elseif strcmpi(criterion, 'speed')
        tempStillCC = bwconncomp(loco(runs).speedDown == 0);
    else
        error('Criterion must be bout, state, or speed')
    end
    tempStillDur = cellfun(@numel, tempStillCC.PixelIdxList)/expt.scanRate;
    runEpochScans = tempStillCC.PixelIdxList(tempStillDur >= minStillDur);
    stillEpoch(runs).Nepoch = numel(runEpochScans);
    % Gather information on each epoch and perform subepoch binning
    for ep = 1:stillEpoch(runs).Nepoch
        E = E+1;
        stillSumm.epoch_run(E) = runs;
        % Full-epoch data
        %stillEpoch(ep).run = runs;
        stillEpoch(runs).scan_full{ep} = runEpochScans{ep}';
        stillEpoch(runs).Nscan_full(ep) = numel(stillEpoch(runs).scan_full{ep});
        stillEpoch(runs).T_full{ep} = T{runs}(stillEpoch(runs).scan_full{ep}); % -Tscan{runs}(stillEpoch(runs).scan(1))
        stillEpoch(runs).Tstart_full(ep) = stillEpoch(runs).T_full{ep}(1);
        stillEpoch(runs).Tmid_full(ep) = median(stillEpoch(runs).T_full{ep});
        stillEpoch(runs).Tstop_full(ep) = stillEpoch(runs).T_full{ep}(end);
        stillEpoch(runs).dur_full(ep) = stillEpoch(runs).Tstop_full(ep) - stillEpoch(runs).Tstart_full(ep);
        % Trim epoch to mitigate post/pre locomotion contamination
        stillEpoch(runs).scan_trim{ep} = runEpochScans{ep}(trimScan):runEpochScans{ep}(end-trimScan); %runEpochScans{ep}([trimScan,end-trimScan]);
        stillEpoch(runs).Nscan_trim(ep) = numel(stillEpoch(runs).scan_trim{ep});
        stillEpoch(runs).T_trim{ep} = T{runs}(stillEpoch(runs).scan_trim{ep}); % -Tscan{runs}(stillEpoch(runs).scan(1))
        stillEpoch(runs).Tstart_trim(ep) = stillEpoch(runs).T_trim{ep}(1);
        stillEpoch(runs).Tmid_trim(ep) = median(stillEpoch(runs).T_trim{ep});
        stillEpoch(runs).Tstop_trim(ep) = stillEpoch(runs).T_trim{ep}(end);
        stillEpoch(runs).dur_trim(ep) = stillEpoch(runs).Tstop_trim(ep) - stillEpoch(runs).Tstart_trim(ep);
        % When did this (trimmed) epoch start, relative to CSD bout
        if ~isempty(csdBout)
            stillEpoch(runs).Tcsd(ep) = (stillEpoch(runs).Tstart_trim(ep) - csdBout.Tstart)/60; % minutes
            if stillEpoch(runs).Tcsd(ep) < preLimit
                stillEpoch(runs).csd_phase{ep} = 'pre';
            elseif stillEpoch(runs).Tcsd(ep) >= preLimit && stillEpoch(runs).Tcsd(ep) < acuteLimit
                stillEpoch(runs).csd_phase{ep} = 'acute';
            else
                stillEpoch(runs).csd_phase{ep} = 'post';
            end
        else
            stillEpoch(runs).Tcsd(ep) = NaN;
            stillEpoch(runs).csd_phase{ep} = 'n/a';
        end
        
        % Gather fluor and deformation data
        if ~isempty(deform)
            for d = 1:NdefVars
                stillEpoch(runs).(defVars(d)){ep} = deform(runs).(defVars(d))(stillEpoch(runs).scan_trim{ep},:);
            end
        end

        % Break the (trimmed) epoch into bins of predetermined size
        [tempChunkLims, ~, tempChunkSize] = MakeChunkLims( stillEpoch(runs).scan_trim{ep}(1), stillEpoch(runs).scan_trim{ep}(end), stillEpoch(runs).scan_trim{ep}(end), 'size',stillChunkNscan);
        if tempChunkSize(end) < stillChunkNscan  % Absorb remainder into last bin
            tempChunkLims(end-1,2) = tempChunkLims(end,2);
            tempChunkLims(end,:) = [];
            tempChunkSize = diff(tempChunkLims,1,2)+1;
        end
        
        % Record bin-level data
        bin = 0;
        for chunk = find(tempChunkSize >= stillChunkNscan)'
            S = S+1; bin = bin+1;
            stillBin(S).run = runs;
            stillBin(S).epoch = ep;
            stillBin(S).bin = bin; % bin index within epoch
            stillBin(S).scan = tempChunkLims(chunk,1):tempChunkLims(chunk,2);
            stillBin(S).Nscan = numel(stillBin(S).scan);
            stillBin(S).T = T{runs}(stillBin(S).scan); % -Tscan{runs}(stillBin(S).scan(1))
            stillBin(S).Tstart = stillBin(S).T(1);
            stillBin(S).Tmid = median(stillBin(S).T);
            stillBin(S).Tstop = stillBin(S).T(end);
            stillBin(S).dur = stillBin(S).T(end) - stillBin(S).T(1);
            if ~isempty(fluor)
                stillBin(S).raster = false(stillBin(S).Nscan, expt.Nroi);
                % get fluor/deformation data during each epoch
                stillBin(S).Froi = median( fluor(runs).Froi.(ROItype)(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).Fnp = median( fluor(runs).Fnp.(ROItype)(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).Fo = median( fluor(runs).Fo.(ROItype)(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).dFF = fluor(runs).dFF.(ROItype)(stillBin(S).scan,:);
                % account for missing (poorly registered or otherwise excluded) (NaN) data
                stillBin(S).Nmissing = sum(isnan( stillBin(S).dFF), 1);
                stillBin(S).missingFrac = stillBin(S).Nmissing/stillBin(S).Nscan;
                stillBin(S).missingFracMed = median(stillBin(S).missingFrac);
                stillBin(S).roiDur = stillBin(S).dur*(1-stillBin(S).missingFrac);
                stillBin(S).roiDurMed = median(stillBin(S).roiDur);
            else
                stillBin(S).raster = [];
                stillBin(S).Froi = [];
                stillBin(S).Fnp = [];
                stillBin(S).Fo = [];
                stillBin(S).dFF = [];
                stillBin(S).Nmissing = [];
                stillBin(S).missingFrac = [];
                stillBin(S).missingFracMed = [];
                stillBin(S).roiDur = [];
                stillBin(S).roiDurMed = [];
            end
            if ~isempty(deform)
                if expt.Nplane > 1
                    stillBin(S).scale = mean(deform(runs).scaleMag(stillBin(S).scan,4:end-3), 2, 'omitnan');
                else
                    stillBin(S).scale = deform(runs).scaleMag(stillBin(S).scan,:);
                end
            else
                stillBin(S).scale = [];
            end

            if show
                line(stillBin(S).scan([1,end]), -2*[1,1], 'color','k', 'lineWidth',2 );
                %pause;
                %plot(-2*vertcat(loco(runs).bout), 'color',[0.6,0.6,0.6], 'LineWidth',2 ); % T{runs}/60,
            end
        end
        stillEpoch(runs).Nchunk(ep) = chunk;
    end
end
if ~isempty(stillBin)
    stillBin([stillBin.missingFracMed] > maxMissingFrac) = []; % Remove any epochs where fluor data is mostly missing
end
%{
if show
    figure;
    for runs = 1:expt.Nruns
        plot(T{runs}/60, loco(runs).Vdown); hold on;
        plot(T{runs}/60, -1*vertcat(loco(runs).bout), 'color',[0.6,0.6,0.6], 'LineWidth',2 );
        for s = 1:S
            line( T{runs}(stillBin(s).scan([1,end]))'/60, -2*[1,1], 'color','k', 'LineWidth',2 )
        end
%{
        for s = 1:S
            plot(T{runs}(stillBin(s).scan(1))/60, loco(runs).Vdown(stillBin(s).scan(1)), 'ko');
            plot(T{runs}(stillBin(s).scan(end))/60, loco(runs).Vdown(stillBin(s).scan(end)), 'kx');
        end
%}
        %pause;
    end
end
%}


% Summarize epoch-level data
stillSumm.Nepoch = sum([stillEpoch.Nepoch]);
stillSumm.dur_full = [stillEpoch.dur_full];
stillSumm.dur_full_tot = sum(stillSumm.dur_full);
stillSumm.dur_trim = [stillEpoch.dur_trim];
stillSumm.dur_trim_tot = sum(stillSumm.dur_trim);


end