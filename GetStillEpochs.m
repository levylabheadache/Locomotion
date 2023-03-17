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
addParameter( IP, 'chunk', 60, @isnumeric ) % divide each epoch into chunk seconds intervals
addParameter( IP, 'maxMissing', 0.1, @isnumeric ) % exclude epochs where the median ROI is missing at least this fraction of data
addParameter( IP, 'acute', 15, @isnumeric ) % define acute phase as first minutes after beginning of CSD run
addParameter( IP, 'post', 60, @isnumeric ) % post-CSD phase ends post minutes after beginning of CSD run
addParameter( IP, 'show', false, @islogical )
parse( IP, expt, T, loco, fluor, deform, csdBout, varargin{:} );
criterion = IP.Results.criterion;
minStillDur = IP.Results.minStillDur;
trim = IP.Results.trim;
stillChunkDur = IP.Results.chunk;
maxMissingFrac = IP.Results.maxMissing; %0.1;
acuteLimit = IP.Results.acute;
postLimit = IP.Results.post;
show = IP.Results.show;
trimScan = round(trim*expt.scanRate);
stillChunkNscan = round(stillChunkDur*expt.scanRate);
stillEpoch = repmat(struct('scan_full',[], 'Nscan_full',NaN, 'T_full',[], 'Tstart_full',NaN, 'Tmid_full',NaN, 'Tstop_full',NaN, 'dur_full',NaN, 'raster_full',[], 'Froi',[], 'Fnp',[], 'Fo',[], 'dFF',[], 'scale',[],...
    'Nmissing',[], 'missingFrac',[], 'missingFracMed',[], 'roiDur',[], 'roiDurMed',[] ), 1, expt.Nruns); % 'run',NaN, 
stillBin = repmat(struct('run',NaN, 'scan',[], 'Nscan',NaN, 'T',[], 'Tstart',NaN, 'Tmid',NaN, 'Tstop',NaN, 'dur',NaN, 'raster',[], 'Froi',[], 'Fnp',[], 'Fo',[], 'dFF',[], 'scale',[],...
    'Nmissing',[], 'missingFrac',[], 'missingFracMed',[], 'roiDur',[], 'roiDurMed',[] ), 0, 1);
S = 0;
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
    else
        error('Criterion must be state or bout')
    end
    tempStillDur = cellfun(@numel, tempStillCC.PixelIdxList)/expt.scanRate;
    runEpochScans = tempStillCC.PixelIdxList(tempStillDur >= minStillDur);
    stillEpoch(runs).Nepoch = numel(runEpochScans);
    % Gather information on each epoch and perform subepoch binning
    for ep = 1:stillEpoch(runs).Nepoch
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
        % Break each (trimmed) epoch of stillness into bins of predetermined size
        [tempChunkLims, ~, tempChunkSize] = MakeChunkLims( stillEpoch(runs).scan_trim{ep}(1), stillEpoch(runs).scan_trim{ep}(end), stillEpoch(runs).scan_trim{ep}(end), 'size',stillChunkNscan);
        if tempChunkSize(end) < stillChunkNscan  % Absorb remainder into last chunk
            tempChunkLims(end-1,2) = tempChunkLims(end,2);
            tempChunkLims(end,:) = [];
            tempChunkSize = diff(tempChunkLims,1,2)+1;
        end
        bin = 0;
        % Record bin-level data
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
                stillBin(S).Froi = median( fluor(runs).Froi.ROI(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).Fnp = median( fluor(runs).Fnp.ROI(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).Fo = median( fluor(runs).Fo.ROI(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).dFF = fluor(runs).dFF.ROI(stillBin(S).scan,:);
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
        stillEpoch(runs).Nbin(ep) = bin;
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
stillSumm.dur = [stillBin.dur];
stillSumm.TotDur = sum(stillSumm.dur);
%{
% Divide up pre-CSD, acute and post-CSD phase epochs
if ~isnan(expt.csd)
    stillBin([stillBin.run] == expt.csd & [stillBin.Tstart] < csdBout(expt.csd).Tstop) = []; % Remove any epochs that include the CSD wave itself
    stillSumm.Trel = ([stillBin.Tmid]-csdBout(expt.csd).Tstart)/60; % epoch timing, relative to CSD, in minutes
else
    stillSumm.Trel = ([stillBin.Tmid]-T{end}(end))/60; %Trel = [stillBin.Tmid]/60; % epoch timing, relative to CSD, in minutes
end
stillSumm.sPre = find(stillSumm.Trel < 0); %stillSumm.sPre = find([stillBin.Tmid] < csdBout(expt.csd).Tstart);
stillSumm.sAcute = find(stillSumm.Trel >= 0 & stillSumm.Trel < acuteLimit);
stillSumm.sPost = find(stillSumm.Trel >= acuteLimit & stillSumm.Trel <= postLimit);
stillSumm.PreDur = sum(stillSumm.dur(stillSumm.sPre), 2);
stillSumm.AcuteDur = sum(stillSumm.dur(stillSumm.sAcute), 2);
stillSumm.PostDur = sum(stillSumm.dur(stillSumm.sPost), 2);
if ~isempty(fluor)
    stillSumm.roiDur = vertcat(stillBin.roiDur);
    stillSumm.roiTotDur = sum(stillSumm.roiDur, 1);
    stillSumm.roiPreDur = sum(stillSumm.roiDur(stillSumm.sPre,:), 1);
    stillSumm.roiAcuteDur = sum(stillSumm.roiDur(stillSumm.sAcute,:), 1);
    stillSumm.roiPostDur = sum(stillSumm.roiDur(stillSumm.sPost,:), 1);
else
    stillSumm.roiDur = [];
    stillSumm.roiTotDur = [];
    stillSumm.roiPreDur = [];
    stillSumm.roiAcuteDur = [];
    stillSumm.roiPostDur = [];
end
%}

    %stillFluorPoolZ = normalize(vertcat(stillBin(:).dFF));
    %stillScalePoolZ = normalize(vertcat(stillBin(:).scale));

    % Unpool the normalized data
    %{
for s = 1:stillSumm.Nepoch
    tempScan = 1:stillBin(s).Nscan;
    stillBin(s).z = stillFluorPoolZ(tempScan,:);
    stillBin(s).zScale = stillScalePoolZ(tempScan,:);
    stillFluorPoolZ(tempScan,:) = [];
    stillScalePoolZ(tempScan,:) = [];
end
    %}

    % Get distribution of deformation variables pre/acute/post-CSD
    %{
if ~isempty(stillSumm.sPre)
    tempPreDeform = vertcat( stillBin(stillSumm.sPre).scale ); % Pre
    [stillSumm.scale.CDFf{1}, stillSumm.scale.CDFx{1}] = ecdf(tempPreDeform);
end
if ~isempty(stillSumm.sAcute)
    tempAcuteDeform = vertcat( stillBin(stillSumm.sAcute).scale ); % Acute
    [stillSumm.scale.CDFf{2}, stillSumm.scale.CDFx{2}] = ecdf(tempAcuteDeform);
end
if ~isempty(stillSumm.sPost)
    tempPostDeform = vertcat( stillBin(stillSumm.sPost).scale ); % Post
    [stillSumm.scale.CDFf{3}, stillSumm.scale.CDFx{3}] = ecdf(tempPostDeform);
end    
    %}

        %{
    for s = find(tempStillDur >= minStillDur)
        % Break each epoch of stillness into bins of predetermined size
        tempStillScanRange = tempStillCC.PixelIdxList{s}([trimScan,end-trimScan]);
        [tempChunkLims, ~, tempChunkSize] = MakeChunkLims( tempStillScanRange(1), tempStillScanRange(end), tempStillScanRange(end), 'size',stillChunkNscan);
        % Absorb remainder into last chunk
        if tempChunkSize(end) < stillChunkNscan
            tempChunkLims(end-1,2) = tempChunkLims(end,2);
            tempChunkLims(end,:) = [];
            tempChunkSize = diff(tempChunkLims,1,2)+1;
        end
        for chunk = find(tempChunkSize >= stillChunkNscan)'
            S = S+1;
            stillBin(S).run = runs;
            stillBin(S).epoch = s;
            stillBin(S).scan = tempChunkLims(chunk,1):tempChunkLims(chunk,2);
            stillBin(S).Nscan = numel(stillBin(S).scan);
            stillBin(S).T = T{runs}(stillBin(S).scan); % -Tscan{run}(stillBin(S).scan(1))
            stillBin(S).Tstart = stillBin(S).T(1);
            stillBin(S).Tmid = median(stillBin(S).T);
            stillBin(S).Tstop = stillBin(S).T(end);
            stillBin(S).dur = stillBin(S).T(end) - stillBin(S).T(1);
            if ~isempty(fluor)
                stillBin(S).raster = false(stillBin(S).Nscan, expt.Nroi);
                % get fluor/deformation data during each epoch
                stillBin(S).Froi = median( fluor(runs).Froi.ROI(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).Fnp = median( fluor(runs).Fnp.ROI(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).Fo = median( fluor(runs).Fo.ROI(stillBin(S).scan,:), 1, 'omitnan' );
                stillBin(S).dFF = fluor(runs).dFF.ROI(stillBin(S).scan,:);
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
    end
    %}
end