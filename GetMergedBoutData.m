function periBout = GetMergedBoutData( periParam, T, boutScans, filtState, loco, varargin ) % , isCSD , deformVars deform, fluor
checkOptionalInputs = @(x)(isempty(x) || isstruct(x) || isnumeric(x));
IP = inputParser;
addRequired( IP, 'periParam', @isstruct )
addRequired( IP, 'T', @isnumeric )
addRequired( IP, 'boutScans', @iscell )
addRequired( IP, 'filtState', @islogical )
addRequired( IP, 'loco', @isstruct )
addOptional( IP, 'deform', [], checkOptionalInputs )
addOptional( IP, 'fluor', [], checkOptionalInputs )
parse( IP, periParam, T, boutScans, filtState, loco, varargin{:} ); %defVars, deform, fluor,
deform = IP.Results.deform;
fluor = IP.Results.fluor;
Nbout = numel(boutScans);
%Nscan = numel(T);
Trun = T-T(1); % timing, within the run
dT = Trun(2);
[~,zero_scan] = min(abs(Trun - periParam.base));
Ton = Trun - Trun(zero_scan);
% Initialize the periBout structure
periBout = struct('Nbout',Nbout, 'scan',[], 'Nscan',nan(1,Nbout), 'dur',nan(1,Nbout), 'iso',[], 'T',[], 'Tstart',nan(1,Nbout), 'Tstop',nan(1,Nbout), 'preScan',[], 'boutScan',[], 'postScan',[], ...
    'velocity',[], 'speed',[], 'state',[], 'fluor',[], 'on',[], 'off',[] );  %
periBout.on.T = Ton(Ton <= periParam.on);
Nscan_Ton = numel(periBout.on.T);
periBout.on.velocity = nan(Nscan_Ton, Nbout); periBout.on.speed = nan(Nscan_Ton, Nbout);
if ~isempty(deform)
    deform_vars = fieldnames(deform);
    Ndeform_vars = numel(deform_vars);
    for v = 1:Ndeform_vars
        periBout.(deform_vars{v}) = [];
        periBout.on.(deform_vars{v}) = nan(Nscan_Ton, Nbout, size(deform.(deform_vars{v}), 2)); %periBout.(deformVars{v}){bout}(periBout.on.ind,:);
    end
end
if ~isempty(fluor)
    Nroi = size(fluor,2);
    periBout.on.fluor = nan(Nscan_Ton, Nbout, Nroi);
end

for bout = flip(1:Nbout)
    % Get full bout data
    periBout.scan{bout} = boutScans{bout}(1):boutScans{bout}(end);  % -1 -periParam.NbaseScan  +periParam.NbaseScan
    %periBout.scan{bout}(periBout.scan{bout} < 1 | periBout.scan{bout} > Nscan) = [];
    periBout.Nscan(bout) = numel(periBout.scan{bout});
    periBout.T{bout} = T( periBout.scan{bout} ); %

    % within the merged bout, when did the mouse actually start/stop running?
    firstRunFrame = find(filtState(boutScans{bout}), 1);
    lastRunFrame = find(filtState(boutScans{bout}), 1, 'last');
    periBout.Tstart(bout) = periBout.T{bout}(firstRunFrame); % T(boutScans{bout}(1))
    periBout.Tstop(bout) = periBout.T{bout}(lastRunFrame); % T(boutScans{bout}(end));
    periBout.dur(bout) = periBout.Tstop(bout) - periBout.Tstart(bout)+dT;
    % Break up the bout into pre/bout/post scan sets
    periBout.preScan{bout} = 1:firstRunFrame-1; %find( periBout.T{bout} < periBout.Tstart(bout) ); % scans within bout
    periBout.boutScan{bout} = firstRunFrame:lastRunFrame; %find(periBout.T{bout} >= periBout.Tstart(bout) & periBout.T{bout} <= periBout.Tstop(bout)); % scans within bout
    periBout.postScan{bout} = lastRunFrame+1:periBout.Nscan(bout); %find( periBout.T{bout} > periBout.Tstop(bout) ); % scans within bout

    % Grab peri-bout locomotion, fluor, and deformation data
    periBout.velocity{bout} = loco.Vdown(periBout.scan{bout});
    periBout.speed{bout} = loco.speedDown(periBout.scan{bout});
    periBout.state{bout} = filtState(periBout.scan{bout});
    if ~isempty(deform)
        for v = 1:Ndeform_vars
            periBout.(deform_vars{v}){bout} = deform.(deform_vars{v})(periBout.scan{bout},:);
            %periBout.effect.(deform_vars{v}){bout} = periBout.(deform_vars{v}){bout}
        end
    end
    if ~isempty(fluor), periBout.fluor{bout} = fluor(periBout.scan{bout},:);  end

    % Concatenate onset data into arrays (optional)
    if periParam.on > 0
        Ton_bout = periBout.T{bout} - periBout.Tstart(bout); % time relative to onset of bout
        on_bout_logical = Ton_bout <= periParam.on; % indices within the current bout, allow for the possibiltiy the bout might be cutoff at the beginning
        on_bout_scans = find(on_bout_logical) + (Nscan_Ton-sum(on_bout_logical));  % indices within the standardized bout.on arrays
        %on_bout_scans(on_bout_scans < 1) = [];
        if any(periBout.velocity{bout}(on_bout_logical) >= periParam.min_vel_on)
            periBout.on.velocity(on_bout_scans,bout) =  periBout.velocity{bout}(on_bout_logical); % periBout.on.ind
            periBout.on.speed(on_bout_scans,bout) = periBout.speed{bout}(on_bout_logical);
            if ~isempty(fluor), periBout.on.fluor(on_bout_scans,bout,:) = periBout.fluor{bout}(on_bout_logical,:); end
            if ~isempty(deform)
                for v = 1:Ndeform_vars
                    periBout.on.(deform_vars{v})(on_bout_scans,bout,:) = periBout.(deform_vars{v}){bout}(on_bout_logical,:);
                end
            end
        end
    end
end
end