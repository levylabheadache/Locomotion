function periBout = GetBoutData( boutScans, T, loco, deform, fluor, periParam ) % , isCSD , deformVars

Nbout = numel(boutScans);
periBout = struct('Nbout',Nbout, 'dur',nan(1,Nbout), 'scan',[], 'Nscan',nan(1,Nbout), 'iso',[], 'T',[], 'Tstart',nan(1,Nbout), 'Tstop',nan(1,Nbout), 'preScan',[], 'boutScan',[], 'postScan',[], ...
    'velocity',[], 'speed',[], 'fluor',[], 'on',[], 'off',[] );  % 
deform_vars = fieldnames(deform);
Ndeform_vars = numel(deform_vars);
Nscan = numel(T);
Nroi = size(fluor,2);
Trun = T-T(1);
[~,zero_scan] = min(abs(Trun - periParam.base));
Ton = Trun - Trun(zero_scan);
periBout.on.T = Ton(Ton <= periParam.on);
Nscan_Ton = numel(periBout.on.T);
periBout.on.velocity = nan(Nscan_Ton, Nbout); periBout.on.speed = nan(Nscan_Ton, Nbout); periBout.on.fluor = nan(Nscan_Ton, Nbout, Nroi);
for v = 1:Ndeform_vars
    periBout.(deform_vars{v}) = [];
    periBout.on.(deform_vars{v}) = nan(Nscan_Ton, Nbout, size(deform.(deform_vars{v}), 2)); %periBout.(deformVars{v}){bout}(periBout.on.ind,:);
end

for bout = flip(1:Nbout)
    periBout.dur(bout) = T(boutScans{bout}(end)) - T(boutScans{bout}(1));
    % Get full bout data
    periBout.scan{bout} = boutScans{bout}(1)-periParam.NbaseScan:boutScans{bout}(end)+periParam.NbaseScan;  % -1
    periBout.scan{bout}(periBout.scan{bout} < 1 | periBout.scan{bout} > Nscan ) = [];
    periBout.Nscan(bout) = numel(periBout.scan{bout});
    periBout.T{bout} = T( periBout.scan{bout} ); %
    periBout.Tstart(bout) = T(boutScans{bout}(1));
    periBout.Tstop(bout) = T(boutScans{bout}(end));
    periBout.preScan{bout} = find( periBout.T{bout} < periBout.Tstart(bout) ); % scans within bout
    periBout.boutScan{bout} = find(periBout.T{bout} >= periBout.Tstart(bout) & periBout.T{bout} <= periBout.Tstop(bout)); % scans within bout
    periBout.postScan{bout} = find( periBout.T{bout} > periBout.Tstop(bout) ); % scans within bout
    
    % Grab peri-bout locomotion, fluor, and deformation data
    periBout.velocity{bout} = loco.Vdown(periBout.scan{bout});
    periBout.speed{bout} = loco.speedDown(periBout.scan{bout});
    for v = 1:Ndeform_vars
        periBout.(deform_vars{v}){bout} = deform.(deform_vars{v})(periBout.scan{bout},:);
    end
    periBout.fluor{bout} = fluor(periBout.scan{bout},:);
    
    % Concatenate onset data into arrays (optional)
    if periParam.on > 0
        Ton_bout = periBout.T{bout} - periBout.Tstart(bout); % time relative to onset of bout
        on_bout_logical = Ton_bout <= periParam.on; % indices within the current bout, allow for the possibiltiy the bout might be cutoff at the beginning
        on_bout_scans = find(on_bout_logical) + (Nscan_Ton-sum(on_bout_logical));  % indices within the standardized bout.on arrays
        %on_bout_scans(on_bout_scans < 1) = [];
        if any(periBout.velocity{bout}(on_bout_logical) >= periParam.min_vel_on)
            periBout.on.velocity(on_bout_scans,bout) =  periBout.velocity{bout}(on_bout_logical); % periBout.on.ind
            periBout.on.speed(on_bout_scans,bout) = periBout.speed{bout}(on_bout_logical);
            periBout.on.fluor(on_bout_scans,bout,:) = periBout.fluor{bout}(on_bout_logical,:);
            for v = 1:Ndeform_vars
                periBout.on.(deform_vars{v})(on_bout_scans,bout,:) = periBout.(deform_vars{v}){bout}(on_bout_logical,:);
            end
        end
    end
    % Concatenate offset data into arrays
end
end