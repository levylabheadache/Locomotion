function loco = GetLocoData(infoStruct, varargin) % quadPath,
% Load wheel data and convert to filtered wheel velocity
IP = inputParser;
addRequired( IP, 'infoStruct', @isstruct )
addParameter( IP, 'wheel_diameter', 14, @isnumeric ) % wheel_diameter = 14; % in cm
addParameter( IP, 'wheel_tabs', 44, @isnumeric ) % wheel_tabs = 44;
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'show', false, @islogical )
parse( IP, infoStruct, varargin{:} );
wheel_diameter = IP.Results.wheel_diameter;
wheel_tabs = IP.Results.wheel_tabs;
wheel_circumference = wheel_diameter*pi;
step_size = wheel_circumference/(wheel_tabs*2);
overwrite = IP.Results.overwrite;
show = IP.Results.show;

% Check if the data is already saved in the expt level folder, and if so, load it
dirParts = strsplit(infoStruct.dir, '\');
exptDir = join(dirParts(1:end-2), '\');
exptDir = [exptDir{1},'\'];
[~,locoDataPath] = FileFinder(exptDir, 'type','mat', 'contains','locomotion');
if ~isempty(locoDataPath) && ~overwrite
    fprintf('\nLoading %s', locoDataPath{1})
    load(locoDataPath{1}, 'locoData');
    runStrInd = strfind(infoStruct.exptName, 'run'); %which run is this?
    runNum = str2double(infoStruct.exptName(runStrInd+3:end));
    loco = locoData(runNum);
else
    frameRate = 15.49; %infoStruct.framerate; %15.49; % in some cases (when FixSBX was used to resample data) effective framerate may be different
    % round(frameRate/infoStruct.framerate); % allow for data where many planes were acquired, but only one plane was sampled
    if infoStruct.optotune_used
        Nplane = infoStruct.otparam(3);
    else
        Nplane = 1;
    end
    % Load quad data and calculate velocity
    [~,quadPath] = FileFinder( infoStruct.dir, 'type','mat', 'contains', 'quad' );
    loco = struct('quad',[], 'Vinst',[], 'Ainst',[], 'Vfilt',[], 'Vdown',[], 'Adown',[], 'speed',[], 'speedDown',[], 'state',[], 'stateDown',[], 'stateBinary',[], 'bout',[] );
    if ~isempty(quadPath)
        quadData = load( quadPath{1} );
        loco.quad = quadData.quad_data;
        loco.Vinst = zeros(length(loco.quad), 1);
        if ~isempty(loco.Vinst)
            loco.Vinst(2:end) = diff(loco.quad);
            loco.Vinst(2) = 0;
            loco.Vinst = loco.Vinst*step_size*frameRate; % fps
            loco.Ainst = [NaN; diff(loco.Vinst)];
        end
    else
        return;
    end

    % Filter the instantaneous velocity and downsample to match scan rate
    Nframe = numel(loco.Vinst); % infoStruct.nframes
    gaussWidth = 1;
    gaussSigma = 0.26;
    gaussFilt = MakeGaussFilt( frameRate, gaussWidth, 0, gaussSigma, 'show',false);%MakeGaussFilt( gaussWidth, 0, gaussSigma, frameRate, false );
    loco.Vfilt = filtfilt( gaussFilt, 1, double(loco.Vinst) );
    loco.Vfilt = loco.Vfilt(1:Nframe); % # of frames in quad data not consistent with number in sbx?  infoStruct.nframes
    loco.Vdown = VectorDownsamp( loco.Vfilt(1:Nframe), Nplane ); % Downsample speed to match the sample number to the number of scans
    %loco.Afilt = filtfilt( gaussFilt, 1, double(loco.Ainst) );
    loco.Adown = VectorDownsamp( [NaN; diff(loco.Vfilt)], Nplane );  %VectorDownsamp( loco.Afilt(1:infoStruct.nframes), infoStruct.otlevels );

    % Convert velocity to speed
    loco.speed = abs( loco.Vfilt );
    loco.speedDown = abs( loco.Vdown );
    % Add other fields to be filled later
    loco.bout = nan(size(loco.speedDown)); %[];
    % Save time data (X-axis)
    loco.Tinst = (1/15)*((1:size(loco.Vinst,1))-1)'; 
    loco.Tdown = (Nplane/15)*((1:size(loco.Vdown,1))-1)';  
end

if show
    LocoData = figure;
    sp(1) = subplot(2,1,1);
    plot( loco.Tinst, loco.Vinst ); hold on;
    plot( loco.Tdown, loco.Vdown );
    ylabel('Velocity (cm/s)');
    sp(2) = subplot(2,1,2);
    plot( loco.Tinst, loco.Ainst ); hold on;
    plot( loco.Tdown, loco.Adown );
    ylabel('Acceleration (cm/s^2)');
    linkaxes(sp,'x');
    xlim([0,Inf]);
    savefig(LocoData, fullfile(exptDir,['LocoData.fig']));
end
end