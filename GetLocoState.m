function loco = GetLocoState(expt, loco, varargin)
% Determine locomotive state using a Hidden Markov Model
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'loco', @isstruct )
addParameter( IP, 'var', 'velocity', @ischar )
addParameter( IP, 'Nstate', 2, @isnumeric ) % Nstate = 2; % still/minor running, major running
addParameter( IP, 'dir', 'D:\2photon\', @ischar )
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'show', false, @islogical )
addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, loco, varargin{:} ); 
loadDir = IP.Results.dir;
loadName = IP.Results.name;
modelVar = IP.Results.var;
Nstate = IP.Results.Nstate;
show = IP.Results.show;
overwrite = IP.Results.overwrite;
modelPath = sprintf('%s%s_LocoHMM_%s_%istate.mat', loadDir, loadName, modelVar, Nstate);
%Nexpt = numel(loco);

% Load master HMM
%modelPath = FileFind( loadDir, 'mat', false, @(x)(contains(x, modelVar) && contains(x, strcat(num2str(Nstate),'state'))) ); 
%modelPath = modelPath{1,2};
%loadPath = sprintf('%sLocoHMM_%istate.mat', loadDir, Nstate);
if exist(modelPath, 'file')
    fprintf('\nLoading %s... \n', modelPath);% tic;
    load(modelPath, 'emitRange', 'transEst', 'transGuess', 'emitEst', 'emitGuessPDF', 'stateName'); %   %, 'Nemit' 
    %toc
else
    warning('model does not exist! - returning');
    return;
end

% Apply the HMM to this data
for runs = 1:numel(loco)
    % Fit individual speed/velocity data sets to the HMM trained from the pooled data
    if strcmpi(modelVar, 'speed')
        loco(runs).state = hmmviterbi( imquantize(loco(runs).speed, emitRange)', transEst, emitEst )';
        %loco(r).stateDown = hmmviterbi(imquantize( loco(r).speedDown , emitRange)', transEst, emitEst )'; 
    elseif strcmpi(modelVar, 'velocity')
        try
            loco(runs).state = hmmviterbi( imquantize(loco(runs).Vfilt, emitRange)', transEst, emitEst )';
            Nplane = round(15.49/expt.scanRate);
            scanStates = reshape(loco(runs).state(1:expt.Nscan(runs)*Nplane), Nplane, [])';  % 1:expt.Nscan(r)*expt.Nplane  1:floor(length(loco(runs).state)/expt.Nplane)
            % imshow( scanStates', [] )
            loco(runs).stateDown = mode(scanStates,2); % max( scanStates, [], 2); % down sample for volume scan using modal or maximal state within the scan
        catch
            fprintf('Failed to calculated loco(%i).state, calculating stateDown instead', runs)
            loco(runs).stateDown = hmmviterbi( imquantize(loco(runs).Vdown , emitRange)', transEst, emitEst )';
        end
    end   
    loco(runs).stateBinary = false( size(loco(runs).stateDown, 1), Nstate );
    for n = 1:Nstate  
        loco(runs).stateBinary(:,n) = loco(runs).stateDown == n;  
    end      
end

% Compare the estimated model parameters to the guess
if show
    close all;
    HMM_model_details = figure('Units','normalized','OuterPosition',[0,0,0.5,1]);
    subplot(3,1,1); plot( emitGuessPDF', '--' ); hold on; 
    plot( emitRange', emitGuessPDF', '--' ); hold on; 
    plot( emitRange', emitEst' ); 
    xlabel(modelVar); ylabel('Emission Probability Density');
    xlim([-3,10]); 
    
    subplot(3,1,2); imagesc( transGuess); axis square; title('Guessed Transition Probabilities'); set(gca,'Xtick',1:Nstate, 'XtickLabel',stateName,'Ytick',1:Nstate, 'YtickLabel',stateName);
    text(1,1, sprintf('%2.4f',transGuess(1,1)), 'color','k', 'HorizontalAlignment','center' );
    text(2,1, sprintf('%2.4f',transGuess(2,1)), 'color','w', 'HorizontalAlignment','center' );
    text(1,2, sprintf('%2.4f',transGuess(1,2)), 'color','w', 'HorizontalAlignment','center' );
    text(2,2, sprintf('%2.4f',transGuess(2,2)), 'color','k', 'HorizontalAlignment','center' );
    
    subplot(3,1,3); imagesc( transEst); axis square; title('Estimated Transition Probabilities'); set(gca,'Xtick',1:Nstate, 'XtickLabel',stateName,'Ytick',1:Nstate, 'YtickLabel',stateName);
    text(1,1, sprintf('%2.4f',transEst(1,1)), 'color','k', 'HorizontalAlignment','center' );
    text(2,1, sprintf('%2.4f',transEst(2,1)), 'color','w', 'HorizontalAlignment','center' );
    text(1,2, sprintf('%2.4f',transEst(1,2)), 'color','w', 'HorizontalAlignment','center' );
    text(2,2, sprintf('%2.4f',transEst(2,2)), 'color','k', 'HorizontalAlignment','center' );
    colormap('gray'); 
    impixelinfo;

    figPath = sprintf('%s%s_HMM_model_details', expt.dir, expt.name); %.pdf
    if ~exist(figPath, 'file') || overwrite 
        fprintf('\nSaving %s', figPath);
        saveas(HMM_model_details, figPath)
        %exportgraphics(HMM_model_details, figPath, 'Resolution',300); 
    end
     
    HMM_model_results = figure('Units','normalized','OuterPosition',[0.5,0,0.5,1]);  
    sp(1) = subplot(2,1,1);
    yyaxis left; 
    if strcmpi(modelVar, 'speed')
        plot(vertcat(loco.Tinst), vertcat(loco.speed) );
    elseif strcmpi(modelVar, 'velocity')
        plot(vertcat(loco.Tinst), vertcat(loco.Vfilt) ); 
    end
    ylabel(modelVar); 
    xlabel('Time (frames)');
    
    yyaxis right; 
    plot(vertcat(loco.Tinst), vertcat(loco.state) -1, 'color',[1,0,0,0.2] ); %plot( loco(r).state -1, 'color',[1,0,0,0.4] );  % , 'LineStyle','--'
    set(gca,'Ytick',0:Nstate-1, 'box','off');
    ylim([0,Nstate-1+0.2]);
    ylabel('Locomotive State');
    xlim([-Inf,Inf]);
    
    sp(2) = subplot(2,1,2);
    yyaxis left; 
    if strcmpi(modelVar, 'speed')
        plot(vertcat(loco.Tdown), vertcat(loco.speedDown) );
    elseif strcmpi(modelVar, 'velocity')
        plot(vertcat(loco.Tdown), vertcat(loco.Vdown) ); 
    end
    ylabel(modelVar);  xlabel('Time (scans)');
    yyaxis right; 
    plot(vertcat(loco.Tdown), vertcat(loco.stateDown) -1, 'color',[1,0,0,0.2] ); %plot( loco(r).state -1, 'color',[1,0,0,0.4] );  % , 'LineStyle','--'
    set(gca,'Ytick',0:Nstate-1, 'box','off');
    ylim([0,Nstate-1+0.2]);
    ylabel('Downsampled Locomotive State');
    xlim([-Inf,Inf]);
    linkaxes(sp,'xy');
    
    % save the figure
    figPath = sprintf('%s%s_HMM_model_results', expt.dir, expt.name);
    if ~exist(figPath, 'file') || overwrite
        fprintf('\nSaving %s', figPath);
        saveas(HMM_model_results, figPath)
    end
end
end