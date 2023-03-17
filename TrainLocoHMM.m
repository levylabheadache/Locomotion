function [transEst, emitEst, discreteExptData] = TrainLocoHMM(loco, varargin)
% Determine locomotive state using a Hidden Markov Model
IP = inputParser;
addRequired( IP, 'loco', @iscell )
addParameter( IP, 'overwrite', false, @islogical )
addParameter( IP, 'int', 0.1, @isnumeric )
addParameter( IP, 'var', 'velocity', @ischar ) % train the model using velocity or speeed?
addParameter( IP, 'Nstate', 2, @isnumeric ) % Nstate = 2; % still/minor running, major running
addParameter( IP, 'stateName', {'Still','Running'}, @iscell )
addParameter( IP, 'transGuess', [0.9, 0.1; 0.5, 0.5], @isnumeric ) %transGuess = [0.9, 0.1; 0.5, 0.5];
addParameter( IP, 'emitGuess', [0,5; 0.5,5], @isnumeric ) % row one: mean emission values for each state. row two: std of emission for each state (in cm/s, not bins)   [0,1;2,5]
addParameter( IP, 'dir', 'D:\2photon\', @ischar )
addParameter( IP, 'name', '', @ischar )
addParameter( IP, 'show', false, @islogical )
parse( IP, loco, varargin{:} ); 
overwrite = IP.Results.overwrite;
discreteInt = IP.Results.int;
modelVar = IP.Results.var; %
Nstate = IP.Results.Nstate;
stateName = IP.Results.stateName; % stateName = {'Still','Running'};
transGuess = IP.Results.transGuess;
emitGuess = IP.Results.emitGuess;
show = IP.Results.show;
saveDir = IP.Results.dir;
saveName = IP.Results.name;
savePath = sprintf('%s%s_LocoHMM_%s_%istate.mat', saveDir, saveName, modelVar, Nstate);

tic
if ~exist(savePath, 'file') || overwrite
    % Discretize speed or velocity data (matlab can't handle continuous HMM) and pool across experiments/runs
    Nexpt = numel(loco);
    discreteExptData = cell(1,Nexpt);
    if strcmpi(modelVar, 'speed')
        fprintf('\nDiscretizing speed data ');
        emitRange = 0:discreteInt:40; 
        for x = find(~cellfun(@isempty, loco))
            for r = find(~cellfun(@isempty, {loco{x}.quad})) %1%:
                loco{x}(r).discrete = imquantize( vertcat( loco{x}(r).speed ), emitRange ); % Down  speedDiscreteExpt 
            end
            discreteExptData{x} = vertcat(loco{x}.discrete); % pool runs
        end
    elseif strcmpi(modelVar, 'velocity')
        fprintf('\nDiscretizing velocity data ');
        emitRange = -10:discreteInt:40;
        for x = find(~cellfun(@isempty, loco))
            for r = find(~cellfun(@isempty, {loco{x}.quad})) %1%:
                loco{x}(r).discrete = imquantize( vertcat( loco{x}(r).Vfilt ), emitRange ); % Down  speedDiscreteExpt 
            end
            discreteExptData{x} = vertcat(loco{x}.discrete); % pool runs
        end
    else
        error('invalid training variable %s', modelVar)
    end
    Nemit = numel(emitRange);
    discretePool = vertcat(discreteExptData{:})';  % pool experiments = vertcat(discreteExptData{:})';  % pool experiments
    toc

    % Generate initial guess emission distributions
    fprintf('\nGenerating initial guess emission distributions');
    emitGuessPDF = nan(Nstate, Nemit); 
    %figure;
    for n = flip(1:Nstate)
        tempGauss = normpdf(emitRange, emitGuess(1,n), emitGuess(2,n) );
        emitGuessPDF(n,:) = tempGauss/sum(tempGauss); 
        %plot( emitRange, emitGuessPDF(n,:) ); hold on;
    end
    toc

    % Train HMM on pooled data
    fprintf('\nTraining hidden Markov model: %i samples from %i experiments', numel(discretePool), Nexpt);
    [transEst, emitEst] = hmmtrain( discretePool, transGuess, emitGuessPDF, 'verbose',true );
    toc

    % Save the results
    fprintf('\nSaving %s', savePath);
    save(savePath); % ,  'transEst', 'emitEst', 'emitGuessPDF', 'transGuess'
    toc
else
    fprintf('\nLoading %s', savePath);
    load(savePath);
    toc
end

% Compare the estimated model parameters to the guess
if show
    poolState = hmmviterbi( discretePool, transEst, emitEst )';
    
    close all;
    figure('Units','normalized','OuterPosition',[0,0,1,1]);   
    subplot(2,3,1:2); 
    plot( emitRange', emitGuessPDF', '--' ); hold on; 
    plot( emitRange', emitEst' ); 
    xlabel(modelVar); ylabel('Emission Probability Density');
    xlim([-3,10]); %xlim([0,20]);
    
    subplot(2,3,4:5);
    yyaxis left; 
    plot( emitRange(discretePool) ); %plot( loco(r).speed ); 
    ylabel(modelVar); xlabel('Frame');
    yyaxis right; 
    plot( poolState-1, 'color',[1,0,0,0.2] ); %plot( loco(r).state -1, 'color',[1,0,0,0.4] );  % , 'LineStyle','--'
    set(gca,'Ytick',0:Nstate-1, 'box','off');
    ylim([0,Nstate-1+0.2]);
    ylabel('Locomotive State');
    xlim([-Inf,Inf]);
    
    subplot(2,3,3); imagesc( transGuess); axis square; title('Guessed Transition Probabilities'); set(gca,'Xtick',1:Nstate, 'XtickLabel',stateName,'Ytick',1:Nstate, 'YtickLabel',stateName);
    text(1,1, sprintf('%2.4f',transGuess(1,1)), 'color','k', 'HorizontalAlignment','center' );
    text(2,1, sprintf('%2.4f',transGuess(2,1)), 'color','w', 'HorizontalAlignment','center' );
    text(1,2, sprintf('%2.4f',transGuess(1,2)), 'color','w', 'HorizontalAlignment','center' );
    text(2,2, sprintf('%2.4f',transGuess(2,2)), 'color','k', 'HorizontalAlignment','center' );
    
    subplot(2,3,6); imagesc( transEst); axis square; title('Estimated Transition Probabilities'); set(gca,'Xtick',1:Nstate, 'XtickLabel',stateName,'Ytick',1:Nstate, 'YtickLabel',stateName);
    text(1,1, sprintf('%2.4f',transEst(1,1)), 'color','k', 'HorizontalAlignment','center' );
    text(2,1, sprintf('%2.4f',transEst(2,1)), 'color','w', 'HorizontalAlignment','center' );
    text(1,2, sprintf('%2.4f',transEst(1,2)), 'color','w', 'HorizontalAlignment','center' );
    text(2,2, sprintf('%2.4f',transEst(2,2)), 'color','k', 'HorizontalAlignment','center' );
    colormap('gray'); 
    impixelinfo;

    %{
    figure('Units','normalized','OuterPosition',[0,0,0.5,1])
    subplot(2,1,1);
    
    subplot(2,1,2);
    yyaxis left; 
    plot( vertcat(loco.speedDown) ); %plot( loco(r).speed ); 
    ylabel('Downsampled Speed'); xlabel('Frame');
    yyaxis right; 
    plot( vertcat(loco.stateDown) -1, 'color',[1,0,0,0.3] ); %plot( loco(r).state -1, 'color',[1,0,0,0.4] );  % , 'LineStyle','--'
    set(gca,'Ytick',0:Nstate-1, 'box','off');
    ylim([0,Nstate-1+0.2]);
    ylabel('Downsampled Locomotive State');
    xlim([-Inf,Inf]);
    
    for r = 1:numel(loco)
        yyaxis left; 
        plot( vertcat(loco.speed) ); %plot( loco(r).speed ); 
        ylabel('Downsampled Speed'); xlabel('Frame');

        yyaxis right; 
        plot( vertcat(loco.state) -1, 'color',[1,0,0,0.4] ); %plot( loco(r).state -1, 'color',[1,0,0,0.4] );  % , 'LineStyle','--'
        set(gca,'Ytick',0:Nstate-1, 'box','off');
        ylim([0,Nstate-1+0.2]);
        ylabel('Locomotive State');
        xlim([-Inf,Inf]);
        pause; cla;
    end
    %}
end

end