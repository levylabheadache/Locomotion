respStruct = struct('Nbout',[], 'meanEffect',[], 'exc',[], 'Nexc',NaN, 'excFrac',NaN, 'excEffect',[], 'neut',[], 'Nneut',NaN, 'neutFrac',NaN, 'neutEffect',[], 'inh',[], 'Ninh',NaN, 'inhFrac',NaN, 'inhEffect',[]); % 'normMean',[], 'valence',[],
boutRespStruct = struct('pre',respStruct, 'post',respStruct, 'change',struct('e2e',[], 'Ne2e',NaN, 'e2n',[], 'Ne2n',NaN, 'e2i',[], 'Ne2i',NaN, ...
    'n2e',[], 'Nn2e',NaN, 'n2n',[], 'Nn2n',NaN, 'n2i',[], 'Nn2i',NaN, 'i2e',[], 'Ni2e',NaN, 'i2n',[], 'Ni2n',NaN, 'i2i',[], 'Ni2i',NaN));
boutResponse = repmat(boutRespStruct, 1, Nexpt);  periBout = cell(1,Nexpt); periParam = cell(1,Nexpt); periStat = cell(1,Nexpt);
plotYlim = [-0.3,0.5]; TL = 0.001;
%BoutResponse_Classification = figure('WindowState','maximized', 'Color','w');   % 'Units','inches', 'Position',[6, 5, 2.1, 1.1] , 'Color','w'
for x = xPresent %find(~cellfun(@isempty, Tscan)) %xPresxent %x3Dcsd %
    savePath = sprintf('%s%s_bouts.mat', expt{x}.dir, expt{x}.name);
    if exist(savePath,'file')
        fprintf('\nLoading %s', savePath)
        load(savePath, 'periData','boutResp'); % ,'minIso','minMeanEffect','pMax','minBoutDur','baseDur','initIso',
        boutResponse(x) = boutResp;  periBout{x} = periData.bout;  periStat{x} = periData.stat;  periParam{x} = periData.param;
    else
        % Set parameters
        minBoutDur = 2;
        baseDur = 10;
        initIso = [0,0];
        minIso = 10;
        minVelOn = 4;
        minMeanEffect = 0.05; % 0.05;
        pMax = 0.05;
        % Get peri-locomotion, stillness, and peri-CSD bout data
        for runs = flip(1:expt{x}.Nruns)
            [periBout{x}(runs), periParam{x}(runs), loco{x}(runs)] = PeriLoco3D(expt{x}, Tscan{x}{runs}, loco{x}(runs), deform{x}(runs), fluor{x}(runs).F.ROI, defVars, 'base',baseDur, 'run',minBoutDur, 'iso',initIso); %periParam(x)
            %InspectPeriDeform3D( expt{x}, periBout{x}(runs), {'fluor','scaleMag'} ); %defVars
        end
        % Determine effect of loco bouts pre-CSD
        preEffect_speed = []; 
        preEffect_fluor_onset = []; preEffect_fluor_offset = [];  preEffect_trans = []; preEffect_scale = []; preEffect_shear = []; preEffect_shift = [];
        for runs = expt{x}.preRuns %1:2
            if periBout{x}(runs).Nbout > 0
                % get effect sizes
                tempEffect_speed = periBout{x}(runs).stat.speed.run;
                tempEffect_fluor_onset = periBout{x}(runs).stat.fluor.effect(:,:,1); % effect of bout on fluor, relative to pre-running period
                tempEffect_fluor_offset = periBout{x}(runs).stat.fluor.effect(:,:,2); % effect of bout offset on fluor, relative to pre-running period
                tempEffect_trans = periBout{x}(runs).stat.transMag.effect(:,:,1); % effect of bout on scaling, relative to pre-running period
                tempEffect_scale = periBout{x}(runs).stat.scaleMag.effect(:,:,1); % effect of bout on scaling, relative to pre-running period
                tempEffect_shear = periBout{x}(runs).stat.shearMag.effect(:,:,1); % effect of bout on shearing, relative to pre-running period
                tempEffect_shift = periBout{x}(runs).stat.shiftZ.effect(:,:,1); % % effect of bout on scaling, relative to pre-running period
                % suppress poorly isolated bouts
                bBad = periBout{x}(runs).iso(:,1) < minIso; % any(periBout{x}(run).iso < minIso, 2);  % find poorly isolated bouts
                tempEffect_speed(bBad) = NaN;
                tempEffect_fluor_onset(bBad,:) = NaN; % suppress poorly isolated bouts
                tempEffect_fluor_offset(bBad,:) = NaN;
                tempEffect_trans(bBad,:) = NaN;
                tempEffect_scale(bBad,:) = NaN; % suppress poorly isolated bouts
                tempEffect_shear(bBad,:) = NaN; % suppress poorly isolated bouts
                tempEffect_shift(bBad,:) = NaN; % suppress poorly isolated bouts
                % join effects across runs
                preEffect_speed = vertcat(preEffect_speed, tempEffect_speed);
                preEffect_fluor_onset = vertcat(preEffect_fluor_onset, tempEffect_fluor_onset);
                preEffect_fluor_offset = vertcat(preEffect_fluor_offset, tempEffect_fluor_offset);
                preEffect_trans = vertcat(preEffect_trans, tempEffect_trans);
                preEffect_scale = vertcat(preEffect_scale, tempEffect_scale);
                preEffect_shear = vertcat(preEffect_shear, tempEffect_shear);
                preEffect_shift = vertcat(preEffect_shift, tempEffect_shift);
            end
        end
        [~, pEffect_pre_onset] = ttest(preEffect_fluor_onset);
        [~, pEffect_pre_offset] = ttest(preEffect_fluor_offset);
        meanPreEffect_onset = mean(preEffect_fluor_onset, 1, 'omitnan');
        meanPreEffect_offset = mean(preEffect_fluor_offset, 1, 'omitnan');
        %JitterPlot(preEffect_fluor_onset);

        % Save the results to boutRespStruct
        boutResponse(x).pre.Nbout = mode(sum(~isnan(preEffect_fluor_onset),1));
        boutResponse(x).pre.speedEffect = preEffect_speed;
        % fluor
        boutResponse(x).pre.boutEffect = cat(3, preEffect_fluor_onset, preEffect_fluor_offset); %preEffect_fluor_onset; %
        boutResponse(x).pre.meanEffect = cat(1, meanPreEffect_onset, meanPreEffect_offset);
        % deformation
        boutResponse(x).pre.transEffect = preEffect_trans;
        boutResponse(x).pre.scaleEffect = preEffect_scale;
        boutResponse(x).pre.shearEffect = preEffect_shear;
        boutResponse(x).pre.shiftEffect = preEffect_shift;
        % subtypes
        boutResponse(x).pre.exc = find( pEffect_pre_onset < pMax & meanPreEffect_onset > minMeanEffect );
        [~, tempSortInd] = sort( meanPreEffect_onset(boutResponse(x).pre.exc), 'descend' ); % sort from most positive to most negative
        boutResponse(x).pre.exc = boutResponse(x).pre.exc(tempSortInd); % plot( meanEffect(boutRespStruct(x).exc) )
        boutResponse(x).pre.Nexc = numel(boutResponse(x).pre.exc);
        boutResponse(x).pre.excFrac = boutResponse(x).pre.Nexc/expt{x}.Nroi;
        boutResponse(x).pre.excEffect = boutResponse(x).pre.meanEffect(:,boutResponse(x).pre.exc); % sorted

        boutResponse(x).pre.inh = find( pEffect_pre_onset < pMax & meanPreEffect_onset < -minMeanEffect );
        [~, tempSortInd] = sort( meanPreEffect_onset(boutResponse(x).pre.inh), 'descend' ); % sort from most positive to most negative
        boutResponse(x).pre.inh = boutResponse(x).pre.inh(tempSortInd); % plot( meanEffect(boutRespStruct(x).inh) )
        boutResponse(x).pre.Ninh = numel(boutResponse(x).pre.inh);
        boutResponse(x).pre.inhFrac = boutResponse(x).pre.Ninh/expt{x}.Nroi;
        boutResponse(x).pre.inhEffect = boutResponse(x).pre.meanEffect(:,boutResponse(x).pre.inh); % sorted

        boutResponse(x).pre.neut = setdiff( 1:expt{x}.Nroi, [boutResponse(x).pre.exc, boutResponse(x).pre.inh]); %find( pEffect_pre > pMax );
        [~, tempSortInd] = sort( meanPreEffect_onset(boutResponse(x).pre.neut), 'descend' ); % sort from most positive to most negative
        boutResponse(x).pre.neut = boutResponse(x).pre.neut(tempSortInd); % plot( meanEffect(boutRespStruct(x).neut) )
        boutResponse(x).pre.Nneut = numel(boutResponse(x).pre.neut);
        boutResponse(x).pre.neutFrac = boutResponse(x).pre.Nneut/expt{x}.Nroi;
        boutResponse(x).pre.neutEffect = boutResponse(x).pre.meanEffect(:,boutResponse(x).pre.neut ); % sorted

        boutResponse(x).pre.offset = find( pEffect_pre_offset < pMax & meanPreEffect_offset > minMeanEffect ); % which units are activated on offset?
        boutResponse(x).pre.Noffset = numel(boutResponse(x).pre.offset);

        % 



        % Determine effect of loco bouts post-CSD
        % {
        if ~isnan(expt{x}.csd) && sum([periBout{x}(expt{x}.csd:end).Nbout]) > 0
            postEffect_speed = []; 
            postEffect_fluor_onset = []; postEffect_fluor_offset = []; postEffect_trans = []; postEffect_scale = []; postEffect_shear = []; postEffect_shift = [];
            for runs = expt{x}.csd+1:expt{x}.Nruns %3:4
                if periBout{x}(runs).Nbout > 0
                    bBad = periBout{x}(runs).iso(:,1) < minIso; % find poorly isolated bouts
                    tempEffect_speed = periBout{x}(runs).stat.speed.run;
                    tempEffect_fluor_onset = periBout{x}(runs).stat.fluor.effect(:,:,1); % effect of bout on fluor, relative to pre-running period
                    tempEffect_trans = periBout{x}(runs).stat.transMag.effect(:,:,1); % effect of bout on scaling, relative to pre-running period
                    tempEffect_scale = periBout{x}(runs).stat.scaleMag.effect(:,:,1); % effect of bout on scaling, relative to pre-running period
                    tempEffect_shear = periBout{x}(runs).stat.shearMag.effect(:,:,1); % effect of bout on shearing, relative to pre-running period
                    tempEffect_shift = periBout{x}(runs).stat.shiftZ.effect(:,:,1); % effect of bout on scaling, relative to pre-running period
                    % suppress poorly isolated bouts
                    tempEffect_speed(bBad) = NaN;
                    tempEffect_fluor_onset(bBad,:) = NaN; % suppress poorly isolated bouts
                    tempEffect_fluor_offset(bBad,:) = NaN;
                    tempEffect_trans(bBad,:) = NaN;
                    tempEffect_scale(bBad,:) = NaN; % suppress poorly isolated bouts
                    tempEffect_shear(bBad,:) = NaN; % suppress poorly isolated bouts
                    tempEffect_shift(bBad,:) = NaN; % suppress poorly isolated bouts
                    % join effects across runs
                    postEffect_speed = vertcat(postEffect_speed, tempEffect_speed);
                    postEffect_fluor_onset = vertcat(postEffect_fluor_onset, tempEffect_fluor_onset);
                    postEffect_fluor_offset = vertcat(postEffect_fluor_offset, tempEffect_fluor_offset);
                    postEffect_trans = vertcat(postEffect_trans, tempEffect_trans);
                    postEffect_scale = vertcat(postEffect_scale, tempEffect_scale);
                    postEffect_shear = vertcat(postEffect_shear, tempEffect_shear);
                    postEffect_shift = vertcat(postEffect_shift, tempEffect_shift);

                end
            end
            [~, pEffect_post] = ttest(postEffect_fluor_onset);
            meanPostEffect = mean(postEffect_fluor_onset, 1, 'omitnan');
            boutResponse(x).post.Nbout = mode(sum(~isnan(postEffect_fluor_onset),1));
            boutResponse(x).post.speedEffect = postEffect_speed;
            boutResponse(x).post.boutEffect = postEffect_fluor_onset;
            boutResponse(x).post.meanEffect = meanPostEffect;

            % deformation
            boutResponse(x).post.transEffect = postEffect_trans;
            boutResponse(x).post.scaleEffect = postEffect_scale;
            boutResponse(x).post.shearEffect = postEffect_shear;
            boutResponse(x).post.shiftEffect = postEffect_shift;

            boutResponse(x).post.exc = find( pEffect_post < pMax & meanPostEffect > minMeanEffect );
            [~, tempSortInd] = sort( meanPostEffect(boutResponse(x).post.exc), 'descend' ); % sort from most positive to most negative
            boutResponse(x).post.exc = boutResponse(x).post.exc(tempSortInd); % plot( meanEffect(boutRespStruct(x).exc) )
            boutResponse(x).post.Nexc = numel(boutResponse(x).post.exc);
            boutResponse(x).post.excFrac = boutResponse(x).post.Nexc/expt{x}.Nroi;
            boutResponse(x).post.excEffect = boutResponse(x).post.meanEffect( boutResponse(x).post.exc ); % sorted

            boutResponse(x).post.inh = find( pEffect_post < pMax & meanPostEffect < -minMeanEffect );
            [~, tempSortInd] = sort( meanPostEffect(boutResponse(x).post.inh), 'descend' ); % sort from most positive to most negative
            boutResponse(x).post.inh = boutResponse(x).post.inh(tempSortInd); % plot( meanEffect(boutRespStruct(x).inh) )
            boutResponse(x).post.Ninh = numel(boutResponse(x).post.inh);
            boutResponse(x).post.inhFrac = boutResponse(x).post.Ninh/expt{x}.Nroi;
            boutResponse(x).post.inhEffect = boutResponse(x).post.meanEffect( boutResponse(x).post.inh ); % sorted

            boutResponse(x).post.neut = setdiff( 1:expt{x}.Nroi, [boutResponse(x).post.exc, boutResponse(x).post.inh]); %find( pEffect_post > pMax );
            [~, tempSortInd] = sort( meanPostEffect(boutResponse(x).post.neut), 'descend' ); % sort from most positive to most negative
            boutResponse(x).post.neut = boutResponse(x).post.neut(tempSortInd); % plot( meanEffect(boutRespStruct(x).neut) )
            boutResponse(x).post.Nneut = numel(boutResponse(x).post.neut);
            boutResponse(x).post.neutFrac = boutResponse(x).post.Nneut/expt{x}.Nroi;
            boutResponse(x).post.neutEffect = boutResponse(x).post.meanEffect( boutResponse(x).post.neut ); % sorted

            boutResponse(x).change.e2e = intersect(boutResponse(x).pre.exc, boutResponse(x).post.exc);
            boutResponse(x).change.Ne2e = numel(boutResponse(x).change.e2e);
            boutResponse(x).change.e2eFrac = boutResponse(x).change.Ne2e/expt{x}.Nroi;

            boutResponse(x).change.e2n = intersect(boutResponse(x).pre.exc, boutResponse(x).post.neut);
            boutResponse(x).change.Ne2n = numel(boutResponse(x).change.e2n);
            boutResponse(x).change.e2nFrac = boutResponse(x).change.Ne2n/expt{x}.Nroi;

            boutResponse(x).change.e2i = intersect(boutResponse(x).pre.exc, boutResponse(x).post.inh);
            boutResponse(x).change.Ne2i = numel(boutResponse(x).change.e2i);
            boutResponse(x).change.e2iFrac = boutResponse(x).change.Ne2i/expt{x}.Nroi;

            boutResponse(x).change.n2e = intersect(boutResponse(x).pre.neut, boutResponse(x).post.exc);
            boutResponse(x).change.Nn2e = numel(boutResponse(x).change.n2e);
            boutResponse(x).change.n2eFrac = boutResponse(x).change.Nn2e/expt{x}.Nroi;

            boutResponse(x).change.n2n = intersect(boutResponse(x).pre.neut, boutResponse(x).post.neut);
            boutResponse(x).change.Nn2n = numel(boutResponse(x).change.n2n);
            boutResponse(x).change.n2nFrac = boutResponse(x).change.Nn2n/expt{x}.Nroi;

            boutResponse(x).change.n2i = intersect(boutResponse(x).pre.neut, boutResponse(x).post.inh);
            boutResponse(x).change.Nn2i = numel(boutResponse(x).change.n2i);
            boutResponse(x).change.n2iFrac = boutResponse(x).change.Nn2i/expt{x}.Nroi;

            boutResponse(x).change.i2e = intersect(boutResponse(x).pre.inh, boutResponse(x).post.exc);
            boutResponse(x).change.Ni2e = numel(boutResponse(x).change.i2e);
            boutResponse(x).change.i2eFrac = boutResponse(x).change.Ni2e/expt{x}.Nroi;

            boutResponse(x).change.i2n = intersect(boutResponse(x).pre.inh, boutResponse(x).post.neut);
            boutResponse(x).change.Ni2n = numel(boutResponse(x).change.i2n);
            boutResponse(x).change.i2nFrac = boutResponse(x).change.Ni2n/expt{x}.Nroi;

            boutResponse(x).change.i2i = intersect(boutResponse(x).pre.inh, boutResponse(x).post.inh);
            boutResponse(x).change.Ni2i = numel(boutResponse(x).change.i2i);
            boutResponse(x).change.i2iFrac = boutResponse(x).change.Ni2i/expt{x}.Nroi;

            boutResponse(x).change.fracMat = [boutResponse(x).change.e2eFrac, boutResponse(x).change.n2eFrac, boutResponse(x).change.i2eFrac; ...
                boutResponse(x).change.e2nFrac, boutResponse(x).change.n2nFrac, boutResponse(x).change.i2nFrac; ...
                boutResponse(x).change.e2iFrac, boutResponse(x).change.n2iFrac, boutResponse(x).change.i2iFrac];
        end
        % Save the results 
        periData.bout = periBout{x}; periData.stat = periStat{x}; periData.param = periParam{x}; boutResp = boutResponse(x);
        fprintf('\nSaving %s', savePath)
        save(savePath, 'initIso','baseDur','minBoutDur', 'periData', 'minIso','minMeanEffect','pMax','boutResp'); % 'periBout','periParam','periStat'
    end

    %}

    %{
    if ~isnan(expt{x}.csd), sp(1) = subplot(2,2,1); end
    line([0,expt{x}.Nroi+1], [0,0], 'color','k'); hold on;
    JitterPlot(preEffect_fluor_onset(:,[boutResponse(x).pre.exc, boutResponse(x).pre.neut, boutResponse(x).pre.inh]), 'monochrome',0.7);
    %tempYlim = get(gca,'Ylim');
    line([0,expt{x}.Nroi+1], minMeanEffect*[1,1], 'color','r', 'lineStyle','--'); hold on;
    line([0,expt{x}.Nroi+1], -minMeanEffect*[1,1], 'color','r', 'lineStyle','--'); hold on;
    line((boutResponse(x).pre.Nexc+0.5)*[1,1], plotYlim, 'color','g');
    line((boutResponse(x).pre.Nexc+boutResponse(x).pre.Nneut+0.5)*[1,1], plotYlim, 'color','r');
    xlim([0,expt{x}.Nroi+1]); 
    ylim(plotYlim);
    set(gca, 'Xtick',[], 'XtickLabel',[], 'Ytick',-0.25:0.25:0.5, 'YtickLabel',[], 'TickDir','out', 'TickLength',[TL,0]) %  -0.4:0.2:0.8 , 'Position',[0.02,0.02,0.96,0.96]
    ylabel('Normalized Bout Effect: (Fbout-Fbase)/Fbase'); xlabel('ROI (sorted, within group, by pre-CSD mean effect)');
    title( sprintf('x = %i: n = %i well-isolated bouts (>= %2.1f s running, >= %2.1f s isolation)', x, boutResponse(x).pre.Nbout, periParam(x).run, minIso) );  

    if ~isnan(expt{x}.csd) 
        sp(2) = subplot(2,2,3); 
        line([0,expt{x}.Nroi+1], [0,0], 'color','k'); hold on;
        JitterPlot(postEffect(:,[boutResponse(x).post.exc, boutResponse(x).post.neut, boutResponse(x).post.inh]), 'monochrome',0.7);
        %JitterPlot(postEffect(:,[boutResponse(x).pre.exc, boutResponse(x).pre.neut, boutResponse(x).pre.inh]), 'monochrome',0.7);
        %tempYlim = get(gca,'Ylim');
        line([0,expt{x}.Nroi+1], minMeanEffect*[1,1], 'color','r', 'lineStyle','--'); hold on;
        line([0,expt{x}.Nroi+1], -minMeanEffect*[1,1], 'color','r', 'lineStyle','--'); hold on;
        line((boutResponse(x).post.Nexc+0.5)*[1,1], plotYlim, 'color','g');
        line((boutResponse(x).post.Nexc+boutResponse(x).post.Nneut+0.5)*[1,1], plotYlim, 'color','r');
        set(gca, 'Xtick',[], 'XtickLabel',[], 'Ytick',-0.25:0.25:0.5, 'YtickLabel',[], 'TickDir','out', 'TickLength',[TL,0]) %  -0.4:0.2:0.8  , 'Position',[0.02,0.02,0.96,0.96]
        ylabel('Normalized Bout Effect: (Fbout-Fbase)/Fbase'); xlabel('ROI (sorted, within group, by post-CSD mean effect)'); % pre-
        title( sprintf('x = %i: n = %i well-isolated bouts (>= %2.1f s running, >= %2.1f s isolation)', x, boutResponse(x).post.Nbout, periParam(x).run, minIso) );
        linkaxes(sp, 'xy');
        xlim([0,expt{x}.Nroi+1]);
        ylim(plotYlim);
        
        subplot(2,2,[2,4]);
        imagesc(boutResponse(x).change.fracMat); axis square;
        caxis([0,1]);
        set(gca,'Xtick',1:3, 'Ytick',1:3, 'XtickLabel',{'Exc','Neut','Inh'}, 'YtickLabel',{'Exc','Neut','Inh'});
        xlabel('Pre-CSD'); ylabel('Post-CSD'); title('Transition frequencies');
        impixelinfo;
    end
    pause; clf;
    %figPath = sprintf('%sBoutResponse_Classification_%s.tif', fig2Dir, expt{x}.name );
    %exportgraphics(BoutResponse_Classification, figPath, 'Resolution',300); fprintf('\nSaved %s', figPath);
    %}
end
tempPre = [boutResponse(xPresent).pre]; %tempPost = [boutResponse(xPresent).post];
boutResponseTypePct = 100*[[tempPre.excFrac]', [tempPre.inhFrac]'];
[~,pBoutResponseFrac] = ttest(boutResponseTypePct(:,1), boutResponseTypePct(:,2));
%{
figDir = 'D:\MATLAB\Dura\Figures\';
savePath = sprintf('%sBoutResponse_Ctrl.mat', figDir);
save(savePath, 'boutResponse', 'boutResponseTypePct', 'expt','minIso','minMeanEffect','pMax','periBout','periStat');
fprintf('\nSaving %s', savePath);
%}

%%
reboundFrac = nan(Nexpt,1);
for x = xPresent
    reboundFrac(x) = numel(intersect(boutResponse(x).pre.offset, boutResponse(x).pre.inh))/boutResponse(x).pre.Ninh;
end

%% Show onset/offset effects for each subype
figure;
for x = xPresent
    subplot(1,3,1); cla;
    JitterPlot(boutResponse(x).pre.excEffect', 'paired',true, 'title','Activated', 'new',false); hold on;
    line(get(gca,'Xlim'), 0.05*[1,1], 'color','r', 'linestyle','--');
    line(get(gca,'Xlim'), -0.05*[1,1], 'color','r', 'linestyle','--');
    set(gca,'Xtick',[1,2], 'XtickLabel',{'Onset','Offset'})
    ylabel('Pre-CSD Bout effect')
    axis square;
    subplot(1,3,2); cla;
    JitterPlot(boutResponse(x).pre.neutEffect', 'paired',true, 'title','Neutral', 'new',false)
    line(get(gca,'Xlim'), 0.05*[1,1], 'color','r', 'linestyle','--');
    line(get(gca,'Xlim'), -0.05*[1,1], 'color','r', 'linestyle','--');
    axis square;
    subplot(1,3,3); cla;
    JitterPlot(boutResponse(x).pre.inhEffect', 'paired',true, 'title','Suppressed', 'new',false)
    line(get(gca,'Xlim'), 0.05*[1,1], 'color','r', 'linestyle','--');
    line(get(gca,'Xlim'), -0.05*[1,1], 'color','r', 'linestyle','--');
    axis square;
    pause;
end
%% Show pre/postCSD effects for each subype
figure;
for x = xPresent
    cla
    JitterPlot([boutResponse(x).pre.meanEffect(1,boutResponse(x).pre.exc)', boutResponse(x).post.meanEffect(1,boutResponse(x).pre.exc)'], 'paired',true, 'title','Activated', 'new',false); hold on;
    line(get(gca,'Xlim'), 0.05*[1,1], 'color','r', 'linestyle','--');
    line(get(gca,'Xlim'), -0.05*[1,1], 'color','r', 'linestyle','--');
    set(gca,'Xtick',[1,2], 'XtickLabel',{'Onset','Offset'})
    ylabel('Pre-CSD Bout effect')
    title(sprintf('%s: %i pre-CSD bouts,  %i post.  Activated', expt{x}.name, boutResponse(x).pre.Nbout, boutResponse(x).post.Nbout), 'Interpreter','none')
    axis square;
    pause;
end

%%

for x = x3Dcsd % xPresent
    sp(1) = subplot(1,2,1);
    line([0,1], [0,1], 'color','r', 'lineStyle','--'); hold on;
    plot( boutResponse(x).pre.excFrac, boutResponse(x).post.excFrac, '.');

    sp(2) = subplot(1,2,2);
    line([0,1], [0,1], 'color','r', 'lineStyle','--'); hold on;
    plot( boutResponse(x).pre.inhFrac, boutResponse(x).post.inhFrac, '.'); hold on;
end
subplot(sp(1)); xlim([0,1]); ylim([0,1]); axis square; xlabel('Pre-CSD frac'); ylabel('Post-CSD frac'); title('Activated');
subplot(sp(2)); xlim([0,1]); ylim([0,1]); axis square; xlabel('Pre-CSD frac'); ylabel('Post-CSD frac'); title('Dectivated');

%% Chi square scrap

for x = x3Dcsd
    % Expected number of excited ROI
    nullFrac = (boutResponse(x).pre.Nexc + boutResponse(x).post.Nexc)/(2*expt{x}.Nroi);
    Nexc_null = nullFrac*expt{x}.Nroi;
    [~, p_exc(x)] = chi2gof([1,2,3,4], 'freq',[boutResponse(x).pre.Nexc, expt{x}.Nroi-boutResponse(x).pre.Nexc, boutResponse(x).post.Nexc, expt{x}.Nroi-boutResponse(x).post.Nexc], ...
        'expected',[Nexc_null, expt{x}.Nroi-Nexc_null, Nexc_null, expt{x}.Nroi-Nexc_null], 'nparams',2);
end


% Observed data
n1 = 51; N1 = 8193;
n2 = 74; N2 = 8201;
% Pooled estimate of proportion
p0 = (n1+n2) / (N1+N2);
% Expected counts under H0 (null hypothesis)
n10 = N1 * p0;
n20 = N2 * p0;
% Chi-square test, by hand
observed = [n1 N1-n1 n2 N2-n2];
expected = [n10 N1-n10 n20 N2-n20];
[h,p,stats] = chi2gof([1 2 3 4],'freq',observed,'expected',expected,'ctrs',[1 2 3 4],'nparams',1);



ct = randsample([0 1],190,true,[0.49 0.51]);
[h,p,stats] = chi2gof(ct,'Expected',[95 95]);
[h,p,stats] = chi2gof([1,2],'Expected',[95 95], 'freq',[97,93])

%%

for x = xPresent
    boutResponse(x).pre.Ninh %+ boutResponse(x).pre.Ninh
    pause;
end
exptPres = [expt{xPresent}]
exptPres.Nroi

allPre = [boutResponse(xPresent).pre];
[allPre.Ninh]
