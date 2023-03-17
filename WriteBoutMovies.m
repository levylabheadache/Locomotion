function [boutStack, Tbout, boutSpeed, storeFrames] = WriteBoutMovies(expt, sbxInfo, Tscan, loco, periBout, movParam, varargin) % 
IP = inputParser;
addRequired( IP, 'expt', @isstruct )
addRequired( IP, 'sbxInfo', @isstruct )
addRequired( IP, 'Tscan', @iscell )
addRequired( IP, 'loco', @isstruct )
addRequired( IP, 'periBout', @isstruct ) 
addRequired( IP, 'movParam', @isstruct ) 
addOptional( IP, 'ROI', [], @isstruct )
addOptional( IP, 'axon', [], @isstruct )
addParameter( IP, 'run', [], @isnumeric )
addParameter( IP, 'bout', 1:1000, @isnumeric )
addParameter( IP, 'scalebar', [], @isstruct )
addParameter( IP, 'Toffset', 0, @isnumeric )
%addParameter( IP, 'overwrite', false, @islogical )
parse( IP, expt, sbxInfo, Tscan, loco, periBout, movParam, varargin{:} );  % , ROI 
ROI = IP.Results.ROI; Nroi = numel(ROI);
axon = IP.Results.axon;  Naxon = numel(axon);
if Naxon > 1
    axonColor = distinguishable_colors(Naxon);
else
    axonColor = [1,1,1]; % 0.5*
end
setRun = IP.Results.run;
setBout = IP.Results.bout;
%overwrite = IP.Results.overwrite;
GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
RGBopt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',true );
mkdir(movParam.dir);
if isfield(movParam, 'sbx') && ~isempty(movParam.sbx)
    sourceSbxPath = movParam.sbx; %[sbxInfo.dir, sbxInfo.exptName, '.', movParam.sourceSbx]; % [expt.dir, expt.name, '.', movParam.sourceSbx];
%elseif isfield(movParam, 'sourceSbx') && ~isempty(movParam.sourceSbx)
%    sourceSbxPath = [expt.dir, expt.name, '.', movParam.sourceSbx]; %[sbxInfo.dir, sbxInfo.exptName, '.', movParam.sourceSbx]; % 
else
    error('Need to specify sbx file');
end
if isempty(movParam.zProj), movParam.zProj = 1:expt.Nplane; end
Tcat = vertcat(Tscan{:});
if ~isfield(movParam, 'scalebar'), movParam.scalebar = []; end
if ~isfield(movParam, 'Toffset'), movParam.Toffset = 0; end
speedCat = vertcat( loco.speedDown );
stateCat = vertcat( loco.state );
if isempty(setRun)
    if strcmpi(movParam.boutType, 'csd')
        setRun = expt.csd;
    else
        setRun = 1:expt.Nruns;
    end
end
if ~isfield(movParam, 'FontSize') || isempty(movParam.FontSize), movParam.FontSize = 14; end
if ~isfield(movParam, 'tif') || isempty(movParam.tif), movParam.tif = false; end
NsetRun = numel(setRun);
tempFrameFile = strcat(movParam.dir, 'tempFrame.jpg'); %  tif
boutStack = cell(1,NsetRun); Tbout = cell(1,NsetRun); boutSpeed = cell(1,NsetRun); storeFrames = [];
for run = setRun %runs
    if isfield(periBout, 'Nbout')
        Nbout = periBout(run).Nbout;
    else
        Nbout = numel(periBout(run));
    end
    boutStack{run} = cell(1,Nbout); %tifPath = cell(1,Nbout);
    setRunBouts = intersect( setBout, 1:Nbout );
    for bout = setRunBouts  %setBout % 1:Nbout
        fprintf('\n[run, bout] = [%i,  %i]', run, bout);
        if ~isempty(movParam.name)
            boutFileName = sprintf('%s%s_%s_%s_run%i_bout%03d', movParam.dir, expt.name, movParam.name, movParam.regType, run, bout); % movParam.regType, 
        else
            boutFileName = sprintf('%s%s_%s_run%i_bout%03d', movParam.dir, expt.name, movParam.regType, run, bout); % movParam.regType, 
        end
        %boutFileName = sprintf('%s%s_run%i_bout%03d_%s_%s', movParam.dir, expt.name, run, bout, movParam.regType); % , movParam.boutType
        stackPath = strcat(boutFileName,'_stack.tif');
        % Determine bout's indices in concatenated movies frame
        [~,firstScan] = min( abs(Tcat - (periBout(run).Tstart(bout) - movParam.Tperi(1))) );  %Tcat(firstScan)
        [~,lastScan] = min( abs(Tcat - (periBout(run).Tstop(bout) + movParam.Tperi(2)) ) ); %Tcat(lastScan)
        NscanBout = lastScan-firstScan+1;
        if expt.Nplane == 1
            boutStack{run}{bout} = WriteSbxPlaneTif(sourceSbxPath, sbxInfo, movParam.zProj, 'firstScan',firstScan, 'Nscan',NscanBout, 'binT',movParam.binT, 'verbose',true, 'overwrite',true, 'monochrome',false); % , 'dir',movParam.dir     
        else
            boutStack{run}{bout} = WriteSbxZproj(sourceSbxPath, sbxInfo, 'z',movParam.zProj, 'firstScan',firstScan, 'Nscan',NscanBout, 'binT',movParam.binT, 'overwrite',true, 'monochrome',false); % , 'dir',movParam.dir
        end
        % Determine the speed and relative timing of each (possibly binned) frame/scan
        if movParam.binT > 1
            [catChunkLims, Nchunk] = MakeChunkLims(firstScan, lastScan, sbxInfo.totScan, 'size',movParam.binT ); % , chunkLength
            boutChunkLims = MakeChunkLims(1, NscanBout, NscanBout, 'size',movParam.binT ); %size(boutProjStack,3)
            for c = 1:Nchunk
                Tbout{run}{bout}(c) = mean( Tcat(catChunkLims(c,1):catChunkLims(c,2)) - periBout(run).Tstart(bout) );
                boutSpeed{run}{bout}(c) = mean( speedCat(catChunkLims(c,1):catChunkLims(c,2)) );% speedCat(firstScan:lastScan);
                boutState{run}{bout}(c) = mean( stateCat(catChunkLims(c,1):catChunkLims(c,2)) );
                %if expt.Nplane > 1,  boutStack{run}{bout}(:,:,c) = mean(boutProjStack(:,:,boutChunkLims(c,1):boutChunkLims(c,2)), 3);  end
            end
        else
            %if expt.Nplane > 1, boutStack{run}{bout} = boutProjStack; end
            Tbout{run}{bout} = Tcat(firstScan:lastScan) - periBout(run).Tstart(bout) - movParam.Toffset; %Toffset; % periBout.T{1}(
            boutSpeed{run}{bout} = speedCat(firstScan:lastScan);
            boutState{run}{bout} = stateCat(firstScan:lastScan);
        end
        boutStack{run}{bout} = uint16(boutStack{run}{bout}(movParam.edges(3)+1:end-movParam.edges(4), movParam.edges(1)+1:end-movParam.edges(2),:,:));
        % figure; plot(Tbout{run}{b}, boutSpeed{run}{b} )
        
        % Write the movie to tiff
        if movParam.tif
            fprintf('\nWriting %s...  ', stackPath); tic;
            saveastiff(boutStack{run}{bout}, stackPath, GrayOpt); toc
        end
        % Compress the data to 8 bit
        displayLims = prctile(boutStack{run}{bout}(:), movParam.displayPct);
        boutStack{run}{bout} = uint8( rescale(boutStack{run}{bout}, 0, 255, 'InputMin',displayLims(1), 'InputMax',displayLims(2)) );
    end
    if strcmpi(movParam.level, 'ind') || strcmpi(movParam.level, 'both')
        for bout = setRunBouts 
            % Make frames of annotated movie for AVI and Tif
            %displayLims = [prctile(boutStack{run}{b}(:),movParam.displayPct(1)), prctile(boutStack{run}{b}(:),movParam.displayPct(2))]; %[min(boutStack{run}{b}(:)), max(boutStack{run}{b}(:))];
            close all; clearvars storeFrames
            NperiScan = size(boutStack{run}{bout},3);
            movFig = figure('color','k','Units', 'normalized', 'Position', [0 0 1 1]); % 'WindowState','maximized',
            for z = flip(1:NperiScan)
                imshow(boutStack{run}{bout}(:,:,z), [0,255], 'border','tight'); hold on; % displayLims
                %MakeScaleBar( round(expt(X).umPerPixel*[50,50]), {get(gca,'Xlim'), get(gca,'Ylim')}, [0.1,0.95], [0,0], 'label',false, 'color','w' );
                if ~isempty(movParam.scalebar),  DrawScaleBar( movParam.scalebar );  end
                if ~isempty(ROI)
                    if ~isempty(axon)
                        for a = 1:Naxon
                            for roi = axon(a).ROI
                                plot( ROI(roi).footprintEdge(:,2)-movParam.edges(1)+1, ROI(roi).footprintEdge(:,1)-movParam.edges(3)+1, '.', 'color',axonColor(a,:), 'MarkerSize',1 );
                                %text(ROI(i).cent(1)-movParam.edges(1)+1, ROI(i).cent(2)-movParam.edges(3)+1, sprintf('%i-%i',i,a), 'HorizontalAlignment','center', 'color',axonColor(a,:), 'FontSize',10  ) 
                            end
                        end
                    else
                        roiColor = distinguishable_colors(Nroi);
                        for roi = 1:Nroi
                            plot( ROI(roi).footprintEdge(:,2)-movParam.edges(1)+1, ROI(roi).footprintEdge(:,1)-movParam.edges(3)+1, '.', 'color',roiColor(roi,:), 'MarkerSize',1 );
                            %text(ROI(i).cent(1)-movParam.edges(1)+1, ROI(i).cent(2)-movParam.edges(3)+1, sprintf('%i',i), 'HorizontalAlignment','center', 'color',roiColor(i,:), 'FontSize',10  ) 
                        end
                        
                    end
                end
                text( 0.9, 0.9, strcat(num2str(Tbout{run}{bout}(z), movParam.fmtSpec), ' s'), 'Units','normalized', 'color','w', 'HorizontalAlignment','center', 'FontSize',movParam.FontSize )
                if strcmpi(movParam.boutType, 'loco')
                    if boutState{run}{bout}(z) == 1
                        text( 0.9, 0.1, strcat( num2str(boutSpeed{run}{bout}(z), '%2.1f'), ' cm/s'), 'Units','normalized', 'color','r', 'HorizontalAlignment','center', 'FontSize',movParam.FontSize )
                    elseif boutState{run}{bout}(z) == 2
                        text( 0.9, 0.1, strcat( num2str(boutSpeed{run}{bout}(z), '%2.1f'), ' cm/s'), 'Units','normalized', 'color','g', 'HorizontalAlignment','center', 'FontSize',movParam.FontSize )
                    else
                        text( 0.9, 0.1, strcat( num2str(boutSpeed{run}{bout}(z), '%2.1f'), ' cm/s'), 'Units','normalized', 'color','w', 'HorizontalAlignment','center', 'FontSize',movParam.FontSize )
                    end
                end
                exportgraphics(gca, tempFrameFile, 'Resolution',300); % save image with set resolution  
                storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
                %pause(0.02); 
                cla;
            end
            close(movFig);
            % Write annotated Tif
            annotatedTif = zeros( [size(storeFrames(1).cdata),NperiScan], 'uint8');
            for z = flip(1:size(boutStack{run}{bout},3))
                annotatedTif(:,:,:,z) = storeFrames(z).cdata;
            end
            annotPath = strcat(boutFileName,'.tif');
            saveastiff(annotatedTif, annotPath, RGBopt )

            % Write mp4 (avi has issues with mac)
            mp4path = sprintf('%s%s_run%i_%s_%s%03d.mp4', movParam.dir, expt.name, run, movParam.boutType, movParam.regType, bout);
            fprintf('\nWriting %s...  ', mp4path); 
            writerObj = VideoWriter(mp4path, 'MPEG-4');
            writerObj.FrameRate = movParam.aviRate; % set the seconds per image
            open(writerObj);
            for z = 1:numel(storeFrames)
                frame = storeFrames(z);    
                writeVideo(writerObj, frame);
            end
            close(writerObj);
            toc
        end
    end
end

if strcmpi(movParam.level, 'cat') || strcmpi(movParam.level, 'both')
    boutStackCat = cell(1,NsetRun); boutSpeedCat = cell(1,NsetRun); boutTimeCat = cell(1,NsetRun);
    for r = 1:NsetRun
        boutStackCat{r} = cat(3, boutStack{setRun(r)}{:});
        boutSpeedCat{r} = cat(1, reshape(boutSpeed{setRun(r)}{:}, [], 1));
        boutTimeCat{r} = cat(1, reshape(Tbout{setRun(r)}{:}, [], 1));
    end
    boutStackCat = cat(3, boutStackCat{:});
    boutSpeedCat = cat(1, boutSpeedCat{:});
    boutTimeCat = cat(1, boutTimeCat{:});
    NcatFrame = size(boutStackCat,3);
    
    % Write the concatenated movie to tiff
    %{
    catTifPath = sprintf('%s%s_%s_%s_concat.tif', movParam.dir, expt.name, movParam.regType, movParam.boutType);% %strcat(movParam.dir, 'CSD_', exptName, '_', movParam.sourceSbx,'.tif');
    GrayOpt = struct('overwrite',true, 'message',false, 'append',false, 'big',true, 'color',false );
    fprintf('\nWriting %s...', catTifPath); tic;
    saveastiff(boutStackCat, catTifPath, GrayOpt); toc
    %}
    
    % Make frames for MP4
    tic
    displayLims = prctile(boutStackCat(:), movParam.displayPct);
    close all; clearvars storeFrames
    movFig = figure('color','k', 'Units','normalized', 'OuterPosition', [0 0 1 1]); % 'WindowState','maximized',
    w = waitbar(0, 'Writing concatenated movie');
    for z = flip(1:NcatFrame) %  313
        imshow(boutStackCat(:,:,z), displayLims); hold on;
        if ~isempty(movParam.scalebar),  DrawScaleBar( movParam.scalebar );  end
        if ~isempty(ROI)
            if ~isempty(axon)
                for a = 1:Naxon
                    for roi = flip(axon(a).ROI)
                        plot( ROI(roi).footprintEdge(:,2)-movParam.edges(1)+1, ROI(roi).footprintEdge(:,1)-movParam.edges(3)+1, '.', 'color',axonColor(a,:), 'MarkerSize',1 );
                        %text(ROI(i).cent(1)-movParam.edges(1)+1, ROI(i).cent(2)-movParam.edges(3)+1, sprintf('%i-%i',i,a), 'HorizontalAlignment','center', 'color',axonColor(a,:), 'FontSize',movParam.FontSize  ) 
                        %pause
                    end
                end
            else
                roiColor = distinguishable_colors(Nroi);
                for roi = 1:Nroi
                    plot( ROI(roi).footprintEdge(:,2)-movParam.edges(1)+1, ROI(roi).footprintEdge(:,1)-movParam.edges(3)+1, '.', 'color',roiColor(roi,:), 'MarkerSize',1 );
                    %text(ROI(i).cent(1)-movParam.edges(1)+1, ROI(i).cent(2)-movParam.edges(3)+1, sprintf('%i',i), 'HorizontalAlignment','center', 'color',roiColor(i,:), 'FontSize',movParam.FontSize  ) 
                end
            end
        end
        text( 0.9, 0.9, strcat(num2str(boutTimeCat(z), movParam.fmtSpec), ' s'), 'Units','normalized', 'color','w', 'HorizontalAlignment','center', 'FontSize',movParam.FontSize )
        if strcmpi(movParam.boutType, 'loco')
            text( 0.9, 0.1, strcat(num2str(boutSpeedCat(z), movParam.fmtSpec), ' cm/s'), 'Units','normalized', 'color','w', 'HorizontalAlignment','center', 'FontSize',movParam.FontSize )
        end
        exportgraphics(gca, tempFrameFile, 'Resolution',100); % save image with set resolution  
        storeFrames(z) = im2frame(imread(tempFrameFile));  % getframe(aviFig);  % convert image to frame
        waitbar((NcatFrame-z+1)/NcatFrame, w)
        %pause%(0.01); 
        cla;
    end
    close(movFig);
    delete(tempFrameFile); delete(w)
    toc

    % Write movie
    mp4Path = sprintf('%s%s_%s_%s_concat.avi', movParam.dir, expt.name, movParam.regType, movParam.boutType);
    fprintf('\nWriting %s...  ', mp4Path); tic;
    writerObj = VideoWriter(mp4Path, 'MPEG-4');
    writerObj.FrameRate = movParam.aviRate; % set the seconds per image
    open(writerObj);
    for z = 1:numel(storeFrames)
        frame = storeFrames(z);    
        writeVideo(writerObj, frame);
    end
    close(writerObj);
    toc 
end
end
