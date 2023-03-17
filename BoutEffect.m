function summStruct = BoutEffect( dataMat, preScan, boutScan, postScan )
    scanMean = @(x,s)( mean(x(s,:,:),1,'omitnan') );
    meanPre = cellfun( scanMean, dataMat, preScan, 'uniformOutput',false );
    meanRun = cellfun( scanMean, dataMat, boutScan, 'uniformOutput',false );
    meanPost = cellfun( scanMean, dataMat, postScan, 'uniformOutput',false );
    summStruct.pre = vertcat(meanPre{:});
    summStruct.run = vertcat(meanRun{:});
    summStruct.post = vertcat(meanPost{:});
    summStruct.still = (summStruct.pre + summStruct.post)/2;
    summStruct.effect(:,:,1) = (summStruct.run - summStruct.pre)./summStruct.pre; % onset effect
    summStruct.effect(:,:,2) = (summStruct.post - summStruct.pre)./summStruct.pre; % offset effect
    %summStruct.effect(:,:,2) = (summStruct.run - summStruct.post)./summStruct.post; % offset effect
    summStruct.effect(:,:,3) = (summStruct.run - summStruct.still)./summStruct.still; % overall effect
end