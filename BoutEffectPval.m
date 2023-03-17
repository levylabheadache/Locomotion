function pEffect = BoutEffectPval(periBout, periStat, periParam, fieldName)
    pEffect = nan(size(periStat.(fieldName).effect));  
    tempEffect = nan([size(periStat.(fieldName).effect), periParam.Nshuff]);
    for s = 1:periParam.Nshuff
        tempStat = BoutEffect(cellfun(@ScrambleColumns, periBout.(fieldName), 'UniformOutput',false), periBout.preScan, periBout.boutScan, periBout.postScan);
        tempEffect(:,:,:,s) = tempStat.effect;
    end
    for b = 1:size(tempEffect,1)
        for r = 1:size(tempEffect,2)
            for t = 1:3
                pEffect(b,r,t) = (sum(abs(tempEffect(b,r,t,:)) > abs(periStat.(fieldName).effect(b,r,t)))+1)/(periParam.Nshuff+1);
            end
        end
    end
end