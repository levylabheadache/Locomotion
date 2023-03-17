function [effectCI, effectSig] = BoutEffectConfInt(periBout, periStat, periParam, fieldName)

scrambleEffect = cell(1,periParam.Nshuff); % nan([size(periStat.(fieldName).effect), periParam.Nshuff]);
parfor s = 1:periParam.Nshuff
    scrambleStat = BoutEffect(cellfun(@ScrambleColumns, periBout.(fieldName), 'UniformOutput',false), periBout.preScan, periBout.boutScan, periBout.postScan);
    scrambleEffect{s} = scrambleStat.effect; % (:,:,:,s)
end
scrambleEffect = cat(4,scrambleEffect{:});
% Calculate confidence intervals (5-95%) on effect sizes
effectCI = prctile(scrambleEffect, [5,95], 4); % nan([size(periStat.(fieldName).effect), 2]);
effectSig = zeros(size(periStat.(fieldName).effect)); % [-1 if negative significant effect, 0 if no significant effect, +1 if positive significant effect]
effectSig(periStat.(fieldName).effect < effectCI(:,:,:,1)) = -1;
effectSig(periStat.(fieldName).effect > effectCI(:,:,:,2)) = 1;
toc;
%{
for bout = 1:size(scrambleEffect,1)
    for roi = 1:size(scrambleEffect,2)
        for Type = 1:3
            effectCI(bout,roi,Type,1) = prctile( scrambleEffect(bout,roi,Type,:), 5 );
            effectCI(bout,roi,Type,2) = prctile( scrambleEffect(bout,roi,Type,:), 95 );
            if periStat.(fieldName).effect(bout,roi,Type) < effectCI(bout,roi,Type,1)
                effectSig(bout,roi,Type) = -1;
            elseif periStat.(fieldName).effect(bout,roi,Type) > effectCI(bout,roi,Type,1)
                effectSig(bout,roi,Type) = 1;
            end
        end
    end
end
%}
end