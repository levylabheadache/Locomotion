%% Train mouse-specific models

% Clear any previous variables in the Workspace and Command Window to start fresh
clear; clc; close all; 

% TODO --- Set the directory of your 'sbx' files
mainDir = 'D:\2photon\Simone\Simone_Macrophages\'; % 'D:\2photon\Simone\'; %'D:\2photon\';
mouse = {'Pf4Ai162-17'}; %  'Pf4Ai162-1' {'DL102','DL89','DL115','DL118','DL68','DL75','DL112','DL117','DL72','DL67','DL122'}; %  D:\2photon\DL102
Nmouse = numel(mouse);
mouseLoco = cell(1,Nmouse);
tic;
for m = 1:Nmouse
    mouseDir = strcat(mainDir, mouse{m},'\');
    [~,mouseExptDir] = FileFinder(mouseDir, 'type',0);
    c = 0;
    for f = 1:size(mouseExptDir,1)
        [mouseExptRunName, mouseExptRunDir] = FileFinder(mouseExptDir{f}, 'type',0, 'contains','run' );  % 
        if ~isempty(mouseExptRunDir)
            tempExptName = mouseExptRunName{1};
            uscoreLoc = strfind(tempExptName,'_');
            tempDate = tempExptName(uscoreLoc(1)+1:uscoreLoc(2)-1);
            fovLoc = strfind(tempExptName, 'FOV');
            if isempty(fovLoc)
                tempFOV = [];     
            else
                tempFOV = tempExptName(fovLoc+3:uscoreLoc(3)-1); % str2double( tempExptName(fovLoc+3:uscoreLoc(3)-1) );
            end       
            for r = 1:size(mouseExptRunDir,1)
                c = c+1;
                tempRun = str2double(tempExptName(strfind(tempExptName,'run')+3:end));
                try
                    runInfoTemp = MakeInfoStruct( mainDir, mouse{m}, tempDate, tempRun, tempFOV );
                    mouseLoco{m}(c) = GetLocoData( runInfoTemp ); % 
                catch 
                    fprintf('\nFailed: f = %i:  %s %s %i', f, mouse{m}, tempDate, tempRun);
                end
            end
        else
            fprintf('\nf = %i: %s has no runs', f, mouseExptDir{f});
        end
    end
    TrainLocoHMM(mouseLoco(m), 'Nstate',2, 'int',1, 'var','velocity', 'dir',mouseDir, 'name',mouse{m}, 'show',true, 'overwrite',true);  %  'speed'
    toc
end

%% Fit a master HMM using all available data
[transEst, emitEst, speedDiscreteExpt] = TrainLocoHMM(loco);