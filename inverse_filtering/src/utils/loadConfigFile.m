function C = loadConfigFile(cfgFile)
% function C = loadSingingCfg(cfgFile)
% Load Configuration file
%   inputs.
%       cfgFile: ini configuration file
%   outputs. 
%       C: struct where configuration is loaded

    %% LOAD INI FILE
    ini=IniConfig();
    ini.ReadFile(cfgFile);
    iniSect = ini.GetSections();
    iniSect = strrep(iniSect,'[',''); iniSect = strrep(iniSect,']','');
    for i=1:length(iniSect)
        [keys, count_keys] = ini.GetKeys(iniSect{i});
        for j=1:count_keys
            eval(['C.' iniSect{i} '.' keys{j} ' = ini.GetValues(''' iniSect{i} ''',''' keys{j} ''');']);
        end
    end
    %% ADDED FOR BACKWARD COMPATIBILITY
    if all(isfield(C,{'gci','ola','gif','error','optim'}))
        if ~isfield(C.error,'selectedPercent'); C.error.selectedPercent=C.error.te_select; C.error=rmfield(C.error,'te_select');end
        if ~isfield(C.error,'relativeErrorRMSE'); C.error.relativeErrorRMSE=C.error.flagNormalizedError; C.error=rmfield(C.error,'flagNormalizedError');end
        if ~isfield(C,'analysis'); C.analysis.plotResult=false;end
    end
end
