%
% Copyright (C) 2023 HER Speech, La Salle - Universitat Ramon Llull
%
% This file is part of the project Inverse Filtering Evaluation repository
% https://github.com/SpeechSalleBcn/inverse-filtering-evaluation
%
% License
% This file is under the LGPL license,  you can
% redistribute it and/or modify it under the terms of the GNU Lesser General 
% Public License as published by the Free Software Foundation, either version 3 
% of the License, or (at your option) any later version. This file is
% distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
% details.
% https://www.gnu.org/licenses/lgpl-3.0.html
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% Authors
% Marc Freixes - marc.freixes@salle.url.edu
% Luis Joglar-Ongay - luis.joglar@salle.url.edu
% Joan Claudi Socoró - joanclaudi.socoro@salle.url.edu
%

function executeInverseFilterAnalysis(experimentConfigPath)
% Script that performs optimization (or analysis from a given set of parameters) of GIF methods on a speech corpus.
% Input:
% experimentConfigPath: path to the config.ini file for the experiment.
%                       Details about this config file can be found the readme.

    currentPath = pwd;
    cd '../'
    startup;
    cd(currentPath);

    initConfig = loadConfigFile(experimentConfigPath);
    experimentConfigPathParts = strsplit(experimentConfigPath, filesep);
    filesepIdx = strfind(experimentConfigPath, filesep);
    initConfig.paths.mainPath = experimentConfigPath(1 : filesepIdx(end - 1)-1);
    initConfig.paths.experimentDir = experimentConfigPathParts{end - 1};
    logPath = [initConfig.paths.mainPath filesep ...
                    initConfig.paths.experimentDir filesep ...
                    initConfig.paths.logDir];
    if ~exist(logPath, 'dir')
       mkdir(logPath)
    end
    diary(fullfile(logPath, initConfig.paths.logFileName))

    configFiles = dir([initConfig.paths.mainPath filesep ...
                       initConfig.paths.experimentDir filesep ...
                       initConfig.paths.invFilterConfigDir]);
    configFiles = configFiles([configFiles.isdir] == 0);

    for n = 1 : length(configFiles)
        C = loadConfigFile([configFiles(n).folder filesep configFiles(n).name]);
        C.experiment = initConfig.experiment;
        C.paths = initConfig.paths;
        C.gif.GIFMethod = C.gif.method;
        if strcmp(C.gif.GIFMethod, 'QCP') && (sum(C.gif.STcompensation) > 0); C.gif.GIFMethod = 'ST-QCP'; end
        if C.experiment.optimization
            disp(['OPTIMIZATION of ', [C.paths.mainPath filesep ...
                                       C.paths.corpusDir], ' using ', configFiles(n).name,' ...'])
            inverseFilterAnalysisCorpus(C);
        else
            disp(['ANALYSIS of ', [C.paths.mainPath filesep ...
                                   C.paths.corpusDir], ' using ', configFiles(n).name,' ...'])
            paramFiles = dir([C.paths.mainPath filesep ...
                              C.paths.experimentDir filesep ...
                              C.paths.invFilterConfigDir filesep ...
                              C.paths.paramFilesDir]);
            paramFiles = paramFiles([paramFiles.isdir] == 0);
            paramFilesMethod = paramFiles(contains({paramFiles.name}, ['_' C.gif.GIFMethod '_' ]));
            C.paramsTable = table;
            for m = 1 : length(paramFilesMethod)
                paramsTable = readtable([C.paths.mainPath filesep ...
                                           C.paths.experimentDir filesep ... 
                                           C.paths.invFilterConfigDir filesep ...
                                           C.paths.paramFilesDir filesep ...
                                           paramFilesMethod(m).name]);
                C.paramsTable = [C.paramsTable; paramsTable];
            end
            inverseFilterAnalysisCorpus(C);
        end
        disp(n + " of " + length(configFiles) + ", GIF Methods processed.")
    end
    diary off
end