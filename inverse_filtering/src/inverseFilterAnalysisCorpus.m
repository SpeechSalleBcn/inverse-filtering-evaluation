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

function inverseFilterAnalysisCorpus(C)
% inverseFilterAnalysisCorpus(path, fileNameFormat, cfgFile, fileSaveResults, plotResult)
%
% This function perfoms optimization of a input corpus (variable path is the input corpus 
% foldername) of speech signals that contain its corresponding ground-truth glottal waveform 
% signals following the configuration given by cfgFile.
%
% The optimization hyperparameters and assessment results are saved into to XLS files 
% fileSaveResults and filnameP in Table format.
%
% If plotResult = true results of optimized GIF method are viewed comparing the
% GT glottal waveform and its estimation for the optimized GIF method.

    tStart = tic;
    corpusPath = [C.paths.mainPath filesep C.paths.corpusDir];
    DIR = dir(corpusPath);
    for n = 1:length(DIR)
        if (isfolder([corpusPath filesep DIR(n).name]))&&(~strcmp(DIR(n).name,'.'))&&(~strcmp(DIR(n).name,'..'))
            disp(['Processing folder ',[corpusPath filesep DIR(n).name],' ...'])
            speechFiles = dir ([corpusPath filesep DIR(n).name filesep '*.wav']);
            for k = 1:length(speechFiles)
                speechFile = [DIR(n).name filesep speechFiles(k).name];
                disp(['   ',speechFiles(k).name])
                try
                    [speechData, timeMarks, C] = getSpeechFileData(corpusPath, speechFile, C);

                    if C.experiment.optimization
                        GIFparams = optimization(speechData, timeMarks, C);
                    else
                        GIFparams = getGIFparams(C);
                    end

                    [invFilterResults] = applyInverseFilter(GIFparams, speechData, timeMarks, C);

                    %% Plot GIF result
                    if C.analysis.plotResult
                        plotGIFResults(speechData.glotFlowGT , invFilterResults.glotFlow, C, GIFparams, 'Optimization');
                    end

                    errors = GIFErrors(invFilterResults, speechData, timeMarks, C);
                    errors.computeAllErrors();
                    saveResultsInCSV(GIFparams, errors, timeMarks, invFilterResults, C);

                catch ME
                    msgText = getReport(ME,'extended','hyperlinks','on');
                    disp(msgText);
                end
                disp(k + " of " + length(speechFiles) + ", speech files processed with " + C.gif.GIFMethod+ " from " + DIR(n).name)
            end
        end
    end
    tEnd = toc(tStart);
    disp(['Total analysis time: ' num2str(tEnd) ' seconds'])
end
