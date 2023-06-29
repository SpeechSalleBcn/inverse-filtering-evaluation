function [] = saveResultsInCSV(GIFparams, GIFErrorsClass, timeMarks, invFilterResults, C)
% function [] = saveResultsInCSV(GIFparams, speechFile, GIFErrorsClass, timeMarks, a_gif, b_gif, C)
	% get filename
	% [~, speechFileName]=fileparts(speechFile);
	% Add analysis information to the table
	T = table;
	TP = table;

	% T.FileName = speechFile;
	T.FileName = C.filesData.speechFileName;
	TP.pulseId = find(timeMarks.isGCISelectedArray');
	% TP.FileName = repmat(speechFile, timeMarks.numGCISelected, 1);
	TP.FileName = repmat(C.filesData.speechFileName, timeMarks.numGCISelected, 1);

	% Parse information from the filename
	fileInfo = C.filesData.fileInfo;
	fileFields = C.filesData.fileFields;
	for i=1:length(fileFields)
	    if ~isnan(str2double(fileInfo.(fileFields{i})))
	        T.(fileFields{i}) = str2double(fileInfo.(fileFields{i}));
	        TP.(fileFields{i}) = repmat(str2double(fileInfo.(fileFields{i})), timeMarks.numGCISelected, 1);
	    elseif ~isempty(fileInfo.(fileFields{i}))
	        T.(fileFields{i}) = fileInfo.(fileFields{i});
	        TP.(fileFields{i}) = repmat(fileInfo.(fileFields{i}), timeMarks.numGCISelected, 1);
	    else
	        T.(fileFields{i}) = "";
	        TP.(fileFields{i}) = repmat("", timeMarks.numGCISelected, 1);
	    end
	end

	T.F0 = timeMarks.F0value;
	T.semitonesDifferenceGT = timeMarks.semitonesDifferenceGT;
	TP.F0value = repmat(timeMarks.F0value, timeMarks.numGCISelected, 1);
	TP.F0Pulse = timeMarks.F0Pulse(timeMarks.isGCISelectedArray)';
	TP.semitonesDifference = timeMarks.semitonesDifferencePulse(timeMarks.isGCISelectedArray)';
	TP.semitonesDifferenceGT = timeMarks.semitonesDifferencePulseGT(timeMarks.isGCISelectedArray)';

	T.GIFmethod = C.gif.method;
	TP.GIFmethod = repmat(C.gif.method, timeMarks.numGCISelected, 1);
	T.OptimMethod = C.optim.method;
	TP.OptimMethod = repmat(C.optim.method, timeMarks.numGCISelected, 1);
	if strcmp(C.optim.method, 'surrogate')
		T.NumEval = C.optim.MaxFunctionEvaluations;
		TP.NumEval = repmat(C.optim.MaxFunctionEvaluations, timeMarks.numGCISelected, 1);
    end
	T.VTorder = GIFparams.VTorder;
	T.GSorder = GIFparams.GSorder;
	T.LipRad = GIFparams.LipRad;
	TP.VTorder = repmat(GIFparams.VTorder, timeMarks.numGCISelected, 1);
	TP.GSorder = repmat(GIFparams.GSorder, timeMarks.numGCISelected, 1);
	TP.LipRad = repmat(GIFparams.LipRad, timeMarks.numGCISelected, 1);
	if (strcmp(C.gif.method,'Original-IAIF')|| strcmp(C.gif.method,'IOP-IAIF'))
	    T.HPflag = GIFparams.HPflag;
	    TP.HPflag = repmat(GIFparams.HPflag,timeMarks.numGCISelected,1);
	end
	if strcmp(C.gif.method,'IOP-IAIF')
	    T.Niter_IOP_IAIF = invFilterResults.numIters;
	    TP.Niter_IOP_IAIF = repmat(invFilterResults.numIters, timeMarks.numGCISelected, 1);
	end
	if contains(C.gif.method,'QCP')
	    T.DQ = GIFparams.DQ;
	    T.PQ = GIFparams.PQ;
	    T.RQ = GIFparams.RQ;
	    T.STcompensation = GIFparams.STcompensation;
	    TP.DQ = repmat(GIFparams.DQ, timeMarks.numGCISelected, 1);
	    TP.PQ = repmat(GIFparams.PQ, timeMarks.numGCISelected, 1);
	    TP.RQ = repmat(GIFparams.RQ, timeMarks.numGCISelected, 1);
	    TP.STcompensation = repmat(GIFparams.STcompensation, timeMarks.numGCISelected, 1);
	end

	[tableErrorsFile, tableErrorsPulse] = GIFErrorsClass.getAllErrorsTables();
	normParams = GIFErrorsClass.getNormParams();
	T.a = median(normParams.scale);
	T.b = median(normParams.dcOffset);
	T.timeLag = median(normParams.timeLag);
	TP.a = normParams.scale(timeMarks.isGCISelectedArray);
	TP.b = normParams.dcOffset(timeMarks.isGCISelectedArray);
	TP.timeLag = repmat(normParams.timeLag, timeMarks.numGCISelected, 1);

	T = [T, tableErrorsFile];
	TP = [TP, tableErrorsPulse];

	fileNameBase = [C.experiment.repository '_' C.gif.GIFMethod '_'];
	if C.experiment.optimization
		fileNameBase = [fileNameBase C.optim.method];
	else
		fileNameBase = [fileNameBase 'analysis'];
	end
	saveResultsFileName = [C.paths.mainPath filesep C.paths.experimentDir filesep C.paths.resultsDir filesep fileNameBase];

    writetable(T,[saveResultsFileName '_file.csv'], 'Delimiter', ',', 'WriteMode', 'append');
    writetable(TP,[saveResultsFileName '_pulse.csv'], 'Delimiter', ',', 'WriteMode', 'append');

end
