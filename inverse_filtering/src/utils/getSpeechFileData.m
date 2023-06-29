function [speechData, timeMarks, C] = getSpeechFileData(path, speechFile, C)

    [x,fs] = audioread([path filesep speechFile]);
    speechData.speech = x(:, 1);
    speechData.glotFlowGT = x(:, 2);
    speechData.fs = fs;

    % Parse information from the filename
    [~, speechFileName]=fileparts(speechFile);
    C.filesData.fileInfo = regexp(speechFileName, C.experiment.filesFormat, 'names');
    C.filesData.speechFileName = speechFileName;
    C.filesData.fileFields = fieldnames(C.filesData.fileInfo);

    %% Time marks computation
    timeMarks = computeTimeMarks(speechData.speech, speechData.fs, C);
end
