% function [times, pmarks, ne_select, gci_t2, gci_sample, gci_pos] = computeTimeMarks(speech, fs, C)
function [timeMarks] = computeTimeMarks(speech, fs, C)

% Returns
% timeMarks.analysisFramesSamples   - (times)
% timeMarks.pitchMarks              - (pmarks)
% timeMarks.selectedSamplesIdxs     - (ne_select)
% timeMarks.gcisSamples             - (gci_sample)
% timeMarks.isGCISelectedArray      - (gci_pos)
% timeMarks.gcisTimes               - (gci_t)
% timeMarks.numGCISelected          - (length(gci_t2))

% Other Variables
% selectedPercent                   - (te_select)
% selectedGCITimes                  - (gci_t2) - deleted

% Usages
% gci_sample -> gcis
% nPulses -> length(gci_t2)
% gci_pos -> pulsesSel

if isfield(C.error, 'selectedPercent')
    selectedPercent = C.error.selectedPercent;
else
    selectedPercent = [0 100];
end

timeMarks.f0Min = C.gci.F0min;
timeMarks.f0Max = C.gci.F0max;
if isfield(C.gci, 'F0MarginFactor')
    if isfield(C.filesData.fileInfo, 'f0')
        timeMarks.f0Min = str2double(C.filesData.fileInfo.f0) * 2 ^ (- C.gci.F0MarginFactor / 12);
        timeMarks.f0Max = str2double(C.filesData.fileInfo.f0) * 2 ^ (C.gci.F0MarginFactor / 12);
    else
        warning('Error: File name does not contain F0 value');
    end
end

%% Pitch comuptation
f0 = pitch(speech, fs, 'Range', [timeMarks.f0Min, timeMarks.f0Max], 'Method', C.gci.F0_method);
timeMarks.F0value = median(f0);

%% GCIs computation
switch C.gci.GCImethod
    case 'SEDREAMS'
        timeMarks.gcisTimes = gci_sedreams(speech, fs, timeMarks.F0value, 1);
        timeMarks.gcisSamples = round(timeMarks.gcisTimes * fs) + 1;
    case 'GLOAT'        %(Drugman, T and Dutoit, T: "Glottal Closure andOpening Instant Detection from Speech Signals", 2009)
        [timeMarks.gcisSamples, goi] = gci(speech, timeMarks.F0value, fs);
        timeMarks.gcisTimes = timeMarks.gcisSamples / fs;
end

periodPulses = [0 timeMarks.gcisTimes];
timeMarks.F0Pulse = 1./(periodPulses(2:end) - periodPulses(1:end-1));
timeMarks.semitonesDifferenceGT = 12 * log2(timeMarks.F0value / str2num(C.filesData.fileInfo.f0));
timeMarks.semitonesDifferencePulse = 12 * log2(timeMarks.F0Pulse / timeMarks.F0value);
timeMarks.semitonesDifferencePulseGT = 12 * log2(timeMarks.F0Pulse / str2num(C.filesData.fileInfo.f0));

timeMarks.pitchMarks = [1 ...
                        round( timeMarks.gcisSamples(1 : end - 1) ...  % period between gcis in samples
                            + (C.ola.phi / 100) ...
                            * (timeMarks.gcisSamples(2 : end) ...
                               - timeMarks.gcisSamples(1 : end - 1) ...
                              ) ...
                        ) ...
                        length(speech) ...
                        ];
timeMarks.pitchMarks = [timeMarks.pitchMarks(1 : end - 1); timeMarks.pitchMarks(2 : end)];

% GCI adapted to the selected percentage interval (selectedPercent)
timeMarks.selectedSamplesIdxs = (round((selectedPercent(1) / 100) * length(speech)) + 1) ...
                                : (round((selectedPercent(2) / 100) * length(speech)));
timeMarks.isGCISelectedArray = ((timeMarks.gcisTimes > ((timeMarks.selectedSamplesIdxs(1) - 1) / fs)) ...
                            & (timeMarks.gcisTimes < (timeMarks.selectedSamplesIdxs(end) - 1) / fs));

timeMarks.numGCISelected = length(timeMarks.gcisTimes(timeMarks.isGCISelectedArray));

switch C.ola.WindType
    case 'Constant Frame Rate'
        wind = round(C.ola.WinLength * fs);
        step = round(C.ola.WinShift * fs);
        np = ceil((length(speech * fs) - wind) / step);   % Number of windows
        timeMarks.analysisFramesSamples = [1 : step : (step * np + 1)];
        timeMarks.analysisFramesSamples = [timeMarks.analysisFramesSamples; wind : step : (wind + step * np + 1)];
        % set last time to the length of the speech signal
        timeMarks.analysisFramesSamples(2, end)=length(speech);
    case 'Pitch Synchronous'
        % each frame includes 2 pulses
        timeMarks.analysisFramesSamples = [timeMarks.pitchMarks(1 : end - 2); timeMarks.pitchMarks(3 : end)];
end

end
