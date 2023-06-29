function [invFilterResults] = timeAlignmentGlotFlow(invFilterResults, speechData, C)
    % TODO: adjust maxlag to fs and f0. Maybe by maxlag = fs / f0 * param (param=?0.5)
    maxlag = 30;
    [xc, lag] = xcorr(invFilterResults.glotFlow, speechData.glotFlowGT, maxlag, 'unbiased');
    [~, idx] = max(xc);
    invFilterResults.timeLag = lag(idx);
    invFilterResults.glotFlowTimeAligned = invFilterResults.glotFlow;
    invFilterResults.glotFlowDTimeAligned = invFilterResults.glotFlowD;
    if ~iscolumn(invFilterResults.glotFlowTimeAligned)
        invFilterResults.glotFlowTimeAligned = transpose(invFilterResults.glotFlowTimeAligned);
    end
    if ~iscolumn(invFilterResults.glotFlowDTimeAligned)
        invFilterResults.glotFlowDTimeAligned = transpose(invFilterResults.glotFlowDTimeAligned);
    end
    if invFilterResults.timeLag < 0
        invFilterResults.glotFlowTimeAligned = [zeros(-invFilterResults.timeLag, 1); ...
                                                invFilterResults.glotFlowTimeAligned(1 : end + invFilterResults.timeLag)];
        invFilterResults.glotFlowDTimeAligned = [zeros(-invFilterResults.timeLag, 1); ...
                                                invFilterResults.glotFlowDTimeAligned(1 : end + invFilterResults.timeLag)];
    elseif invFilterResults.timeLag > 0
        invFilterResults.glotFlowTimeAligned = [invFilterResults.glotFlowTimeAligned(invFilterResults.timeLag + 1 : end); ...
                                                zeros(invFilterResults.timeLag, 1)];
        invFilterResults.glotFlowDTimeAligned = [invFilterResults.glotFlowDTimeAligned(invFilterResults.timeLag + 1 : end); ...
                                                zeros(invFilterResults.timeLag, 1)];
    end
end

