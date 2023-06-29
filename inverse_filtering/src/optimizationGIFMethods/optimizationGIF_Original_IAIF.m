function error = optimizationGIF_Original_IAIF(VTorder, GSorder, LipRad, HPflag, speechData, timeMarks, C)
%        error = optimizationGIF_Original_IAIF(VTorder, GSorder, LipRad, HPflag, speechData, timeMarks, C)
%
% Function that allows evaluation of GIF technique IAIF in order to provide
% its optimization via arrayfun.m (grid-search) and solve.m (surrogate) methods. Their 
% corresponding hyperparameters are passed as input arguments:
%
%   - VTorder: Vocal-Tract order (e.g. 12)
%   - GSorder: Glottal-Source order (e.g. 3)
%   - LipRad: lip-radiation coefficient (e.g. 0.99)
%   - HPflag: High-Pass filter flag (e.g. 0 or 1)
%
% As output (error) a given metric is returned.
%

    %% Inverse Filtering
    [invFilterResults.glotFlow, invFilterResults.glotFlowD] = iaif_olaf(speechData.speech,...
                                                                      speechData.fs,...
                                                                      VTorder,...
                                                                      GSorder,...
                                                                      LipRad,...
                                                                      HPflag,...
                                                                      timeMarks.analysisFramesSamples...
                                                                     );

    errors = GIFErrors(invFilterResults, speechData, timeMarks, C);
    if (C.error.pulseLevel)
        errors.computeRMSEPulse();
        errorRMSEPulse = errors.getRMSEPulse();
        error = errorRMSEPulse.file.rmseNormMedian;
    else
        errors.computeRMSE();
        errorRMSE = errors.getRMSE();
        error = errorRMSE.file.rmse;
    end
end
