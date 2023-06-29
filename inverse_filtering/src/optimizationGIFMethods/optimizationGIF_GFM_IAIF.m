function error = optimizationGIF_GFM_IAIF(VTorder, GSorder, LipRad, speechData, timeMarks, C)
%error = optimizationGIF_GFM_IAIF(VTorder,GSorder,LipRad)
%
% Function that allows evaluation of GIF technique Glottal-Flow Model (GFM) in order to provide
% its optimization via arrayfun.m (grid-search) and solve.m (surrogate) methods. Their 
% corresponding hyperparameters are passed as input arguments:
%
%   - VTorder: Vocal-Tract order (e.g. 12)
%   - GSorder: Glottal-Source order (e.g. 3)
%   - LipRad: lip-radiation coefficient (e.g. 0.99)

    %% Inverse Filtering
    [invFilterResults.glotFlow, invFilterResults.glotFlowD] = gfmiaif_olaf(speechData.speech,...
                                                                      speechData.fs,...
                                                                      VTorder,...
                                                                      GSorder,...
                                                                      LipRad,...
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

