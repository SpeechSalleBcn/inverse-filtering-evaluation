function error = optimizationGIF_QCP(VTorder, GSorder, LipRad, DQ, PQ, RQ, STcompensation, speechData, timeMarks, C)
%error = optimizationGIF_QCP(VTorder, GSorder, LipRad, DQ, PQ, RQ, STcompensation)
%
% Function that allows evaluation of GIF technique Quasi-Closed Phase (QCP) in order to provide
% its optimization via arrayfun.m (grid-search) and solve.m (surrogate) methods. Their 
% corresponding hyperparameters are passed as input arguments:
%
%   - VTorder: Vocal-Tract order (e.g. 12)
%   - GSorder: Glottal-Source order (e.g. 3)
%   - LipRad: lip-radiation coefficient (e.g. 0.99)
%   - DQ: duration quotient (e.g. 0.4)
%   - PQ: position quotient (e.g. 0.1)
%   - RQ: ramp quotient (e.g. 0.1)
%   - STcompensation: spectral tilt compensation flag (e.g. 1 or 0)


    debug = false;

    %% Inverse Filtering
    [invFilterResults.glotFlow, invFilterResults.glotFlowD]  = qcp_olaf(speechData.speech,...
                                                                      speechData.fs,...
                                                                      VTorder,...
                                                                      GSorder,...
                                                                      LipRad,...
                                                                      DQ,...
                                                                      PQ,...
                                                                      RQ,...
                                                                      'causal',...
                                                                      STcompensation,...
                                                                      timeMarks.analysisFramesSamples,...
                                                                      timeMarks.gcisSamples...
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

    if debug
        GIFparams.VTorder = VTorder;
        GIFparams.GSorder = GSorder;
        GIFparams.LipRad = LipRad;
        GIFparams.DQ = DQ;
        GIFparams.PQ = PQ;
        GIFparams.RQ = RQ;
        GIFparams.STcompensation = STcompensation;
        plotGIFResults(speechData, invFilterResults.glotFlow, C, GIFparams, 'Optimization')
    end
end
