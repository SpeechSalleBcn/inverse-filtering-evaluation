classdef GIFErrors < handle
    properties
        invFilterResults
        speechData
        timeMarks
        config
        rmse
        rmsePulse
        rmsePulseNorm
        % rmsePulseNormInterp
        naqInvFilter
        naqGroundTruth
        qoqInvFilter
        qoqGroundTruth
        h1h2InvFilter
        h1h2GroundTruth
        hrfInvFilter
        hrfGroundTruth
        pspInvFilter
        pspGroundTruth
        spectralTiltGroundTruth
        spectralTiltInvFilter
        peakSlopeGroundTruth
        peakSlopeInvFilter
        spectralDist
        spectralDistComp
    end
    methods
        function obj = GIFErrors(invFilterResults, speechData, timeMarks, C)
            obj.timeMarks = timeMarks;
            obj.config = C;
            obj.speechData = speechData;
            obj.speechData.glotFlowGTD = filter([1 -C.gif.derivLipRad], 1, speechData.glotFlowGT); % Derivative of Glottal Flow Ground Truth
            obj.invFilterResults = invFilterResults;
            obj.conditionData();
        end

        function conditionData(obj)
            obj.invFilterResults = timeAlignmentGlotFlow(obj.invFilterResults, obj.speechData, obj.config);
            obj.invFilterResults = normalizationGlotFlow(obj.invFilterResults, obj.speechData, obj.timeMarks, obj.config);
            if ~iscolumn(obj.speechData.glotFlowGT); obj.speechData.glotFlowGT = transpose(obj.speechData.glotFlowGT); end
            if ~iscolumn(obj.speechData.glotFlowGTD); obj.speechData.glotFlowGTD = transpose(obj.speechData.glotFlowGTD); end
            if ~iscolumn(obj.invFilterResults.glotFlow)
                obj.invFilterResults.glotFlow = transpose(obj.invFilterResults.glotFlow);
            end
            if ~iscolumn(obj.invFilterResults.glotFlowD)
                obj.invFilterResults.glotFlowD = transpose(obj.invFilterResults.glotFlowD);
            end
            if ~iscolumn(obj.invFilterResults.glotFlowTimeAligned)
                obj.invFilterResults.glotFlowTimeAligned = transpose(obj.invFilterResults.glotFlowTimeAligned);
            end
        end

        function computeAllErrors(obj)
            obj.computeVqParamsGroundTruth();
            obj.computeVqParamsInvFilter();
            obj.computePeakSlope();
            obj.computeSpectralDistortion();
            obj.computeRMSEPulse();
            obj.computeRMSE();
        end
        
        function computeVqParamsGroundTruth(obj)
            [NAQ, QOQ, H1H2, HRF, PSP, ST] = get_vq_params(...
                                        obj.speechData.glotFlowGT,...
                                        obj.speechData.glotFlowGTD,...
                                        obj.speechData.fs,...
                                        obj.timeMarks.gcisTimes...
                                        );
            obj.naqGroundTruth = NAQ(obj.timeMarks.isGCISelectedArray, 2);
            obj.qoqGroundTruth = QOQ(obj.timeMarks.isGCISelectedArray, 2);
            obj.h1h2GroundTruth = H1H2(obj.timeMarks.isGCISelectedArray, 2);
            obj.hrfGroundTruth = HRF(obj.timeMarks.isGCISelectedArray, 2);
            obj.pspGroundTruth = PSP(obj.timeMarks.isGCISelectedArray, 2);
            obj.spectralTiltGroundTruth = ST(obj.timeMarks.isGCISelectedArray, 2);
        end
        
        function computeVqParamsInvFilter(obj)
            [NAQ, QOQ, H1H2, HRF, PSP, ST] = get_vq_params(...
                                        obj.invFilterResults.glotFlow,...
                                        obj.invFilterResults.glotFlowD,...
                                        obj.speechData.fs,...
                                        obj.timeMarks.gcisTimes...
                                        );
            obj.naqInvFilter = NAQ(obj.timeMarks.isGCISelectedArray, 2);
            obj.qoqInvFilter = QOQ(obj.timeMarks.isGCISelectedArray, 2);
            obj.h1h2InvFilter = H1H2(obj.timeMarks.isGCISelectedArray, 2);
            obj.hrfInvFilter = HRF(obj.timeMarks.isGCISelectedArray, 2);
            obj.pspInvFilter = PSP(obj.timeMarks.isGCISelectedArray, 2);
            obj.spectralTiltInvFilter = ST(obj.timeMarks.isGCISelectedArray, 2);
        end

        function computePeakSlope(obj)
            peakSlope = peakslope(obj.speechData.glotFlowGT, obj.speechData.fs);
            peakSlopeInterpolated = interp1td(peakSlope, obj.timeMarks.gcisTimes); % interpolate ST values at pulse level
            obj.peakSlopeGroundTruth = peakSlopeInterpolated(obj.timeMarks.isGCISelectedArray, 2);
            peakSlope = peakslope(obj.invFilterResults.glotFlow, obj.speechData.fs);
            obj.peakSlopeInvFilter = peakSlopeInterpolated(obj.timeMarks.isGCISelectedArray, 2);
        end

        function computeRMSEPulse(obj)
            err = zeros(size(obj.timeMarks.pitchMarks,2),1);
            errNorm = zeros(size(obj.timeMarks.pitchMarks,2),1);
            errNormInterp = zeros(size(obj.timeMarks.pitchMarks,2),1);
            % [scaleInterpolated, dcOffSetInterpolated] = interpolateScaleAndDCOffset(obj.invFilterResults, obj.timeMarks);
            % glotFlowIFNormInterp = obj.invFilterResults.glotFlowTimeAligned .* scaleInterpolated + dcOffSetInterpolated;
            for n = 2:size(obj.timeMarks.pitchMarks, 2)-1
                glotFlowIF_frame = obj.invFilterResults.glotFlow(obj.timeMarks.pitchMarks(1,n):obj.timeMarks.pitchMarks(2,n));
                glotFlowIF_frameTimeAligned = obj.invFilterResults.glotFlowTimeAligned(obj.timeMarks.pitchMarks(1,n):obj.timeMarks.pitchMarks(2,n));
                glotFlowIF_frameNorm = glotFlowIF_frameTimeAligned * obj.invFilterResults.scale(n) + obj.invFilterResults.dcOffset(n);
                % glotFlowIF_frameNormInterp = glotFlowIFNormInterp(obj.timeMarks.pitchMarks(1,n):obj.timeMarks.pitchMarks(2,n));
                glotFlowGT_frame = obj.speechData.glotFlowGT(obj.timeMarks.pitchMarks(1,n):obj.timeMarks.pitchMarks(2,n));
                if (obj.config.error.relativeErrorRMSE)
                    err(n) = 100 * norm(glotFlowGT_frame - glotFlowIF_frameTimeAligned) / norm(glotFlowGT_frame);
                    errNorm(n) = 100 * norm(glotFlowGT_frame - glotFlowIF_frameNorm) / norm(glotFlowGT_frame);
                    % errNormInterp(n) = 100 * norm(glotFlowGT_frame - glotFlowIF_frameNormInterp) / norm(glotFlowGT_frame);
                else
                    err(n) = norm(glotFlowGT_frame -  glotFlowIF_frame);
                    errNorm(n) = norm(glotFlowGT_frame -  glotFlowIF_frameNorm);
                    % errNormInterp(n) = norm(glotFlowGT_frame -  glotFlowIF_frameNormInterp);
                end
            end
            err(1)=err(2);
            err(n+1)=err(n);
            obj.rmsePulse = err;
            errNorm(1)=errNorm(2);
            errNorm(n+1)=errNorm(n);
            obj.rmsePulseNorm = errNorm;
            % errNormInterp(1) = errNormInterp(2);
            % errNormInterp(n+1) = errNormInterp(n);
            % obj.rmsePulseNormInterp = errNormInterp;
        end

        function computeRMSE(obj)
            if (obj.config.error.relativeErrorRMSE)
                obj.rmse = 100 * norm(obj.speechData.glotFlowGT(obj.timeMarks.selectedSamplesIdxs) ...
                                       - obj.invFilterResults.glotFlow(obj.timeMarks.selectedSamplesIdxs)) ...
                                / norm(obj.speechData.glotFlowGT(obj.timeMarks.selectedSamplesIdxs));
            else
                obj.rmse = norm(obj.speechData.glotFlowGT(obj.timeMarks.selectedSamplesIdxs) ...
                                       - obj.invFilterResults.glotFlow(obj.timeMarks.selectedSamplesIdxs));
            end
        end

        % [Drugman2012]
        function [] = computeSpectralDistortion(obj)
            [obj.spectralDist, obj.spectralDistComp] = get_sd(...
                                                    obj.speechData.glotFlowGT,...
                                                    obj.invFilterResults.glotFlow,...
                                                    obj.speechData.fs,...
                                                    obj.timeMarks.gcisTimes...
                                                    );
        end

        function [tableFile, tablePulse] = getAllErrorsTables(obj)
            tableFile = table;
            tablePulse = table;
            errorRMSE = obj.getRMSE();
            errorRMSEPulse = obj.getRMSEPulse();
            errorNAQ = obj.getErrorNAQ();
            errorQOQ = obj.getErrorQOQ();
            errorH1H2 = obj.getErrorH1H2();
            errorHRF = obj.getErrorHRF();
            errorHRFdB = obj.getErrorHRFdB();
            errorPSP = obj.getErrorPSP();
            spectralTilt = obj.getSpectralTilt();
            errorSpectralTilt = obj.getErrorSpectralTilt();
            errorPeakSlope = obj.getErrorPeakSlope();
            errorSpectralDist = obj.getErrorSpectralDistortion();
            tableFile = [
                        struct2table(errorRMSE.file),...
                        struct2table(errorRMSEPulse.file),...
                        struct2table(errorNAQ.file),...
                        struct2table(errorQOQ.file),...
                        struct2table(errorH1H2.file),...
                        struct2table(errorHRF.file),...
                        struct2table(errorHRFdB.file),...
                        struct2table(errorPSP.file),...
                        struct2table(spectralTilt.file),...
                        struct2table(errorSpectralTilt.file),...
                        struct2table(errorPeakSlope.file),...
                        struct2table(errorSpectralDist.file),...
                        ];
            tablePulse = [
                        struct2table(errorRMSEPulse.pulse),...
                        struct2table(errorNAQ.pulse),...
                        struct2table(errorQOQ.pulse),...
                        struct2table(errorH1H2.pulse),...
                        struct2table(errorHRF.pulse),...
                        struct2table(errorHRFdB.pulse),...
                        struct2table(errorPSP.pulse),...
                        struct2table(spectralTilt.pulse),...
                        struct2table(errorSpectralTilt.pulse),...
                        struct2table(errorPeakSlope.pulse),...
                        struct2table(errorSpectralDist.pulse),...
                        ];
        end

        function [normParams] = getNormParams(obj)
            normParams.scale = obj.invFilterResults.scale;
            normParams.dcOffset = obj.invFilterResults.dcOffset;
            normParams.timeLag = obj.invFilterResults.timeLag;
        end

        function [error] = getRMSEPulse(obj)
            error.file.rmseMedian = median(obj.rmsePulse(obj.timeMarks.isGCISelectedArray), 'omitnan');
            error.file.rmseNormMedian = median(obj.rmsePulseNorm(obj.timeMarks.isGCISelectedArray), 'omitnan');
            % error.file.rmseNormInterpMedian = median(obj.rmsePulseNormInterp(obj.timeMarks.isGCISelectedArray), 'omitnan');
            error.pulse.rmse = obj.rmsePulse(obj.timeMarks.isGCISelectedArray);
            error.pulse.rmseNorm = obj.rmsePulseNorm(obj.timeMarks.isGCISelectedArray);
            % error.pulse.rmseNormInterp = obj.rmsePulseNormInterp(obj.timeMarks.isGCISelectedArray);
        end
                errNormInterp
        function [error] = getRMSE(obj)
            error.file.rmse = obj.rmse;
        end

        function [error] = getErrorNAQ(obj)
            error.pulse.naq = abs(obj.naqGroundTruth - obj.naqInvFilter)./obj.naqGroundTruth;
            filteredNAQ = error.pulse.naq(~isinf(error.pulse.naq));
            error.file.naqInf = 0;
            if sum(isinf(error.pulse.naq)) > 0
                error.file.naqInf = sum(isinf(error.pulse.naq)) / length(error.pulse.naq);
            end
            error.pulse.naqSign = (obj.naqGroundTruth - obj.naqInvFilter)./obj.naqGroundTruth;
            filteredNAQSign = error.pulse.naqSign(~isinf(error.pulse.naqSign));
            error.file.naqMean = mean(filteredNAQ, 'omitnan');
            error.file.naqMedian = median(filteredNAQ, 'omitnan');
            error.file.naqMedianSign = median(filteredNAQSign, 'omitnan');
        end

        function [error] = getErrorQOQ(obj)
            error.pulse.qoq = abs(obj.qoqGroundTruth - obj.qoqInvFilter)./obj.qoqGroundTruth;
            filteredQOQ = error.pulse.qoq(~isinf(error.pulse.qoq));
            error.file.qoqInf = 0;
            if sum(isinf(error.pulse.qoq)) > 0
                error.file.qoqInf = sum(isinf(error.pulse.qoq)) / length(error.pulse.qoq);
            end
            error.pulse.qoqSign = (obj.qoqGroundTruth - obj.qoqInvFilter)./obj.qoqGroundTruth;
            filteredQOQSign = error.pulse.qoqSign(~isinf(error.pulse.qoqSign));
            error.file.qoqMean = mean(filteredQOQ, 'omitnan');
            error.file.qoqMedian = median(filteredQOQ, 'omitnan');
            error.file.qoqMedianSign = median(filteredQOQSign, 'omitnan');
        end

        function [error] = getErrorH1H2(obj)
            error.pulse.h1h2 = abs(obj.h1h2GroundTruth - obj.h1h2InvFilter);
            error.pulse.h1h2Sign = (obj.h1h2GroundTruth - obj.h1h2InvFilter);
            error.file.h1h2Mean = mean(abs(obj.h1h2GroundTruth - obj.h1h2InvFilter));
            error.file.h1h2Median = median(abs(obj.h1h2GroundTruth - obj.h1h2InvFilter));
            error.file.h1h2MedianSign = median((obj.h1h2GroundTruth - obj.h1h2InvFilter));
        end

        function [error] = getErrorHRF(obj)
            error.pulse.hrf = abs((obj.hrfGroundTruth - obj.hrfInvFilter)./obj.hrfGroundTruth);
            error.pulse.hrfSign = ((obj.hrfGroundTruth - obj.hrfInvFilter)./obj.hrfGroundTruth);
            error.file.hrfMean = mean(abs((obj.hrfGroundTruth - obj.hrfInvFilter)./obj.hrfGroundTruth), 'omitnan');
            error.file.hrfMedian = median(abs((obj.hrfGroundTruth - obj.hrfInvFilter)./obj.hrfGroundTruth), 'omitnan');
            error.file.hrfMedianSign = median(((obj.hrfGroundTruth - obj.hrfInvFilter)./obj.hrfGroundTruth), 'omitnan');
        end

        function [error] = getErrorHRFdB(obj)
            error.pulse.hrfdB = abs(obj.hrfGroundTruth - obj.hrfInvFilter);
            error.pulse.hrfdBSign = (obj.hrfGroundTruth - obj.hrfInvFilter);
            error.file.hrfMeandB = mean(abs(obj.hrfGroundTruth - obj.hrfInvFilter), 'omitnan');
            error.file.hrfMediandB = median(abs(obj.hrfGroundTruth - obj.hrfInvFilter), 'omitnan');
            error.file.hrfMedianSigndB = median((obj.hrfGroundTruth - obj.hrfInvFilter), 'omitnan');
        end
        
        function [error] = getErrorPSP(obj)
            error.pulse.psp = abs((obj.pspGroundTruth - obj.pspInvFilter)./obj.pspGroundTruth);
            error.pulse.pspSign = ((obj.pspGroundTruth - obj.pspInvFilter)./obj.pspGroundTruth);
            error.file.pspMean = mean(abs((obj.pspGroundTruth - obj.pspInvFilter)./obj.pspGroundTruth), 'omitnan');
            error.file.pspMedian = median(abs((obj.pspGroundTruth - obj.pspInvFilter)./obj.pspGroundTruth), 'omitnan');
            error.file.pspMedianSign = median(((obj.pspGroundTruth - obj.pspInvFilter)./obj.pspGroundTruth), 'omitnan');
        end

        function [spectralTilt] = getSpectralTilt(obj)
            spectralTilt.file.spectralTiltInvFilter = median(obj.spectralTiltInvFilter, 'omitnan');;
            spectralTilt.file.spectralTiltGroundTruth = median(obj.spectralTiltGroundTruth, 'omitnan');;
            spectralTilt.pulse.spectralTiltInvFilter = obj.spectralTiltInvFilter;
            spectralTilt.pulse.spectralTiltGroundTruth = obj.spectralTiltGroundTruth;
        end

        function [error] = getErrorSpectralTilt(obj)
            error.pulse.errorSpectralTilt = abs(obj.spectralTiltGroundTruth - obj.spectralTiltInvFilter);
            error.pulse.errorSpectralTiltSign = (obj.spectralTiltGroundTruth - obj.spectralTiltInvFilter);
            error.file.errorSpectralTiltMean = mean(abs(obj.spectralTiltGroundTruth - obj.spectralTiltInvFilter));
            error.file.errorSpectralTiltMedian = median(abs(obj.spectralTiltGroundTruth - obj.spectralTiltInvFilter));
            error.file.errorSpectralTiltMedianSign = median((obj.spectralTiltGroundTruth - obj.spectralTiltInvFilter));
        end

        function [spectralTilt] = getPeakSlope(obj)
            spectralTilt.file.peakSlopeInvFilter = median(obj.peakSlopeInvFilter, 'omitnan');;
            spectralTilt.file.peakSlopeGroundTruth = median(obj.peakSlopeGroundTruth, 'omitnan');;
            spectralTilt.pulse.peakSlopeInvFilter = obj.peakSlopeInvFilter;
            spectralTilt.pulse.peakSlopeGroundTruth = obj.peakSlopeGroundTruth;
        end

        function [error] = getErrorPeakSlope(obj)
            error.pulse.errorPeakSlope = abs(obj.peakSlopeGroundTruth - obj.peakSlopeInvFilter);
            error.pulse.errorPeakSlopeSign = (obj.peakSlopeGroundTruth - obj.peakSlopeInvFilter);
            error.file.errorPeakSlopeMean = mean(abs(obj.peakSlopeGroundTruth - obj.peakSlopeInvFilter));
            error.file.errorPeakSlopeMedian = median(abs(obj.peakSlopeGroundTruth - obj.peakSlopeInvFilter));
            error.file.errorPeakSlopeMedianSign = median((obj.peakSlopeGroundTruth - obj.peakSlopeInvFilter));
        end

        function [error] = getErrorSpectralDistortion(obj)
            error.pulse.spectralDist = obj.spectralDist(obj.timeMarks.isGCISelectedArray, 2);
            error.pulse.spectralDistComp = obj.spectralDistComp(obj.timeMarks.isGCISelectedArray, 2);
            error.file.spectralDistMedian = median(obj.spectralDist(obj.timeMarks.isGCISelectedArray, 2), 'omitnan');
            error.file.spectralDistCompMedian = median(obj.spectralDistComp(obj.timeMarks.isGCISelectedArray, 2), 'omitnan');
        end
    end
end