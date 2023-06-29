function [invFilterResults] = normalizationGlotFlow(invFilterResults, speechData, timeMarks, C)
% function [g_opt,gd_opt,a_gif,b_gif,err,errNoOpt] = optimizeg2(g,gd,g_GT,gcis,opt)
% g_opt = optimizeg2(g_opt,gd_opt,a,b)

% Performs the normalization of glottal waveform input invFilterResults.glotFlow based on applying
% both scaling and DC removal of each signal period so as to resemble the
% reference ground-truth (GT) glottal waveform speechData.glotFlowGT. 

% For each signal period, based on input timeMarks.pitchMarks (Glottal Closure Instants), the optimum scale and 
% DC coefficients are computed based on RMSE between the reference period
% and the modified signal.
% Normalized glottal wavform and its derivative are returned as glotFlowInvFiltOpt and
% glotFlowDInvFiltOpt, respectively, and the scaling and DC parameters for each signal
% period are returned in a_gif and b_gif, respectively.

% parameters are a parametric linear regression between input glottal waveform invFilterResults.glotFlow
% and the ground-truth glottal one, in order to perform a pairwise
% comparison 
% Input timeMarks.pitchMarks are used to perform the linear regressions

% follwing: glotFlowInvFiltOpt(n) = a·invFilterResults.glotFlow(n)+b, so as to minimize euclidean distance with
% ground truth E = sum(abs(glotFlowInvFiltOpt(n)-speechData.glotFlowGT(n))).

if length(invFilterResults.glotFlow) ~= length(speechData.glotFlowGT)
    error('normalizationGlotFlow: Input vectors must have the same length');
end
% column vectors % TODO: centralize column vector check
if ~iscolumn(invFilterResults.glotFlowTimeAligned); invFilterResults.glotFlowTimeAligned = transpose(invFilterResults.glotFlowTimeAligned); end
if ~iscolumn(speechData.glotFlowGT); speechData.glotFlowGT = transpose(speechData.glotFlowGT); end

invFilterResults.scale = zeros(size(timeMarks.pitchMarks, 2), 1);
invFilterResults.dcOffset = zeros(size(timeMarks.pitchMarks, 2), 1);
for n = 2 : size(timeMarks.pitchMarks, 2) - 1
    g_frame = invFilterResults.glotFlowTimeAligned(timeMarks.pitchMarks(1, n) : timeMarks.pitchMarks(2, n));
    g_GT_frame = speechData.glotFlowGT(timeMarks.pitchMarks(1, n) : timeMarks.pitchMarks(2, n));
    MatrixSys = [sum(abs(g_frame) .^2) sum(g_frame); sum(g_frame) length(g_frame)];
    VectorSys = [sum(g_frame .* g_GT_frame); sum(g_GT_frame)];
    Solution = MatrixSys \ VectorSys;
    if Solution(1) < 0
        Solution(1) = 0;
    end
    invFilterResults.scale(n) = Solution(1);
    invFilterResults.dcOffset(n) = Solution(2);
end
