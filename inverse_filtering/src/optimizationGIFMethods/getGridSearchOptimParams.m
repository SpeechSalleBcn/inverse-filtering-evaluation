function solp = getGridSearchOptimParams(G,fitresult)
% solp = getGridSearchOptimParams(G,fitresult)
%
% Get the optimimum params of a gridsearch according to fitresult minimum

[minval, minidx] = min(fitresult, [], 'all');
solp.VTorder = G.VTs(minidx);
solp.GSorder = G.GSs(minidx);
solp.LipRad = G.LRs(minidx);
if strcmp(G.method, 'Original-IAIF') || strcmp(G.method, 'IOP-IAIF')
    solp.HPflag = G.HPs(minidx);
elseif strcmp(G.method,'QCP')
    solp.DQ = G.DQs(minidx);
    solp.PQ = G.PQs(minidx);
    solp.RQ = G.RQs(minidx);
    solp.STcompensation = G.STs(minidx);
end
