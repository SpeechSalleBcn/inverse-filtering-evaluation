function [solp] = optimization(speechData, timeMarks, C)
    %% OPTIMIZATION
    % In case that C.optim.method = "grid-search", this function calls to optimizations of different GIF
    % methods (optimizationGIF_Original_IAIF, optimizationGIF_IOP_IAIF, optimizationGIF_GFM_IAIF, 
    % optimizationGIF_QCP_par). If using calls to ndgrid.m, arrayfun.m and min.m to
    % evaluate all hyperparameter candidates and select the one with minum error.
    %
    % In case C.optim.method = "surrogate" calls to optimproblem.m,
    % optimoptions.m and solve.m with "surrogate" solver are invoked.
    switch C.optim.method
        case 'grid-search'
            G = getParamsGrid(C);
            G.method = C.gif.method;
            fitresult=zeros(size(G.VTs));
            switch C.gif.method
                case 'Original-IAIF'
                    tic
                    VTs=G.VTs; GSs=G.GSs; LRs=G.LRs; HPs=G.HPs;
                    parfor i=1:length(VTs)
                        fitresult(i) = optimizationGIF_Original_IAIF(VTs(i), GSs(i), LRs(i), HPs(i), speechData, timeMarks,  C);
                    end
                    toc
                case 'IOP-IAIF'
                    tic
                    VTs=G.VTs; GSs=G.GSs; LRs=G.LRs; HPs=G.HPs;
                    parfor i=1:length(VTs)
                        fitresult(i) = optimizationGIF_IOP_IAIF(VTs(i), GSs(i), LRs(i), HPs(i), speechData, timeMarks,  C);
                    end
                    toc
                case 'GFM-IAIF'
                    tic
                    VTs=G.VTs; GSs=G.GSs; LRs=G.LRs;
                    parfor i=1:length(VTs)
                        fitresult(i) = optimizationGIF_GFM_IAIF(VTs(i), GSs(i), LRs(i), speechData,  timeMarks, C);
                    end
                    toc
                case 'QCP'
                    tic
                    VTs=G.VTs;GSs=G.GSs;LRs=G.LRs;DQs=G.DQs;PQs=G.PQs;RQs=G.RQs;STs=G.STs;
                    parfor i=1:length(VTs)
                        fitresult(i)=optimizationGIF_QCP(VTs(i), GSs(i), LRs(i), DQs(i), PQs(i), RQs(i), STs(i), speechData, timeMarks,  C);
                    end
                    toc
            end
            solp = getGridSearchOptimParams(G,fitresult);

            resultsPath = [C.paths.mainPath filesep C.paths.experimentDir filesep  C.paths.resultsDir filesep 'mat'];
            if ~exist(resultsPath, 'dir')
               mkdir(resultsPath)
            end
            save([resultsPath filesep C.filesData.speechFileName '_'  C.gif.GIFMethod '.mat'], 'G', 'fitresult');
            
        case 'surrogate'
            % Get GIF optimization expression
            [errorF, x0] = getGIFoptimExpr(C, speechData, timeMarks);
            % Get GIF optimization problem
            prob = optimproblem("Objective",errorF);
            % Get GIF optimization options
            options = optimoptions("surrogateopt",...
                                "UseParallel", C.optim.UseParallel,...
                                'MaxFunctionEvaluations', C.optim.MaxFunctionEvaluations,...
                                "MinSurrogatePoints", 20,...
                                "PlotFcn",C.optim.PlotFcn...
                                );
            tic
            % GIF optimization solve
            [solp, ~, ~, ~] = solve(prob, x0, "Options", options, "Solver", "surrogateopt");
            toc
    end
end
