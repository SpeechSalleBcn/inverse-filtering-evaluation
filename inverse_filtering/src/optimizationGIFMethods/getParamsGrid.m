function G = getParamsGrid(C)
% G = getParamsGrid(C)
%
% Build the grid of parameters G according to the configuration C.

VTorder = parseArrayParameter(C.gif.VTorder);
GSorder = parseArrayParameter(C.gif.GSorder);
LipRad = parseArrayParameter(C.gif.LipRad);
method = C.gif.method;
if strcmp(method, 'Original-IAIF') || strcmp(method, 'IOP-IAIF')
    HPflag = parseArrayParameter(C.gif.HPflag);
    [VTs, GSs, LRs, HPs] = ndgrid(VTorder, GSorder, LipRad, HPflag);
    G.VTs=VTs(:);G.GSs=GSs(:);G.LRs=LRs(:);G.HPs=HPs(:);
elseif strcmp(method,'GFM-IAIF')
    [VTs, GSs, LRs] = ndgrid(VTorder, GSorder, LipRad);
    G.VTs=VTs(:);G.GSs=GSs(:);G.LRs=LRs(:);
elseif strcmp(method, 'QCP')
    DQ = parseArrayParameter(C.gif.DQ);
    PQ = parseArrayParameter(C.gif.PQ);
    RQ = parseArrayParameter(C.gif.RQ);
    STcompensation = parseArrayParameter(C.gif.STcompensation);
    [VTs,GSs,LRs,DQs,PQs,RQs,STs] = ndgrid(VTorder, GSorder,LipRad,DQ,PQ,RQ,STcompensation);
    G.VTs=VTs(:);G.GSs=GSs(:);G.LRs=LRs(:);G.DQs=DQs(:);G.PQs=PQs(:);G.RQs=RQs(:);G.STs=STs(:);
end
end

function out = parseArrayParameter(in)
    if length(in) == 3 % in corresponds to initial_value: increment: final_value
        if in(2)==0; in(2)=1; end % Matlab cannot handle increment 0
        out = in(1) : in(2) : in(3);
    elseif length(in)==1 || length(in)>3 % in directly contains the value or the array of values
        out = in;
    end
end