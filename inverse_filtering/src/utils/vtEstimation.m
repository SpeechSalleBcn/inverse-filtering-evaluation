function [A,B] = vtEstimation(x,p)
% [A,B] = vtEstimation(x,p)
% Vocal tract estimation through LPC (all pole modelling with p poles) 
    B = 1;
    A = lpc(x,p);
end