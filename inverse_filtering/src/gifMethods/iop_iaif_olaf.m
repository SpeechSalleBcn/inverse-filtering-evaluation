% IOP-IAIF Glottal Inverse Filtering
%   [G,Gd,a,b,ag] = iop_iaif_olaf(s,fs,p_vt,p_gl,d,tmarks,hpfilt)
%
% Description
%   This function estimates vocal tract linear prediction coefficients and
%   the glottal volume velocity waveform from a speech signal frame using
%   Iterative Optimal Pre-emphasis (IOP) Iterative Adaptive Inverse Filtering (IAIF) method.
%
%  New version that incorporates OLA processing (variable tmarks) (Joan
%  Claudi Socoró)
%
% Inputs
%   s      : Speech signal  [samples]
%   fs     : Sampling frequency [Hz]
%   p_vt   : Order of LPC/ARMA analysis for vocal tract
%   p_gl   : Order of LPC analysis for glottal source
%   d      : Leaky integration coefficient (e.g. 0.99)
%   hpfilt : High-pass filter flag (0: do not apply, 1...N: apply N times)
%  tmarks: matrix of Nf x 2 where (tmarks(n,1),tmarks(n,2)) are the
%       interval of samples of frame n for the OLA processing.
%
% Outputs
%   G      : Glottal volume velocity waveform
%   Gd    : Glottal volume velocity derivative waveform
%   a      : matrix of denominator coefficients of vocal tract of [Nf x p_vt] where
%       a(n,1:p_vt) are the LPC/ARMA denominator vocal tract coeficients of frame n
%   b      : matrix of  numerator coefficients of vocal tract of [Nf x z_vt] where
%       b(n,1:z_vt) are the ARMA numerator vocal tract coeficients of frame n
%   ag     : matrix LPC coefficients of source spectrum of [Nf x p_gl] where
%       ag(n,1:p_gl) are the LPC glottal source coeficients of frame n
%
% Notes
%   This function does not perform pitch synchronous analysis. This ensures
%   the robustness of the method regardless of the GCI estimation
%   performance.
%
% Example
%   Simplest, type g = iaif(x,fs) to estimate the glottal flow of a speech
%   frame.
%
% References
%  [1] P. Alku, "Glottal wave analysis with pitch synchronous iterative
%      adaptive inverse filtering", Speech Communication, vol. 11, no. 2-3,
%      pp. 109–118, 1992.
%
% Copyright (c) 2013 Aalto University
%
% License
%  This file is under the LGPL license,  you can
%  redistribute it and/or modify it under the terms of the GNU Lesser General 
%  Public License as published by the Free Software Foundation, either version 3 
%  of the License, or (at your option) any later version. This file is
%  distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
%  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
%  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
%  details.
%
% This function is part of the Common Speech Processing Repository 
% http://covarep.github.io/covarep/
%
% Octave Compatible
% 
% Author
%  Tuomo Raitio tuomo.raitio@aalto.fi
%
% $Id <info set by the versioning system> $
%
% Modifed version by
% Joan Claudi Socoró joanclaudi.socoro@salle.url.edu
%

function [G,Gd,a,ag,hpfilter_out,Niter] = iop_iaif_olaf(s,fs,p_vt,p_gl,d,hpfilt,tmarks,hpfilter_in)
    flagDebug = false;  % Only to show intermediate results of frame-based processing

    % Set default parameters
    if nargin < 8
        hpfilter_in = [];
        if nargin < 6
            hpfilt = 1;
            if nargin < 5
                d = 0.99;
                if nargin < 4
                    p_gl = 2*round(fs/4000);
                    if nargin < 3
                        p_vt = 2*round(fs/2000)+4;
                        if nargin < 2
                            disp('Error: Not enough input arguments.');
                        end
                    end
                end
            end
        end
    end
    preflt = p_vt+1;

    % Ensure column vector form
    s = s(:);

    % High-pass filter speech in order to remove possible low frequency
    % fluctuations (Linear-phase FIR, Fc = 70 Hz)
    hpfilter_out = [];
    if hpfilt > 0
        Fstop = 40;                 % Stopband Frequency
        Fpass = 70;                 % Passband Frequency
        Nfir = round(300/16000*fs); % FIR numerator order
        if mod(Nfir,2) == 1
            Nfir = Nfir + 1;
        end

        % it is very very expensive to calculate the firls filter! However, as 
        % long as the fs does not change, the firls filter does not change.
        % Therefore, the computed filter is returned and can be passed to this
        % function later on to avoid the calculated of the (same) filter.
        if ~isempty(hpfilter_in)
            B = hpfilter_in;
        else
            B = hpfilter_fir(Fstop,Fpass,fs,Nfir);
        end
        hpfilter_out = B;

        for i = 1:hpfilt
            x = filter(B,1,[s ; zeros(round(length(B)/2)-1,1)]);
            x = x(round(length(B)/2):end);
        end
    end

    % Estimate the combined effect of the glottal flow and the lip radiation
    % (Hg1) and cancel it out through inverse filtering. Note that before
    % filtering, a mean-normalized pre-frame ramp is appended in order to
    % diminish ripple in the beginning of the frame. The ramp is removed after
    % filtering.

    if length(s)<p_vt
        error('Signal is too short')
    end

    % Global Glotal flow signal
    G = zeros(size(s));
    Gd = G; % Derivative
    % Number of frames
    nf = size(tmarks,2);
    % Frame processing time marks
    tcenter = round(tmarks(1,:)+tmarks(2,:))/2;   % Center Window Marks
    tmc = round((tcenter(2:end)+tcenter(1:end-1))/2);  % Mid Points between Center window Marks
    tmc = [tmarks(1,1) tmc tmarks(2,end)];  % Complete for frame-based processing with filters

    % Initialization of output matrices
    a = cell(1,nf);
    ag = cell(1,nf);
    % Initialization of FIR filters memories (inverse filters)
    Z_Hg1 = 0;  % Memory of 1st inverse glottal filter
    Z_Hg2 = zeros(1,p_gl);  % Memory of 2nd inverse glottal filter
    Z_Hvt1 = zeros(1,p_vt); % Memory of 1st inverse vocal tract filter
    Z_Hvt2 = zeros(1,p_vt); % Memory of 2nd inverse vocal tract filter
    Z_int1 = 0; % Memory of 1st integrator
    Z_int2 = 0; % Memory of 2nd integrator
    Z_int3 = 0; % Memory of 3st integrator

    Z_pre = cell(1,50); % Memory of 3st integrator
    for i = 1:49
        Z_pre{i} = zeros(1,i);
    end

    for n = 1:nf

        % Signal frame for LPC computations
        x = s(tmarks(1,n):tmarks(2,n));
        % Signal frame for inverse filtering computations
        if n == 1
            signal = [linspace(-s(tmc(1)),s(tmc(1)),preflt)' ;s(tmc(n):(tmc(n+1)-1))];
            idx =  (preflt+1):length(signal);
        else
            signal = s(tmc(n):(tmc(n+1)-1));
            idx = 1:length(signal);
        end

        win = hanning(length(x));
        %signal = [linspace(-x(1),x(1),preflt)' ; x];
        %idx = preflt+1:numel(signal);

        % Traiem el tema pre-ramp ja que ara és important donar continuïtat
        % entre trames adjacents
        %signal = x;
        %idx = 1:numel(signal);
        %disp(['Frame n = ',num2str(n)])
        ag1_it = lpc(x.*hann(length(x)),1);         % First 1st order LPC estimation
        ag1_cum = ag1_it;
        i = 1;
         while (abs(ag1_it(2))>0.001)
            % Cancel current estimate of glottis contribution from speech signal
            [x_v1x,Z_pre{i}] = filter(ag1_cum,1,signal,Z_pre{i});	% Inverse filtering
            s_v1x = x_v1x(idx);      % Remove pre-ramp
        
            % Next 1st order LPC estimation
            ag1_it = lpc(s_v1x.*hann(length(s_v1x)),1);   % 1st order LPC
        
            % Update gross estimate of glottis contribution
            ag1_cum = conv(ag1_cum,ag1_it);	% Combine 1st order estimation with previous
            i = i + 1;
        end
        [y,Z_pre{i}] = filter(ag1_cum,1,signal,Z_pre{i});
        y = y(idx);

        Niter(n) = i;

        % Estimate the effect of the vocal tract (Hvt1) and cancel it out through
        % inverse filtering. The effect of the lip radiation is canceled through
        % integration. Signal g1 is the first estimate of the glottal flow.
        %Hvt1 = lpc(y.*hanning(length(y)),p_vt);
        [Avt1,Bvt1] = vtEstimation(y.*hanning(length(y)),p_vt);

        [g1,Z_Hvt1] = filter(Avt1,Bvt1,signal,Z_Hvt1);
        [g1,Z_int1] = filter(1,[1 -d],g1,Z_int1);
        g1 = g1(idx);

        % Re-estimate the effect of the glottal flow (Hg2). Cancel the contribution
        % of the glottis and the lip radiation through inverse filtering and
        % integration, respectively.
        Hg2 = lpc(g1.*hanning(length(g1)),p_gl);
        [y,Z_Hg2] = filter(Hg2,1,signal,Z_Hg2);
        [y,Z_int2] = filter(1,[1 -d],y,Z_int2);
        y = y(idx);

        % Estimate the model for the vocal tract (Hvt2) and cancel it out through
        % inverse filtering. The final estimate of the glottal flow is obtained
        % through canceling the effect of the lip radiation.
        %Hvt2 = lpc(y.*hanning(length(y)),p_vt);
        [Avt2,Bvt2] = vtEstimation(y.*hanning(length(y)),p_vt);

        [dg,Z_Hvt2] = filter(Avt2,Bvt2,signal,Z_Hvt2);
        [g,Z_int3] = filter(1,[1 -d],dg,Z_int3);
        g = g(idx);
        dg = dg(idx);

        G(tmc(n):(tmc(n+1)-1)) = g;
        Gd(tmc(n):(tmc(n+1)-1)) = dg;

        % Set vocal tract model to 'a' and glottal source spectral model to 'ag'
        a{n} = Avt2;
        ag{n} = Hg2;

        if flagDebug
            figure(1);clf
            subplot(311);
            plot(s);hold on;plot(tmarks(1,n):tmarks(2,n),s(tmarks(1,n):tmarks(2,n)),'r'); hold off; axis tight;grid;
            title('LPC frame')
            subplot(312);
            plot(s);hold on;plot(tmc(n):(tmc(n+1)-1),s(tmc(n):(tmc(n+1)-1)),'r'); hold off; axis tight;grid;
            title('Filter frame')
            subplot(313);
            plot(G);axis tight;grid;
            title('Output')
            pause
        end
    end
    Niter = median(Niter);
end

function B = hpfilter_fir(Fstop,Fpass,fs,N)
% FIR least-squares Highpass filter design using the FIRLS function
% Tuomo Raitio
% 20.8.2013
% Calculate the coefficients using the FIRLS function.
    B  = firls(N, [0 Fstop Fpass fs/2]/(fs/2), [0 0 1 1], [1 1]);
end
