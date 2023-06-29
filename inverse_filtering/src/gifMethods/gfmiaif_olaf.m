% GFM-IAIF
% Glottal Flow Model-based Iterative Adaptive Inverse Filtering
%
% Description
%   This function estimates the linear prediction coefficients of both
%   vocal tract and glottis filters from a speech signal frame with the
%   GFM-IAIF method [1].
%   The latter is an extension of IAIF [2], with an improved pre-emphasis
%   step, that allows to extract a wide-band glottis response,
%   incorporating both glottal formant and spectral tilt characteristics.
%   This function is based on the iaif.m implementation from the COVAREP
%   toolbox [3].
%
%  New version that incorporates OLA processing (variable tmarks) (Joan
%  Claudi Socoró)
%
% Inputs
% (s) :  [Nx1]	Speech signal 
% (nv)	:  [1x1]	Number of poles of LP analysis for vocal tract (def. 48)
% (ng)	:  [1x1]	Order of LP analysis for glottal source (def. 3)
% (d)  	:  [1x1]	Leaky integration coefficient (def. 0.99)
% (tmarks): [Nf x 2] matrix of Nf x 2 where (tmarks(n,1),tmarks(n,2)) are the
%       interval of samples of frame n for the OLA processing
% (win)	:  [Nx1]	Window used before LPC (def. Hanning)
%
% Outputs
%  Av  	:  [nfxnv]	 Av(n,1:nv) are the denominator coefficients of vocal tract
%               contribution of frame n
%  Ag 	:  [nfxng]	ag(n,1:ng) are the LP coefficients of glottis
%               contribution of frame n
%  al 	:  [1x2]	LP coefficients of lip radiation contribution
%
%
% Examples
%  [av,ag,al] = gfmiaif(x) provides the LP coefficients of vocal tract,
%               glottis and lip radiation with default parameters
%  [av,ag,al] = gfmiaif(x,nv,ng,d,win) allows to choose parameters
%
% GFM-IAIF has been designed on the assumption that a third order filter
% allows to describe most of the glottis-related timbre variations (e.g.,
% tenseness, effort) with a compact set of parameters.
% Thus, the use of ng = 3 is highly encouraged.
%
%
% References
%  [1] O. Perrotin and I. V. McLoughlin (2019)
%      "A spectral glottal flow model for source-filter separation of
%      speech", in IEEE International Conference on Acoustics, Speech, and
%      Signal Processing (ICASSP), Brighton, UK, May 12-17, pp. 7160-7164.
%
%  [2] P. Alku (1992)
%      "Glottal wave analysis with pitch synchronous iterative adaptive
%      inverse filtering", Speech Communication, 11(2-3), pp. 109-118.
%
%  [3] G. Degottex, J. Kane, T. Drugman, T. Raitio and S. Scherer (2014)
%      "COVAREP - A collaborative voice analysis repository for speech
%      technologies", in IEEE International Conference on Acoustics,
%      Speech and Signal Processing (ICASSP), Florence, Italy, May 4-9,
%      pp. 960-964.
%
% How to cite
%   Cite reference [1] above
%
%
% Copyright (c) 2019 Univ. Grenoble Alpes, CNRS, Grenoble INP, GIPSA-lab
%
% License
%   This file is free software; you can redistribute it and/or modify it
%   under the terms of the GNU Lesser General Public License as published
%   by the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%   This file is distributed in the hope that it will be useful, but
%   WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
%   General Public License for more details. 
%
% Author
%  Olivier Perrotin olivier.perrotin@gipsa-lab.grenoble-inp.fr
%
% Modifed version by
% Joan Claudi Socoró joanclaudi.socoro@salle.url.edu
%

function [G,Gd,Av,Ag,al] = gfmiaif_olaf(s,fs,nv,ng,d,tmarks,win)

    % ----- Set default parameters -------------------------------------------

    if nargin < 7
        % Window for LPC estimation
        win = hann(length(s));
        if nargin < 5
            % Lip radiation leaky integration coefficient
            d = 0.99;
            if nargin < 4
                % Glottis LPC order
                ng = 3;
                if nargin < 3
                    % Vocal tract LPC order
                    nv = 48;
                end
            end
        end
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
    Z_ag1 = zeros(1,ng);  % Memory of 1st inverse glottal filter
    Z_ag = zeros(1,ng);  % Memory of 2nd inverse glottal filter
    Z_av1 = zeros(1,nv); % Memory of 1st inverse vocal tract filter
    Z_av = zeros(1,nv); % Memory of 2nd inverse vocal tract filter
    Z_int1 = 0; % Memory of 1st integrator
    Z_int2 = 0; % Memory of 2nd integrator
    Z_int3 = cell(1,ng-1); % Memory of 3st integrator
    Z_int4 = 0; % Memory of 4rd integrator
    for i = 1:ng-1
        Z_int3{i} = zeros(1,i);
    end

    % ----- Cancel lip radiation contribution --------------------------------

    % Define lip radiation filter
    al = [1 -d];

    Lpf = nv+1;                         % Pre-frame length

    for n = 1:nf
        % Frame-based Processing

    % For the successive removals of the estimated LPC envelopes, a
    % mean-normalized pre-frame ramp is added at the beginning of the frame
    % in order to diminish ripple. The ramp is removed after each filtering.

        % Signal frame for LPC computations
        s_gvl = s(tmarks(1,n):tmarks(2,n));
        % Signal frame for inverse filtering computations
        if n == 1
            % ----- Addition of pre-frame --------------------------------------------
            x_gvl = [linspace(-s(tmc(1)),s(tmc(1)),Lpf)' ;s(tmc(n):(tmc(n+1)-1))];    % Prepend
            idx_pf = (Lpf+1):length(x_gvl);     % Indexes that exclude the pre-frame
        else
            x_gvl =s(tmc(n):(tmc(n+1)-1));
            idx_pf = 1:length(x_gvl);
        end

        % Integration of signal using filter 1/[1 -d z^(-1)]
        % - Input signal (for LPC estimation)
        [s_gv,Z_int1] = filter(1,al,s_gvl,Z_int1);
        % - Pre-framed input signal (for LPC envelope removal)
        [x_gv,Z_int2] = filter(1,al,x_gvl,Z_int2);

        % ----- Gross glottis estimation -----------------------------------------

        % Iterative estimation of glottis with ng first order filters
        %disp(['Frame n = ',num2str(n)]);
        ag1 = lpc(s_gv.*hann(length(s_gv)),1);         % First 1st order LPC estimation

        for i = 1:ng-1
            % Cancel current estimate of glottis contribution from speech signal
            [x_v1x,Z_int3{i}] = filter(ag1,1,x_gv,Z_int3{i});	% Inverse filtering
            s_v1x = x_v1x(idx_pf);      % Remove pre-ramp

            % Next 1st order LPC estimation
            ag1x = lpc(s_v1x.*hann(length(s_v1x)),1);   % 1st order LPC

            % Update gross estimate of glottis contribution
            ag1 = conv(ag1,ag1x);	% Combine 1st order estimation with previous
        end

        % ----- Gross vocal tract estimation -------------------------------------

        % Cancel gross estimate of glottis contribution from speech signal
        [x_v1,Z_ag1] = filter(ag1,1,x_gv,Z_ag1);      % Inverse filtering
        % JC: Lip radiation compensation (integration)
        %x_v1 = filter(1,al,x_v1);
        s_v1 = x_v1(idx_pf);            % Remove pre-ramp

        % Gross estimate of the vocal tract filter
        %av1 = lpc(s_v1.*hann(length(s_v1)),nv);        % nv order LPC estimation
        [avt1,bvt1] = vtEstimation(s_v1.*hann(length(s_v1)),nv);

        % ----- Fine glottis estimation ------------------------------------------

        % Cancel gross estimate of vocal tract contribution from speech signal
        [x_g1,Z_av1] = filter(avt1,bvt1,x_gv,Z_av1);      % Inverse filtering
        % JC: Lip radiation compensation (integration)
        %x_g1 = filter(1,al,x_g1);
        s_g1 = x_g1(idx_pf);            % Remove pre-ramp

        % Fine estimate of the glottis filter
        ag = lpc(s_g1.*hann(length(s_g1)),ng);         % ng order estimation

        % ----- Fine vocal tract estimation --------------------------------------

        % Cancel fine estimate of glottis contribution from speech signal
       [x_v,Z_ag] = filter(ag,1,x_gv,Z_ag);        % Inverse filtering
        % JC: Lip radiation compensation (integration)
        %x_v = filter(1,al,x_v);
        s_v = x_v(idx_pf);              % Remove pre-ramp

        % Fine estimate of the vocal tract filter
        %av = lpc(s_v.*hann(length(s_v)),nv);          % nv order LPC estimation
         [av,bv] = vtEstimation(s_v.*hann(length(s_v)),nv);

        % Calcular el pols glotal fent el filtre invers del tracte vocal i
        % cancel·lant l'efecte de radiació dels llavis aplicant filtre d'integració
        % amb "al"
        [g,Z_av] = filter(av,bv,x_gv,Z_av);
        [gd,Z_int4] = filter(al,1,g,Z_int4);
        gd = gd(idx_pf);
        g = g(idx_pf);

        Av{n} = av;
        Ag{n} = ag;

        G(tmc(n):(tmc(n+1)-1)) = g;
        Gd(tmc(n):(tmc(n+1)-1)) = gd;

    end
end