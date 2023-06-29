function [G,Gd,H_VTA,H_G]=qcp_olaf(s,fs,VTorder,GSorder,lipRad,DQ,PQ,RQ,causality,STcompens,tmarks,gcis)
%[sg,Hvt2,e_ar,Hg,dg]=qcp_olaf(s,times,VTorder,GSorder,lipRad,DP,PQ,RQ,causality,STcompensation)
%   
% Quasi Closed-Phase glottal inverse filtering.
% New version that incorporates OLA processing (variable tmarks) (Joan Claudi Socoró)
%
% Inputs:
%
% s:        [Nx1]	Speech signal 
% fs:        [1x1] Sampling frequency (in Hz)
% VTorder:  [1x1]	Order of LP analysis for vocal tract (def. 48)
% GSorder:  [1x1]	Order of LP analysis for glottal source (def. 3)
% lipRad:   [1x1]	Leaky integration coefficient (def. 0.99)
% DQ:       [1x1]   Duration Quotient (between [0.4,1])
% PQ:       [1x1]   Position Quotient (between [0,0.1 or 0.2])
% RQ:       [1x1]   Ramp Quotient (between 0.14??)
% causality [1x1]   Flag that determines the type of FIR filters
%                   ('causal','noncausal')
% STcompens [1x1]   Flag that enables Spectral Tilt compensation (ref. ...)
% (tmarks): [Nf x 2] matrix of Nf x 2 where (tmarks(n,1),tmarks(n,2)) are the
%       interval of samples of frame n for the OLA processing
%
% Outputs:
%
% G:  [Nx1] Glottal flow signal
% Gd:  [Nx1] Glottal flow derivative signal
% H_VTA:  [nf x VTorder] H_VTA(n,1:VTorder) are the LP denominator coefficients of vocal tract contribution of frame n
% H_G:  [nf x GSorder]	H_G(n,1:ng) are the LP coefficients of glottis contribution of frame n
%
% Modifed version by
% Joan Claudi Socoró joanclaudi.socoro@salle.url.edu
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

% JC. Comented
% retsig=1;
% if ~isa(s,'signal')
%   fs = argin{1} ; argin = argin(2:end);
%   s=signal(s,fs);
%   retsig=0;
% end
% fs = s.fs;

% p = get_if_order(s.fs);
% g = 4;

% AR default options

% aropt = struct;
% options = struct;

% integrator coeff
% rho = 0.99;
% causality = 'noncausal';

%DQ and PQ
%DQ = 0.7;
%PQ = 0.1;

% window function
% winfunc = @hamming;
%winfunc = @hanning;

% remove positive real poles
remove_real_poles = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the options

rho = lipRad;
p = VTorder;automodel = false;
g = GSorder;
%T0 = (pmarks(2,n) - pmarks(1,n) + 1)/2;

%f0 = fs/T0;
winfunc = 'hanning'; % Idem que per IAIF i GFM-IAIF

% if ~isempty(argin)
%   options = argin{1};
%   if isfield(options,'p')
%     p = options.p;
%     automodel = false;
%   end
%   if isfield(options,'automodel')
%       automodel = options.automodel;
%   end
%   if isfield(options,'g')
%     g = options.g;
%   end
%   if isfield(options,'rho')
%     rho = options.rho;
%   end
%   if isfield(options,'dq')
%     DQ = options.dq;
%   end
%   if isfield(options,'pq')
%     PQ = options.pq;
%   end
%   if isfield(options,'nramp')
%     Nramp = options.nramp;
%   end
%   if isfield(options,'winfunc')
%     winfunc = options.winfunc;
%   end
%   if isfield(options,'aropt')
%     aropt = mergestruct(aropt,options.aropt);
%   end
%   if isfield(options,'diffout')
%     diffout = options.diffout;
%   end
%   if isfield(options,'f0')
%     aropt.f0 = options.f0;
%   else
%     aropt.f0 = find_f0(s);
%   end
%   if isfield(options,'flist') 
%     aropt.flist = options.flist;
%   end
%   if isfield(options,'remove_real_poles') 
%     remove_real_poles = options.remove_real_poles;
%   end
%   if isfield(options,'causality')
%       if options.causality == "causal"
%           causality = 'causal';
%       else
%           causality = 'noncausal';
%       end
%   end
%          
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% processing begins here
%x = s;
%f0 = aropt.f0;
% note that the highpass filtering (block 1) needs to be done in advance!
%[gci_ins goi_ins] = gci(x,aropt.f0,fs); % GCI estimation

% [gci_ins, ~, ~] = gci_sedreams(x, fs, aropt.f0,1);

% figure
% plot(x)
% hold on
% stem(gci_ins,ones(1,length(gci_ins))*-.1,'r');
% hold off

wmin = 0.00001; % Minimum value of weighting function
%DQ = 0.7;  Duration Quotient (relative to fundamental period) 0.4 -- 1
%PQ = 0.1; % Position Quotient (relative to fundamental period) 0.0, 0.05, 0.1
%Nramp = round(fs/8000*7); % Length of linear ramp (in samples)

%W = makeW(x,p,DQ,PQ,wmin,Nramp,gci_ins,fs); % Use this to manually select PQ and DQ values
% w = makeW_adapt(x,p,gci_ins,f0); % Use this for predefined, pitch-adaptive values for weighting function


if iscolumn(s)
    s = transpose(s);
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

% Initialization of FIR filters memories (inverse filters)
Z_vt = zeros(1,p); % Memory of inverse WLP vocal tract filter

Z_pre = 0; % Memory of pre-emphasis
Z_int = 0; % Memory of 1st integrator
Z_int2 = 0; % Memory of 1st integrator
Z_st = 0; % Memory of spectral tilt
Z_der = 0;  % Memory of derivative for Glottal flow derivative computation

% Initialization of output matrices
H_G = cell(1,nf);
H_VTA = cell(1,nf);

Lpf = p + 1;
for n = 1:nf
    % Frame-based Processing
    % Signal frame for LPC computations
    x = s(tmarks(1,n):tmarks(2,n));

    % Serach for AMD function of the current frame
    pos_gci_frame = find((gcis>=tmarks(1,n))&(gcis<=tmarks(2,n)));
    gci_ins = gcis(pos_gci_frame) - tmarks(1,n) + 1;
    T0 = mean(gci_ins(2:end)-gci_ins(1:end-1));
    Nramp = round(RQ*T0);
    w = makeW(x,p,DQ,PQ,wmin,Nramp,gci_ins,fs); % Use this to manually select PQ and DQ values
    
    % Signal frame for inverse filtering computations
    if n == 1
        % ----- Addition of pre-frame --------------------------------------------
        signal = [linspace(-s(tmc(n)),s(tmc(n)),Lpf)  s(tmc(n):(tmc(n+1)-1))];    % Prepend
        idx = (Lpf+1):length(signal);     % Indexes that exclude the pre-frame
    else
        signal =s(tmc(n):(tmc(n+1)-1));
        idx= 1:length(signal);
    end

    [s2,Z_pre] = filter([1 -1],1,x,Z_pre); % Pre-emphasis
    sw = win(s2,winfunc);
    [Hvt, e_ar] = wlp(sw,w,p);
    Hvt = Hvt(:)';

    if remove_real_poles
        [z,p,k] = tf2zp(1,Hvt);
        p_ = p(p<0 | abs(imag(p)) > 1e-15);
        [dummy,Hvt] = zp2tf(z,p_,k);
    end

    Hvt2 = Hvt;

    [dg,Z_vt] = filter(Hvt2,1,signal,Z_vt);
    dg = dg(idx);
    [sg,Z_int] = integrate(dg,rho,causality,Z_int);
    Hg = lpc(win(sg,winfunc),g);   % JC. Abans posava Hg = lpc(win(sg,winfunc),g) perÃ² aleshores Hg era una matriu, pel fet que win(sg',winfunc) era un vector columna (??)

    if STcompens
        % Spectral Tilt Compensation procedure for the QCP
        % method

        % get VT frequency response from prediction filter coefficients
        [Hgs,ws]=freqz(1,Hg);
        HVT=freqz(1,Hvt);

        % estimate tilt with a first order LP filter
        [bvtt,avtt]=invfreqz(HVT,ws,0,1);

        % spectral tilt to freq response
        nfft=512;
        HVTT=freqz(bvtt,avtt,nfft);
        
        % Frame filtering and assignation to the filtered excitation using
        % inverse low-order filter of Hvt trend
        %[g_frame,Z] = filter(bvtt,avtt,g_frame,Z);
        [sg,Z_st] = filter(bvtt,avtt,sg,Z_st);

        % Derivative for Glottal flow derivative computation (JC)
        [dg, Z_der] = filter([1 -rho],1,sg,Z_der);

        % JC: Recalcular Hg
        Hg = lpc(win(sg,winfunc),g);

        % JC: Recalcular Hvt ??
    end

    H_VTA{n} = Hvt;
    H_G{n} = Hg;
    G(tmc(n):(tmc(n+1)-1)) = sg;
    Gd(tmc(n):(tmc(n+1)-1)) = dg;    

% Only for debugging
%     figure(1);clf
%     subplot(311)
%     plot(1:length(x),x,1:length(x),w(1:length(x))*max(abs(x))/max(abs(w(1:length(x)))),gci_ins,x(gci_ins),'ro');
%     axis tight;
%     legend('Input signal frame','Weighting function')
%     subplot(312)
%     plot(1:length(sg),sg,1:length(sg),dg);axis tight;
%     legend('Glottal flow','Glottal flow derivative')
%     subplot(313);
%     plot(G(1:(tmc(n+1)-1)));axis tight;
%     pause

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y]=win(x,winfunc)
y=feval(winfunc,length(x))'.*x;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,Z]=integrate(x,rho,causality,Z)
if strcmp(causality, 'causal')
   [y,Z] = filter(1,[1 -rho],x,Z);
else
    [y,Z] = filter(1,[1 -rho],flip(x),Z);
    y = -flip(y);
end
end
%y=filter(1,[1 -rho],x,'f0',f0,causality);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = get_if_order(fs)

% GET_IF_ORDER  Get suitable inverse filter order for the sampling frequency
%
%  P = GET_IF_ORDER(FS)
%  Returns a suitable order for DAP/LPC inverse filter for a given
%  FS.

p = round(fs/1000)+2;
end

function [ w ] = makeW( x, p, DQ, PQ, d, Nramp, gci_ins, fs )
%   Create a AME weight function for frame x for LPC order p
%   DQ = duration quotient (from 0 to 1)
%   PQ = Position Quotient (from 0 to 1)
%   d = minimum value of the weight function
%   Nramp = length of the linear ramp (in samples)
%   gci_ins = Glottal Closure Instants of the frame x

N = length(x);
if Nramp > 0
UPramp = linspace(d,1,2+Nramp);
UPramp = UPramp(2:end-1);
DOWNramp = UPramp(end:-1:1);
end

if DQ+PQ > 1
    DQ = 1-PQ;
end

w = d.*ones(1,N+p);

for i = 1:length(gci_ins)-1
   T = gci_ins(i+1)-gci_ins(i);
   T1 = round(DQ*T);
   T2 = round(PQ*T);
   while T1+T2 > T
       T1 = T1-1;
   end
   w(gci_ins(i)+T2:gci_ins(i)+T2+T1-1) = 1;
   if Nramp > 0
       w(gci_ins(i)+T2:gci_ins(i)+T2+Nramp-1) = UPramp;
       if gci_ins(i)+T2+T1-Nramp > 0
           w(gci_ins(i)+T2+T1-Nramp:gci_ins(i)+T2+T1-1) = DOWNramp;
       end
   end
end

if length(gci_ins)==1
    T = gci_ins(1);
    T1 = round(DQ*T);
    T2 = round(PQ*T);
end

if isempty(gci_ins)
    Nend = 0;
    T2 = 0;
    T1 = 1;
else
    Nend = N-(T2+gci_ins(i+1));
end

if T2+gci_ins(i+1) < N
    if T1+T2 < Nend
        w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+T1-1) = 1;
        if Nramp > 0
            w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+Nramp-1) = UPramp;
            w(gci_ins(i+1)+T2+T1-Nramp:gci_ins(i+1)+T2+T1-1) = DOWNramp;
        end
    else
        T1 = Nend-T2;
                w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+T1-1) = 1;
        if Nramp > 0
            w(gci_ins(i+1)+T2:gci_ins(i+1)+T2+Nramp-1) = UPramp;
        end
    end
end


end

