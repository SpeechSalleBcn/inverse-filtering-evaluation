% Description
%
% Function to compute Spectral Distortion between glottal flow from GIF technique and its
% gound truth reference
%
% Inputs
%  gf_GT    : [samples] [Nx1] Glottal flow ground truth
%  gf       : [samples] [Nx1] Glottal flow estimation
%  fs       : [Hz]      [1x1] sampling frequency
%  GCI      : [s]       [Mx1] Glottal closure instants 
%
% Outputs
%  SD      : [s,samples] [Mx2] Spectral distortion
%  SD_comp
%
% References
%       [1] Thomas Drugman, Baris Bozkurt, and Thierry Dutoit,
%       "A comparative study of glottal source estimation techniques,"
%       Computer Speech & Language, vol. 26, no. 1,
%       pp. 20–34, 2012.
%
%       [2] Nordin, F., Eriksson,
%       "A speech spectrum distortion measure with interframe memory,"
%       in: Proc. ICASSP, pp. 717–720, 2001.
%
% Copyright (c) 2022 La Salle Campus Barcelona
% LICENSE
% TODO: add license
%
% Author 
%  Joan Claudi Socoró joanclaudi.socoro@salle.url.edu

function [SD,SD_comp] = get_sd(gf_GT,gf,fs,GCI)

% according to [1] the spectral distortion is computed between 20 and 4000 Hz
NFFT=1024;
freq = fs*(0:(NFFT/2))/NFFT;
fmin=20;
fmax=4000;
selfreqs=false(1,NFFT);
selfreqs(freq>=fmin & freq<=fmax)=true;

GCI = round(GCI*fs)+1;

%% Initial settings
F0min=20;
F0max=500;
SD=zeros(1,length(GCI));
SD_comp=zeros(1,length(GCI));

%% Do processing
for n=1:length(GCI)
    
    % Get glottal pulse compensated for zero-line drift
    if n==1
        start=1;
        stop=GCI(n);
        T0=GCI(n+1)-GCI(n);
    else start=GCI(n-1);
        stop=GCI(n);
        T0=GCI(n)-GCI(n-1);
    end
    F0=fs/T0; 
    
    if isinf(F0)==0 && T0~=0 && F0 > F0min && F0<F0max 
       
        % Compensate for zero-line drift in the glottal flow pulse
        gf_comb=[gf(start) gf(stop)];
        gf_GT_comb=[gf_GT(start) gf_GT(stop)];
        if start~=stop
            if length(gf_comb)>1
                line=interp1(linspace(1,stop-start+1,2),gf_comb, ...
                    1:stop-start+1);
                line_GT=interp1(linspace(1,stop-start+1,2),gf_GT_comb, ...
                    1:stop-start+1);
            else
                line=gf_comb;
                line_GT=gf_GT_comb;
            end
        else
            line=0;
            line_GT=0;
        end
        gf_seg=gf(start:stop);
        gf_seg_comp=gf_seg(:)-line(:);

        gf_GT_seg=gf_GT(start:stop);
        gf_GT_seg_comp=gf_GT_seg(:)-line_GT(:);
        
        S_estimated_comp = fft(gf_seg_comp,NFFT);
        S_reference_comp = fft(gf_GT_seg_comp,NFFT);
        SD_comp(n) = sqrt( (2/(2*fmax)) * sum( (20*log10(abs(S_estimated_comp(selfreqs)./S_reference_comp(selfreqs)))) .^2 ) );

        S_estimated = fft(gf_seg,NFFT);
        S_reference = fft(gf_GT_seg,NFFT);
        SD(n) = sqrt( (2/(2*fmax)) * sum( (20*log10(abs(S_estimated(selfreqs)./S_reference(selfreqs)))) .^2 ) );
        
    end
   
end

%% Add in time to parameters
SD_comp=[GCI(:)/fs SD_comp(:)];
SD=[GCI(:)/fs SD(:)];

