function [ U, u ] = lf( EE, Ra, Rg, Rk, f0, fs )
%LF Synthesizes one LF glottal flow and glottal flow derivative pulse using given parameters
%   Detailed explanation goes here

%EE=1;Ra=0.01;Rg=1.17;Rk=0.34;f0=100; % Modal
%EE=10^(0.7/20);Ra=0.025;Rg=0.88;Rk=0.41;f0=100; % Breathy



fsH = 48000;

wg = 2*pi*Rg*f0;
Ta = Ra/f0;
Tp = 0.5/(Rg*f0);
Tn = Rk*Tp;
Te = (1+Rk)/(2*Rg*f0);
Tb = (1-(Rk+1)/(2*Rg))/f0;


% Estimate epsilon
epsilon0 = 1;
delta = 1;

while delta > 0.01
    epsilon = 1/Ta*(1-exp(-epsilon0*Tb));
    delta = abs(epsilon-epsilon0)/epsilon0;
    epsilon0 = epsilon;
end

% Iterative estimation for E0 and alpha

A2 = -EE/(epsilon^2*Ta)*(1-exp(-epsilon*Tb)*(1+epsilon*Tb));
step = 0.1;
limit = EE/100000;
E0 = 1;
alpha = 1/Te*log(-EE/(E0*sin(wg*Te)));
A1 = E0*exp(alpha*Te)/sqrt(alpha^2+wg^2)*sin(wg*Te-atan(wg/alpha))+E0*wg/(wg^2+alpha^2);
A = A1+A2;

while abs(A) > limit
   if A > 0
       E0 = E0-step;
   else
       E0 = E0+step;
   end
   
   alpha = 1/Te*log(-EE/(E0*sin(wg*Te)));
   
   Aiter = A2 + E0*exp(alpha*Te)/sqrt(alpha^2+wg^2)*sin(wg*Te-atan(wg/alpha))+E0*wg/(wg^2+alpha^2);
   
   if sign(A) ~= sign(Aiter)
       step = step/2;
   end
    A = Aiter;
   
end

% Synthesize the LF-pulse using the given parameters

N1 = floor(fsH*Te);
N2 = floor(fsH*Tb+1);
if N1+N2 < floor(fs/f0)
   N1 = floor(fs/f0)-(N1+N2); 
end

t1 = linspace(0,Te,N1);
t2 = linspace(Te,Te+Tb,N2+1);
t2 = t2(2:end);

u1 = E0*exp(alpha*t1).*sin(wg*t1);
u2 = -EE/(epsilon*Ta).*(exp(-epsilon*(t2-Te))-exp(-epsilon*Tb));

U1 = E0*exp(alpha*t1).*sin(wg*t1-atan(wg/alpha))/sqrt(alpha^2+wg^2) + E0*wg/(alpha^2+wg^2);
U2 = EE/(epsilon^2*Ta)*(exp(-epsilon*(t2-Te))+epsilon*exp(-epsilon*Tb)*(t2-(1/f0+1/epsilon)));

u = [u1 u2];
U = [U1 U2];

% Lowpass-filter and resample the synthesized signals
u = resample(u,fs,fsH);
U = resample(U,fs,fsH);

end

