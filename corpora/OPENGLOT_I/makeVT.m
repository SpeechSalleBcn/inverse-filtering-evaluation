function [ b, a ] = makeVT( f1, f2, f3, fs )
%MAKEVT returns an all-pole transfer function with specified 3 formant
%frequencies
%   Detailed explanation goes here
% Details Gold et. al: "Analysis of Digital and Analog Formant Synthesizers", 1968: 
% [a] :
% f1=730; f2=1090; f3=2440; fs=8000;

f4 = 3500; g4 = 175/2;

g1 = 60/2;g2 = 100/2; g3 = 120/2; % Half-bandwidth frequencies

r1 = exp(-2*pi*g1/fs);
r2 = exp(-2*pi*g2/fs);
r3 = exp(-2*pi*g3/fs);
r4 = exp(-2*pi*g4/fs);

b1 = 2*pi*f1;
b2 = 2*pi*f2;
b3 = 2*pi*f3;
b4 = 2*pi*f4;

K1 = 1+r1^2-2*r1*cos(b1/fs);
K2 = 1+r2^2-2*r2*cos(b2/fs);
K3 = 1+r3^2-2*r3*cos(b3/fs);
K4 = 1+r4^2-2*r4*cos(b4/fs);

[Z1, P1, K1] = tf2zpk(K1,[1, -2*r1*cos(b1/fs), r1^2]);
[Z2, P2, K2] = tf2zpk(K2,[1, -2*r2*cos(b2/fs), r2^2]);
[Z3, P3, K3] = tf2zpk(K3,[1, -2*r3*cos(b3/fs), r3^2]);
[Z4, P4, K4] = tf2zpk(K4,[1, -2*r4*cos(b4/fs), r4^2]);

[b, a] = zp2tf([Z1; Z2; Z3; Z4], [P1; P2; P3; P4], K1*K2*K3*K4);




end

