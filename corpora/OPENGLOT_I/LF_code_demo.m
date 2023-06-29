% Demonstration of how to use LF generation function


% set desired parameter values
fs = 8000;              % sampling frequency
framelength = 0.2;      % length of singal (sec)
f0 = 100;               % fundamental frequency (Hz)
vowel = 'ae';            % 'a','e','i','o','u', or 'ae'
LF_params = [1 0.01 1.17 0.34];  % LF parameters [Ee Ra Rg Rk]

% For reference, LF_params for different phonation types:
%modal [1 0.01 1.17 0.34];
%breathy [1.08 0.025 0.88 0.41];
%whispery [0.59 0.07 0.94 0.32];
%creaky [1.23 0.008 1.13 0.2];


%% Run the generation function

[S, dU, U, a] = synthFrame(vowel,'custom',framelength,f0,fs,LF_params);

%% Plot the results

time = (0:length(S)-1)/fs;
figure
subplot(3,1,1)
plot(time,S); xlabel('Time (s)'); ylabel('Speech pressure signal');
subplot(3,1,2)
plot(time,U); xlabel('Time (s)'); ylabel('Glottal flow')
subplot(3,1,3)
plot(time,dU); xlabel('Time (s)'); ylabel('Derivative of glottal flow')
