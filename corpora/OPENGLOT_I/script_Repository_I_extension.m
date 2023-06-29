% Script to generate OPENGLOT - Repository I vowels with female extension

% output path
outPath = ['..' filesep 'extended'];
% set desired parameter values
fs = 8000;              % sampling frequency
framelength = 0.2;      % length of signal (sec)

gender = {'Male','Female'};
f0_values = {(100:20:360),(180:20:500)};               % fundamental frequency (Hz)
vowel_values = {'a','e','i','o','u','ae'};             % 'a','e','i','o','u', or 'ae'
phonation_values = {'m','b','w','c'};
phonation_label = {'normal','breathy','whispery','creaky'};
phonation_dict = containers.Map(phonation_values,phonation_label);

for genderId=1:length(gender)
    % prepare output subfolder
    outSubfolder= [outPath filesep gender{genderId}];
    if ~exist(outSubfolder,'dir'); mkdir(outSubfolder); end
    % obtain the parameter combinations
    [vowels, phonations, f0s] = ndgrid(vowel_values, phonation_values, f0_values{genderId});
    vowels=vowels(:); phonations=phonations(:); f0s=f0s(:);
    
    % Run the generation function for each combination of parameters
    for i=1:length(vowels)
        filename = [gender{genderId} '_' upper(vowels{i}) '_' phonation_dict(phonations{i}) '_' num2str(f0s(i)) 'Hz.wav'];
        [S, dU, U, a] = synthFrame(vowels{i},phonations{i},framelength,f0s(i),fs,[],gender{genderId});
        % normalize before writing them into a wav
        Snorm = S(:)/max(abs(S));
        Unorm = U(:)/max(abs(U));
        % the normalization factors are written into the wav as Comment
        comment=['S_fact:' num2str(max(abs(S))) ',' 'U_fact:' num2str(max(abs(U)))];
        audiowrite([outSubfolder filesep filename],[Snorm Unorm],fs,'Comment',comment);
    end
end


