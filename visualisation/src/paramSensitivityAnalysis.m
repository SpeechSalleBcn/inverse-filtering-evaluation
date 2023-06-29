path = 'D:\One Drive -Salle\La Salle\DOA-DdR-GTM-FEMVoQ - Documents\Corpus\OPENGLOT\specialIssue\OpenGlotIExtended\optimizationResults\mat';
filename = 'Female_A_breathy_80Hz_QCP.mat';
load([path filesep filename]);

VTs = unique(G.VTs);
LRs = unique(G.LRs);
PQs = unique(G.PQs);
DQs = unique(G.DQs);
RQs = unique(G.RQs);

sensResults = 0;
hnd = waitbar(0,"Computing sensitivity results ...");
N = length(VTs)*length(LRs)*length(PQs)*length(DQs)*length(RQs);
n = 0;
for n1 = 1:length(VTs) 
    for n2 = 1:length(LRs) 
        for n3 = 1:length(PQs) 
            for n4 = 1:length(DQs) 
                for n5 = 1:length(RQs) 
                    n = n+1;
                    waitbar(n/N,hnd);
                    pos1 = (G.VTs == VTs(n1));
                    pos2 = (G.LRs == LRs(n2));
                    pos3 = (G.PQs == PQs(n3));
                    pos4 = (G.DQs == DQs(n4));
                    pos5 = (G.RQs == RQs(n5));
                    pos = pos1.*pos2.*pos3.*pos4.*pos5;
                    sensResults = sensResults + var(fitresult(pos ==1));
                end
            end
        end
    end
end
close(hnd);
sensResults = sensResults/N;

