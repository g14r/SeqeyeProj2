 %% MT
 
out  = se2_pubFigs(Dall , 'MT','RandvsStructCommpare'); 
out  = se2_pubFigs(Dall , 'MT','RandStructAcrossDays' , 'poolDays' , 1); 
out  = se2_pubFigs(Dall , 'MT','LearningEffect' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'MT','compareLearning' , 'poolDays' , 1);

out  = se2_pubFigs(Dall , 'RT','RandvsStructCommpare','poolDays' , 1); 
out  = se2_pubFigs(Dall , 'RT','RandStructAcrossDays' , 'poolDays' , 0); 
out  = se2_pubFigs(Dall , 'RT','BoxAcrossDays' , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'RT','BoxAcrosSeqType' , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'RT','LearningEffectHeat' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'RT','LearningEffectShade' , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'RT','compareLearning' , 'poolDays' , 1);


out  = se2_pubFigs(Dall , 'MT_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeDayz', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeSeqType', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 50);


out  = se2_pubFigs(Dall , 'IPI','IPIFullDispHeat', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','IPIFullDispShade', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','compareLearning', 'poolDays' , 0);

out  = se2_pubFigs(Dall , 'test_MT_asymptote','', 'poolDays' , 1);


%% MT seg test
clear hpval
S = {[0] , [1:2]};
% dayz = {[1] [2 3] [4 5]};
dayz = {[1] [2] [ 3] [4] [ 5]};

for s = 1:length(S)
    for d = 1:length(dayz)
        for h = 1:5
            stats = se2_SigTest(Dall , 'MT'  , 'seqNumb' , S{s} , 'Day' , dayz{d} , 'Horizon' , [h:13],...
                'PoolDays' , 1,'whatIPI','WithinBetweenChunk','PoolSequences' , 0 ,...
                'PoolHorizons' , [6:13]);
            close all
            hpval{s}(d,h) = stats.eff(2).p;
        end
    end
    hpval{s}( hpval{s}<=0.05) = 1;
    hpval{s}( hpval{s}~=1) = 0;
end
%% IPI seg test
se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0:2] , 'Day' , [2 3] , 'Horizon' , [4,6:13],...
     'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
     'PoolHorizons' , [6:13]);
 %% b = []
 b(1) + (b(2) - b(1))*exp(-(x-1)/b(3))