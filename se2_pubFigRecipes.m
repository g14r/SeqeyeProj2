 %% plotting recipes
 
out  = se2_pubFigs(Dall , 'MT','RandvsStructCommpare'); 
out  = se2_pubFigs(Dall , 'MT','RandStructAcrossDays' , 'poolDays' , 0); 
out  = se2_pubFigs(Dall , 'MT','compareLearning' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'MT','LearningEffectShade' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'MT','BoxFirstLastDays' , 'poolDays' , 0);


out  = se2_pubFigs(Dall , 'RT','RandvsStructCommpare','poolDays' , 1); 
out  = se2_pubFigs(Dall , 'RT','RandStructAcrossDays' , 'poolDays' , 0); 
out  = se2_pubFigs(Dall , 'RT','BoxAcrossDays' , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'RT','BoxAcrosSeqType' , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'RT','LearningEffectHeat' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'RT','LearningEffectShade' , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'RT','compareLearning' , 'poolDays' , 1);


out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fitHorz', 'poolDays' , 0, 'MaxIter' , 150);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fitDayz', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeDayzTotalLearning', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeDay2Day', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeSeqType', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 150);
out  = se2_pubFigs(Dall , 'MT_asymptote','', 'poolDays' , 0, 'MaxIter' , 150);

out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fitHorz', 'poolDays' , 0, 'MaxIter' , 300);
out  = se2_pubFigs(Dall , 'IPI_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 200);
out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fit%ChangeDay2Day', 'poolDays' , 0, 'MaxIter' , 300);
out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fit%ChangeDayzTotalLearning', 'poolDays' , 0, 'MaxIter' , 300);


out  = se2_pubFigs(Dall , 'IPI','IPIFullDispHeat', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','IPIFullDispShade', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','compareLearning', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','compareLearning_histogram', 'poolDays' , 0);


out  = se2_pubFigs(Dall , 'Eye', '' , 'isSymmetric' , 0 , 'poolDays' , 0);


%% significance test recipes



%% significance test on MTs
stats = se2_SigTest(Dall , 'MT' , 'seqNumb' , [0:2] , 'Day' , [5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [5:13],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1:13]);

%% significance test on IPIs
stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0:2] , 'Day' , [1:5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0:2] , 'poolIPIs' , 0 , 'subjnum' , [1:13]);
%% single subject horizon significance test
stats = se2_SigTest(Dall , 'PerSubjMTHorz' , 'seqNumb' , [0:2] , 'Day' , [5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0 2] , 'poolIPIs' , 0 , 'subjnum' , [1:13]);

%% significance test for percent change in MTs (Random Structured)
stats = se2_SigTest(Dall , 'PercentseqType' , 'seqNumb' , [0:2] , 'Day' , [1:5] , 'Horizon' , [7:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [6:13],'ipiOfInterest' , [0 : 2] , 'poolIPIs' , 0 , 'subjnum' , [1:13]);


%% significance test for percent change in IPIs (random within between)
stats = se2_SigTest(Dall , 'PercentIPItype' , 'seqNumb' , [0:2] , 'Day' , [5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0 2] , 'poolIPIs' , 0 , 'subjnum' , [1:13]);

%%

        
        
        
        
        
        
        
        
