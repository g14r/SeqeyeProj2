%% plotting recipes

out  = se2_pubFigs(Dall , 'Error','ErrorvsWindow' , 'subjnum' , [1:15]);


out  = se2_pubFigs(Dall , 'MT','RandvsStructCommpare' , 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'MT','RandStructAcrossDays' , 'poolDays' , 0 , 'subjnum' , [1,2,4:15]);
out  = se2_pubFigs(Dall , 'MT','compareLearning' , 'poolDays' , 0 , 'subjnum' , [1,2,4:15]);
out  = se2_pubFigs(Dall , 'MT','LearningEffectShade' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'MT','BoxFirstLastDays' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'MT','subjEffectiveHorizon' , 'poolDays' , 1, 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'MT','subjEffectiveHorizonThresh' , 'poolDays' , 1, 'subjnum' , [1:15]);


Dall = Dall2;
out  = se2_pubFigs(Dall , 'IPI','IPIFullDispsplitDay', 'poolDays' , 1 );
out  = se2_pubFigs(Dall , 'IPI','IPIFullDispsplitseqNumb', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','IPIFullDispsplitHorizon', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','compareLearning', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','compareLearning_histogram', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','IPILearningPlacement', 'poolDays' , 0, 'dayz' , {[1] [4 5]});
out  = se2_pubFigs(Dall , 'IPI','percentTotalLearning_IPIplacement', 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'IPI','subjEffectiveHorizon' , 'poolDays' , 1, 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'IPI','IPIbyTransition' , 'poolDays' , 1, 'subjnum' , [1:15],'dayz' , {[4 5]});
out  = se2_pubFigs(Dall , 'IPI','RTvsInitialIPIs' , 'poolDays' , 1, 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'IPI','initialEyeInitialIPI' , 'poolDays' , 1, 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'IPI','subjEffectiveHorizonThresh' , 'poolDays' , 1, 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'IPI','plotAverage' , 'poolDays' , 1, 'subjnum' , [1:15]);




out  = se2_pubFigs(Dall , 'RT','RandvsStructCommpare');
out  = se2_pubFigs(Dall , 'RT','RandStructAcrossDays' , 'poolDays' , 1 , 'dayz' , {[1] [2 3] [4 5]}, 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'RT','compareLearning' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'RT','LearningEffectShade' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'RT','BoxFirstLastDays' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'RT','subjEffectiveHorizon' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'RT','subjEffectiveHorizonThresh' , 'poolDays' , 1, 'subjnum' , [1:15]);




out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fitHorz', 'poolDays' , 0, 'MaxIter' , 150);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fitDayz', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeDayzTotalLearning', 'poolDays' , 0, 'MaxIter' , 50 , 'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeDay2Day', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','Actual&fit%ChangeSeqType', 'poolDays' , 0, 'MaxIter' , 50);
out  = se2_pubFigs(Dall , 'MT_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 150);
out  = se2_pubFigs(Dall , 'MT_asymptote','', 'poolDays' , 0, 'MaxIter' , 150);

out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fitHorz', 'poolDays' , 0, 'MaxIter' , 300);
out  = se2_pubFigs(Dall , 'IPI_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 200);
out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fit%ChangeDay2Day', 'poolDays' , 0, 'MaxIter' , 300);
out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fit%ChangeDayzTotalLearning', 'poolDays' , 1, 'MaxIter' , 30);



out  = se2_pubFigs(Dall , 'Eye', 'sacDurSplitDay' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'sacDurSplitseqType' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'sacAmpSplitDay' , 'isSymmetric' , 1 , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'Eye', 'sacAmpSplitseqType' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'sacFreqSplitDay' , 'isSymmetric' , 1 , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'Eye', 'sacFreqSplitseqType' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'FixDurSplitipitype' , 'isSymmetric' , 1 , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'Eye', 'FixDurSplitwindow' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'EyePrsTimePos' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'previewSplitipitype' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'previewSplitwindow' , 'isSymmetric' , 1 , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'Eye', 'previewSplitDays' , 'isSymmetric' , 1 , 'poolDays' , 1);
out  = se2_pubFigs(Dall , 'Eye', 'startlookahead' , 'isSymmetric' , 1 , 'poolDays' , 1,'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'Eye', 'presspositionlook_ahead' , 'isSymmetric' , 1 , 'poolDays' , 0,'subjnum' , [1:15]);
out  = se2_pubFigs(Dall , 'Eye', 'initialEyeInitialIPI' , 'isSymmetric' , 1 ,'poolDays' , 1, 'subjnum' , [1:15]);


se2_compareExp(Dall1 , Dall2 , 'MT')
se2_compareExp(Dall1 , Dall2 , 'RT')

for sn = 1:15
    for h  = 1:4
        stats(sn,h) = se2_SigTest(Dall , 'Ttest_MT' , 'seqNumb' , [0] , 'Day' , [4 5] , 'Horizon' , [h,6:13],...
            'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
            'PoolHorizons' , [6:13],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [sn]);
        
    end
    Hoz(sn,1) = find(stats(sn,:)>0.05 ,1 , 'first')-1;
end

%% significance test on MTs

stats = se2_SigTest(Dall , 'MT' , 'seqNumb' , [0] , 'Day' , [1:5] , 'Horizon' , [2],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1,2,4:15]);

%% significance test on RTs
stats = se2_SigTest(Dall , 'RT' , 'seqNumb' , [0] , 'Day' , [4 5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);

%% significance test on RTs
stats = se2_SigTest(Dall , 'RT' , 'seqNumb' , [103 203 303] , 'Day' , [1:5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);

%% significance test on IPIs % {[1] [2] [3] [4] [5] [6] [7] [8] [9] [10] [11] [12] [13]}
stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0] , 'Day' , [2 3] , 'Horizon' , [1],...
    'PoolDays' , 1,'whatIPI','ipistoEachother','PoolSequences' , 0 ,...
    'PoolHorizons' , [6:13],'ipiOfInterest' , num2cell([1:13]) , 'poolIPIs' , 0 , 'subjnum' , [1:15]);
% WHAT I USED IN THE PAPER
% stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0] , 'Day' , [1 : 5] , 'Horizon' , [1:13],...
%     'PoolDays' , 0,'whatIPI','ipistoEachother','PoolSequences' , 0 ,...
%     'PoolHorizons' , [0],'ipiOfInterest' , {[1 2] [5:9] [12:13] } , 'poolIPIs' , 0 , 'subjnum' , [1:15]);
%%
stats = se2_SigTest(Dall , 'PercentIPIplace' , 'seqNumb' , [0] , 'Day' , [1 5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , {[1:4] [5:9]} , 'poolIPIs' , 0 , 'subjnum' , [1:15],'isSymmetric' , 1 , 'prsnumb' , [1]);
%%
stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0:1] , 'Day' , [1:5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0:2] , 'poolIPIs' , 0, 'subjnum' , [1:13]);
%%

%% significance test on IPIs
stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0] , 'Day' , [1] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','ipiOfInterestToSS','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);
%% single subject horizon significance test
stats = se2_SigTest(Dall , 'PerSubjMTHorz' , 'seqNumb' , [0] , 'Day' , [1 3 5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [7:13],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);

%% significance test for percent change in MTs (Random Structured)
stats = se2_SigTest(Dall , 'PercentLearning_CompareSeqType' , 'seqNumb' , [0:1] , 'Day' , [1:5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);


%% significance test for percent change in IPIs (random within between)
stats = se2_SigTest(Dall , 'PercentIPItype' , 'seqNumb' , [0:2] , 'Day' , [1: 5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0:2] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);

%% significance test for percent change within sequence tpe acoss dayz
stats = se2_SigTest(Dall , 'PercentLearning_Overall_withinseqType' , 'seqNumb' , [0] , 'Day' , [ 5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);
%% significance test for percent change within sequence tpe acoss dayz
stats = se2_SigTest(Dall , 'PercentLearning_Partial_withinseqType' , 'seqNumb' , [0] , 'Day' , [1 3 5] , 'Horizon' , [1:13],...
    'PoolDays' ,0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0] , 'poolIPIs' , 0 , 'subjnum' , [1:15]);
%% Saccade Frequency
stats = se2_SigTest(Dall , 'Eye_seq_sacPerSec' , 'seqNumb' , [0:2] , 'Day' , [1:5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0 2] , 'poolIPIs' , 0 , 'subjnum' , [1:15],'isSymmetric' , 1);
%%
stats = se2_SigTest(Dall , 'Eye_seq_sacAmp' , 'seqNumb' , [0] , 'Day' , [1 5] , 'Horizon' , [13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0 2] , 'poolIPIs' , 0 , 'subjnum' , [1:15],'isSymmetric' , 1);
%%
stats = se2_SigTest(Dall , 'Eye_ipi_fixDur' , 'seqNumb' , [0:2] , 'Day' , [1:5] , 'Horizon' , [1:13],...
    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0:3] , 'poolIPIs' , 0 , 'subjnum' , [1:15],'isSymmetric' , 1);
%%
stats = se2_SigTest(Dall , 'Eye_ipi_lookahead' , 'seqNumb' , [0] , 'Day' , [1 4 5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0] , 'poolIPIs' , 0 , 'subjnum' , [1:15],'isSymmetric' , 1);

%%
stats = se2_SigTest(Dall , 'Eye_ipi_lookahead_prsnumb' , 'seqNumb' , [0] , 'Day' , [1 : 5] , 'Horizon' , [1:13],...
    'PoolDays' ,1 ,'whatIPI','WithBetRand','PoolSequences' , 0 ,... 
    'PoolHorizons' , [],'ipiOfInterest' , [0] , 'poolIPIs' , 0 , 'subjnum' , [1:9 , 11,12,14],'isSymmetric' , 1 , ...
    'prsnumb' , [1:14],'poolpress' , {[5:9]});

%%
stats = se2_SigTest(Dall , 'Eye_ipi_lookahead_prsnumb_persubj' , 'seqNumb' , [0] , 'Day' , [1 5] , 'Horizon' , [1:13],...
    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
    'PoolHorizons' , [],'ipiOfInterest' , [0] , 'poolIPIs' , 0 , 'subjnum' , [1:15],'isSymmetric' , 1 , 'prsnumb' , [1]);


