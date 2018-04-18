 %% MT
 
out  = se2_pubFigs(Dall , 'MT','RandvsStructCommpare'); 
out  = se2_pubFigs(Dall , 'MT','RandStructAcrossDays' , 'poolDays' , 0); 
out  = se2_pubFigs(Dall , 'MT','compareLearning' , 'poolDays' , 0);
out  = se2_pubFigs(Dall , 'MT','LearningEffectShade' , 'poolDays' , 0);



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

out  = se2_pubFigs(Dall , 'IPI_asymptote','Actual&fitHorz', 'poolDays' , 0, 'MaxIter' , 150);
out  = se2_pubFigs(Dall , 'IPI_asymptote','plotCoef', 'poolDays' , 0, 'MaxIter' , 200);


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

%% SeqType sig test
clear hpval
horz = {[1] [2] [3] [4] [5] [6:13]};
for h = [1:8,13]
    for d  = 1:5
        stats = se2_SigTest(Dall , 'MT' , 'seqNumb' , [0:2] , 'Day' , d , 'Horizon' , [h],...
            'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
            'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , [1:13]);
        hpval(d,h) = stats.eff(2).p;
        close all
    end
end
pval = hpval;
pval(pval>0.05) = NaN;
%% single subject sig test
clear pval effHorz
seqN = {[0] , [1 2]};
allcount = 1;
for sn = 1:13
    dcount = 1;
    for d  = [1,5]
        for sq = 1:length(seqN)
            effHorz.Day(allcount,1) = d;
            effHorz.SN(allcount,1) = sn;
            effHorz.sq(allcount,1) = sq;
            for h = [1:8]
                stats = se2_SigTest(Dall , 'MT' , 'seqNumb' , seqN{sq} , 'Day' , d , 'Horizon' , [h:13],...
                    'PoolDays' , 0,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
                    'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , sn);
                pval{sq}(sn,h,dcount) = stats(1);
                close all
            end
            temp = squeeze(pval{sq}(sn,:,dcount));
            effHorz.effH(allcount,1) = find(temp>0.05 ,1 , 'first');
            allcount = allcount+1;
        end
        dcount = dcount+1;
    end
end

figure('color' , 'white')
barplot([effHorz.Day effHorz.sq] ,effHorz.effH);
set(gca , 'FontSize' , 18 , 'YLim' , [1 4] , 'YTick' , [1 2 3 4])
title('The Effective Window Size Grows Significantly from First to Last Day in Random Sequences' , 'FontSize' , 24)
ylabel('Viewing Window Size', 'FontSize' , 21)


        
        
        
        
        
        
        
        
        
        
        
