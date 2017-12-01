function se2_TranProb(Dall)
load([baseDir , '/CMB_34_1.mat'])
CMB = CMB_34_1;
snum = length(unique(Dall.SN));



for subjnum = 1:snum
    Dall.isWrong = Dall.AllPress ~=Dall.AllResponse;
    
    ANA_allh = getrow(Dall , ismember(Dall.SN , subjnum));
    ANA1_allh = getrow(Dall ,ismember(Dall.seqNumb , [1:2]) & ismember(Dall.SN , subjnum));
    ANA0_allh = getrow(Dall ,ismember(Dall.seqNumb , 0) & ismember(Dall.SN , subjnum));
    
    
    
    allT2 = [length(ANA1_allh.AllPress)*(size(ANA1_allh.AllPress , 2)-1) length(ANA0_allh.AllPress)*(size(ANA0_allh.AllPress , 2)-1) ...
        sum(ismember(ANA_allh.seqNumb , [0:2]))*(size(ANA_allh.AllPress , 2)-1) + sum(ismember(ANA_allh.seqNumb , [103 203 303]))*2 + sum(ismember(ANA_allh.seqNumb , [104 204 304]))*3];
    t2_Nums_allh(subjnum).Chunked = zeros(length(CMB.comb2) , 1);
    t2_Nums_allh(subjnum).Rand    = zeros(length(CMB.comb2) , 1);
    t2_Nums_allh(subjnum).All     = zeros(length(CMB.comb2) , 1);
    
    ANA1_allh.t2_Nums = zeros(length(ANA1_allh.AllPress) , size(ANA1_allh.AllPress , 2) -1);
    ANA0_allh.t2_Nums = zeros(length(ANA0_allh.AllPress) , size(ANA0_allh.AllPress , 2) -1);
    ANA_allh.t2_Nums = zeros(length(ANA_allh.AllPress) , size(ANA_allh.AllPress , 2) -1);
    for t2 = 1:length(CMB.comb2)
        
        t2_Nums_allh(subjnum).TranNumb(t2 , 1) = t2;
        t2_Nums_allh(subjnum).Transition(t2 , 1:2) = CMB.comb2(t2,:);
        for p = 1:size(ANA1_allh.AllPress , 2) -1
            t2_Nums_allh(subjnum).Chunked(t2,1) =  t2_Nums_allh(subjnum).Chunked(t2,1) + sum(ismember(ANA1_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
            t2_Nums_allh(subjnum).Rand(t2,1) =  t2_Nums_allh(subjnum).Rand(t2,1) + sum(ismember(ANA0_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
            t2_Nums_allh(subjnum).All(t2,1) =  t2_Nums_allh(subjnum).All(t2,1) + sum(ismember(ANA_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
            
            ANA1_allh.t2_Nums(ismember(ANA1_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
            ANA0_allh.t2_Nums(ismember(ANA0_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
            ANA_allh.t2_Nums(ismember(ANA_allh.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
        end
    end
    t2_Nums_allh(subjnum).Chunked = t2_Nums_allh(subjnum).Chunked/allT2(1);
    t2_Nums_allh(subjnum).Rand    = t2_Nums_allh(subjnum).Rand/allT2(2);
    t2_Nums_allh(subjnum).All    = t2_Nums_allh(subjnum).All/allT2(3);
    
    for p1 = 1:5
        for p2 = 1:5
            i = find(ismember(CMB.comb2 , [p1 p2] , 'rows'));
            PoC2(p1 , p2) = t2_Nums_allh(subjnum).All(i);
        end
        PoC2_n(p1, :) = PoC2(p1, :)/sum(PoC2(p1, :)); % sets the sum of every row to 1
        for p2 = 1:5
            i = find(ismember(CMB.comb2 , [p1 p2] , 'rows'));
            t2_Nums_allh(subjnum).All_normalized(i,1) = PoC2_n(p1, p2);
        end
    end
    
    allT3 = [length(ANA1_allh.AllPress)*(size(ANA1_allh.AllPress , 2)-2) length(ANA0_allh.AllPress)*(size(ANA0_allh.AllPress , 2)-2) ...
        sum(ismember(ANA_allh.seqNumb , [0:2]))*(size(ANA_allh.AllPress , 2)-2) + sum(ismember(ANA_allh.seqNumb , [103 203 303])) + sum(ismember(ANA_allh.seqNumb , [104 204 304]))*2];
    t3_Nums_allh(subjnum).Chunked = zeros(length(CMB.comb3) , 1);
    t3_Nums_allh(subjnum).Rand    = zeros(length(CMB.comb3) , 1);
    t3_Nums_allh(subjnum).All    = zeros(length(CMB.comb3) , 1);
    
    ANA1_allh.t3_Nums = zeros(length(ANA1_allh.AllPress) , size(ANA1_allh.AllPress , 2) -2);
    ANA0_allh.t3_Nums = zeros(length(ANA0_allh.AllPress) , size(ANA0_allh.AllPress , 2) -2);
    ANA_allh.t3_Nums = zeros(length(ANA_allh.AllPress) , size(ANA_allh.AllPress , 2) -2);
    for t3 = 1:length(CMB.comb3)
        
        
        t3_Nums_allh(subjnum).TranNumb(t3 , 1) = t3;
        t3_Nums_allh(subjnum).Transition(t3 , :) = CMB.comb3(t3,:);
        for p = 1:size(ANA1_allh.AllPress , 2) -2
            t3_Nums_allh(subjnum).Chunked(t3,1) =  t3_Nums_allh(subjnum).Chunked(t3,1) + sum(ismember(ANA1_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
            t3_Nums_allh(subjnum).Rand(t3,1) =  t3_Nums_allh(subjnum).Rand(t3,1) + sum(ismember(ANA0_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
            t3_Nums_allh(subjnum).All(t3,1) =  t3_Nums_allh(subjnum).All(t3,1) + sum(ismember(ANA_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
            
            ANA1_allh.t3_Nums(ismember(ANA1_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
            ANA0_allh.t3_Nums(ismember(ANA0_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
            ANA_allh.t3_Nums(ismember(ANA_allh.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
        end
    end
    t3_Nums_allh(subjnum).Chunked = t3_Nums_allh(subjnum).Chunked/allT3(1);
    t3_Nums_allh(subjnum).Rand = t3_Nums_allh(subjnum).Rand/allT3(2);
    t3_Nums_allh(subjnum).All = t3_Nums_allh(subjnum).All/allT3(3);
    for p12 = 1:25
        for p3 = 1:5
            i = find(ismember(CMB.comb3 , [CMB.comb2(p12 , :) , p3] , 'rows'));
            PoC3(p12 , p3) = t3_Nums_allh(subjnum).All(i);
        end
        PoC3_n(p12, :) = PoC3(p12, :)/sum(PoC3(p12, :)); % sets the sum of every row to 1
        PoC3_n(isnan(PoC3_n)) = 0;
        for p3 = 1:5
            i = find(ismember(CMB.comb3 , [CMB.comb2(p12 , :) p3] , 'rows'));
            t3_Nums_allh(subjnum).All_normalized(i,1) = PoC3_n(p12, p3);
        end
    end
    
    allT4 = [length(ANA1_allh.AllPress)*(size(ANA1_allh.AllPress , 2)-3) length(ANA0_allh.AllPress)*(size(ANA0_allh.AllPress , 2)-3) ...
        sum(ismember(ANA_allh.seqNumb , [0:2]))*(size(ANA_allh.AllPress , 2)-3)+sum(ismember(ANA_allh.seqNumb , [104 204 304]))];
    t4_Nums_allh(subjnum).Chunked = zeros(length(CMB.comb4) , 1);
    t4_Nums_allh(subjnum).Rand    = zeros(length(CMB.comb4) , 1);
    t4_Nums_allh(subjnum).All    = zeros(length(CMB.comb4) , 1);
    
    ANA1_allh.t4_Nums = zeros(length(ANA1_allh.AllPress) , size(ANA1_allh.AllPress , 2) -3);
    ANA0_allh.t4_Nums = zeros(length(ANA0_allh.AllPress) , size(ANA0_allh.AllPress , 2) -3);
    ANA_allh.t4_Nums = zeros(length(ANA_allh.AllPress) , size(ANA_allh.AllPress , 2) -3);
    for t4 = 1:length(CMB.comb4)
        t4_Nums_allh(subjnum).TranNumb(t4 , 1) = t4;
        t4_Nums_allh(subjnum).Transition(t4 , :) = CMB.comb4(t4,:);
        for p = 1:size(ANA1_allh.AllPress , 2) -3
            t4_Nums_allh(subjnum).Chunked(t4,1) =  t4_Nums_allh(subjnum).Chunked(t4,1) + sum(ismember(ANA1_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
            t4_Nums_allh(subjnum).Rand(t4,1) =  t4_Nums_allh(subjnum).Rand(t4,1) + sum(ismember(ANA0_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
            t4_Nums_allh(subjnum).All(t4,1) =  t4_Nums_allh(subjnum).All(t4,1) + sum(ismember(ANA_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
            
            ANA1_allh.t4_Nums(ismember(ANA1_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
            ANA0_allh.t4_Nums(ismember(ANA0_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
            ANA_allh.t4_Nums(ismember(ANA_allh.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
        end
    end
    t4_Nums_allh(subjnum).Chunked = t4_Nums_allh(subjnum).Chunked/allT4(1);
    t4_Nums_allh(subjnum).Rand = t4_Nums_allh(subjnum).Rand/allT4(2);
    t4_Nums_allh(subjnum).All = t4_Nums_allh(subjnum).All/allT4(3);
    for p123 = 1:125
        for p4 = 1:5
            i = find(ismember(CMB.comb4 , [CMB.comb3(p123 , :) p4] , 'rows'));
            PoC4(p123 , p4) = t4_Nums_allh(subjnum).All(i);
        end
        
        PoC4_n(p123, :) = PoC4(p123, :)/sum(PoC4(p123, :)); % sets the sum of every row to 1
        PoC4_n(isnan(PoC4_n)) = 0;
        for p4 = 1:5
            i = find(ismember(CMB.comb4 , [CMB.comb3(p123 , :) p4] , 'rows'));
            t4_Nums_allh(subjnum).All_normalized(i,1) = PoC4_n(p123, p4);
        end
    end
    t4_Nums_allh(subjnum).All_normalized(isnan(t4_Nums_allh(subjnum).All_normalized)) = 0;
    
    
    
    
    %%
    for h  = [1:8 13]
        id = ANA_allh.seqlength <=4;
        
        A = [ANA_allh.AllPress(id , :) ANA_allh.seqlength(id)];
        A(isnan(A)) = 0;
        A = unique(A , 'rows');
        for sl = 2:4
            id = A(:,end) == sl;
            CMB.Chunks{sl} = A(id,1:4);
        end
        
        ANA1 = getrow(Dall ,ismember(Dall.seqNumb , [1:2]) & ismember(Dall.SN , subjnum) & ismember(Dall.Horizon , h));
        ANA0 = getrow(Dall ,ismember(Dall.seqNumb , 0) & ismember(Dall.SN , subjnum) & ismember(Dall.Horizon , h));
        ANA = getrow(Dall ,ismember(Dall.SN , subjnum) & ismember(Dall.Horizon , h));
        
        allT2 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-1) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-1) length(ANA.AllPress)*(size(ANA.AllPress , 2)-1)];
        t2_Nums(subjnum,h).Chunked = zeros(length(CMB.comb2) , 1);
        t2_Nums(subjnum,h).Rand    = zeros(length(CMB.comb2) , 1);
        t2_Nums(subjnum,h).All    = zeros(length(CMB.comb2) , 1);
        
        ANA1.t2_Nums = zeros(length(ANA1.AllPress) , size(ANA1.AllPress , 2) -1);
        ANA0.t2_Nums = zeros(length(ANA0.AllPress) , size(ANA0.AllPress , 2) -1);
        ANA.t2_Nums = zeros(length(ANA.AllPress) , size(ANA.AllPress , 2) -1);
        for t2 = 1:length(CMB.comb2)
            t2_Nums(subjnum,h).Chunked_IPI{t2,1} = [];
            t2_Nums(subjnum,h).Rand_IPI{t2,1} = [];
            t2_Nums(subjnum,h).All_IPI{t2,1} = [];
            
            t2_Nums(subjnum,h).TranNumb(t2 , 1) = t2;
            t2_Nums(subjnum,h).Transition(t2 , 1:2) = CMB.comb2(t2,:);
            for p = 1:size(ANA1.AllPress , 2) -1
                t2_Nums(subjnum,h).Chunked(t2,1) =  t2_Nums(subjnum,h).Chunked(t2,1) + sum(ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                t2_Nums(subjnum,h).Rand(t2,1)    =  t2_Nums(subjnum,h).Rand(t2,1) + sum(ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                t2_Nums(subjnum,h).All(t2,1)     =  t2_Nums(subjnum,h).Rand(t2,1) + sum(ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows'));
                
                ANA1.t2_Nums(ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                ANA0.t2_Nums(ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                ANA.t2_Nums(ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') , p) = t2;
                
                CorID = ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') & ~sum(ANA1.isWrong(:,p:p+1) , 2);
                t2_Nums(subjnum,h).Chunked_IPI{t2} = [t2_Nums(subjnum,h).Chunked_IPI{t2} ; [ANA1.IPI(CorID , p) t2*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                
                CorID = ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') & ~sum(ANA0.isWrong(:,p:p+1) , 2);
                t2_Nums(subjnum,h).Rand_IPI{t2} = [t2_Nums(subjnum,h).Rand_IPI{t2} ; [ANA0.IPI(CorID , p) t2*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                
                CorID = ismember(ANA.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows') & ~sum(ANA.isWrong(:,p:p+1) , 2);
                t2_Nums(subjnum,h).All_IPI{t2} = [t2_Nums(subjnum,h).All_IPI{t2} ; [ANA.IPI(CorID , p) t2*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
            end
            t2_Nums(subjnum,h).MeanChunked_IPI(t2,1) = nanmean(t2_Nums(subjnum,h).Chunked_IPI{t2}(:,1));
            t2_Nums(subjnum,h).MeanRand_IPI(t2,1) = nanmean(t2_Nums(subjnum,h).Rand_IPI{t2}(:,1));
            t2_Nums(subjnum,h).MeanAll_IPI(t2,1) = nanmean(t2_Nums(subjnum,h).All_IPI{t2}(:,1));
        end
        
        
        t2_Nums(subjnum,h).Chunked = t2_Nums_allh(subjnum).Chunked;%t2_Nums(h).Chunked/allT2(1);
        t2_Nums(subjnum,h).Rand = t2_Nums_allh(subjnum).Rand;%t2_Nums(h).Rand/allT2(2);
        t2_Nums(subjnum,h).All = t2_Nums_allh(subjnum).All;%t2_Nums(h).Rand/allT2(2);
        [~ , t2_Nums(subjnum,h).sort_ID] = sort(t2_Nums(subjnum,h).All , 'descend');
        t2_Nums(subjnum,h).All_normalized = t2_Nums_allh(subjnum).All_normalized;%t2_Nums(h).Rand/allT2(2);
        [~ , t2_Nums(subjnum,h).sort_norm_ID] = sort(t2_Nums(subjnum,h).All_normalized , 'descend');
        
        
        
        t2_Nums(subjnum,h).SortedMeanIPI_Chunked = t2_Nums(subjnum,h).MeanChunked_IPI(t2_Nums(subjnum,h).sort_ID);
        t2_Nums(subjnum,h).SortedIPI_Chunked = t2_Nums(subjnum,h).Chunked_IPI(t2_Nums(subjnum,h).sort_ID);
        t2_Nums(subjnum,h).ReadyToPlot_Chunked = cell2mat(t2_Nums(subjnum,h).SortedIPI_Chunked);
        clear xtick_r xtick_c xticklab_r xticklab_c
        counter = 1;
        counter_n = 1;
        ChunkNum2 = [];
        ChunkNum2_34 = [];
        for i = 1:length(CMB.comb2)
            idd = t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,2) == t2_Nums(subjnum,h).sort_ID(i);
            if sum(ismember(CMB.Chunks{2}(:,1:2) , CMB.comb2(i,:) , 'rows'))
                ChunkNum2 = [ChunkNum2 ; [i , find(t2_Nums(subjnum,h).sort_ID == i)]];
            end
            if sum(ismember(CMB.Chunks{3}(:,1:2) , CMB.comb2(i,:) , 'rows')) | sum(ismember(CMB.Chunks{3}(:,2:3) , CMB.comb2(i,:) , 'rows'))
                ChunkNum2_34 = [ChunkNum2_34 ; [i , find(t2_Nums(subjnum,h).sort_ID == i)]];
            end
            if sum(ismember(CMB.Chunks{4}(:,1:2) , CMB.comb2(i,:) , 'rows')) | sum(ismember(CMB.Chunks{4}(:,2:3) , CMB.comb2(i,:) , 'rows')) | sum(ismember(CMB.Chunks{4}(:,3:4) , CMB.comb2(i,:) , 'rows'))
                ChunkNum2_34 = [ChunkNum2_34 ; [i , find(t2_Nums(subjnum,h).sort_ID == i)]];
            end
            if sum(idd)
                T2(subjnum,h).xticklab_c{counter} = num2str(unique(t2_Nums(subjnum,h).ReadyToPlot_Chunked(idd,2)));
                t2_Nums(subjnum,h).ReadyToPlot_Chunked(idd,7) = t2_Nums(subjnum,h).All(t2_Nums(subjnum,h).sort_ID(i));
                T2(subjnum,h).xtick_c(counter) = i;
                counter = counter + 1;
            end
            t2_Nums(subjnum,h).ReadyToPlot_Chunked(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            
            idd_norm = t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,2) == t2_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T2(subjnum,h).xticklab_a{counter_n} = num2str(unique(t2_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,2)));
                t2_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,9) = t2_Nums(subjnum,h).All_normalized(t2_Nums(subjnum,h).sort_norm_ID(i));
                T2(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t2_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
        end
        
        t2_Nums(subjnum,h).SortedMeanIPI_Rand = t2_Nums(subjnum,h).MeanRand_IPI(t2_Nums(subjnum,h).sort_ID);
        t2_Nums(subjnum,h).SortedIPI_Rand = t2_Nums(subjnum,h).Rand_IPI(t2_Nums(subjnum,h).sort_ID);
        t2_Nums(subjnum,h).ReadyToPlot_Rand = cell2mat(t2_Nums(subjnum,h).SortedIPI_Rand);
        counter = 1;
        counter_n = 1;
        for i = 1:length(CMB.comb2)
            idd = t2_Nums(subjnum,h).ReadyToPlot_Rand(:,2) == t2_Nums(subjnum,h).sort_ID(i);
            if sum(idd)
                T2(subjnum,h).xticklab_r{counter} = num2str(unique(t2_Nums(subjnum,h).ReadyToPlot_Rand(idd,2)));
                t2_Nums(subjnum,h).ReadyToPlot_Rand(idd,7) = t2_Nums(subjnum,h).All(t2_Nums(subjnum,h).sort_ID(i));
                T2(subjnum,h).xtick_r(counter) = i;
                counter = counter + 1;
            end
            t2_Nums(subjnum,h).ReadyToPlot_Rand(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            
            idd_norm = t2_Nums(subjnum,h).ReadyToPlot_Rand(:,2) == t2_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T2(subjnum,h).xticklab_a{counter_n} = num2str(unique(t2_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,2)));
                t2_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,9) = t2_Nums(subjnum,h).All_normalized(t2_Nums(subjnum,h).sort_norm_ID(i));
                T2(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t2_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
            
        end
        
        t2_Nums(subjnum,h).SortedMeanIPI_All= t2_Nums(subjnum,h).MeanAll_IPI(t2_Nums(subjnum,h).sort_ID);
        t2_Nums(subjnum,h).SortedIPI_All = t2_Nums(subjnum,h).All_IPI(t2_Nums(subjnum,h).sort_ID);
        t2_Nums(subjnum,h).ReadyToPlot_All = cell2mat(t2_Nums(subjnum,h).SortedIPI_All);
        counter = 1;
        counter_n = 1;
        for i = 1:length(CMB.comb2)
            idd = t2_Nums(subjnum,h).ReadyToPlot_All(:,2) == t2_Nums(subjnum,h).sort_ID(i);
            if sum(idd)
                T2(subjnum,h).xticklab_a{counter} = num2str(unique(t2_Nums(subjnum,h).ReadyToPlot_All(idd,2)));
                t2_Nums(subjnum,h).ReadyToPlot_All(idd,7) = t2_Nums(subjnum,h).All(t2_Nums(subjnum,h).sort_ID(i));
                T2(subjnum,h).xtick_a(counter) = i;
                counter = counter + 1;
            end
            t2_Nums(subjnum,h).ReadyToPlot_All(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            
            idd_norm = t2_Nums(subjnum,h).ReadyToPlot_All(:,2) == t2_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T2(subjnum,h).xticklab_a{counter_n} = num2str(unique(t2_Nums(subjnum,h).ReadyToPlot_All(idd_norm,2)));
                t2_Nums(subjnum,h).ReadyToPlot_All(idd_norm,9) = t2_Nums(subjnum,h).All_normalized(t2_Nums(subjnum,h).sort_norm_ID(i));
                T2(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t2_Nums(subjnum,h).ReadyToPlot_All(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
        end
        
        
        h1 = figure;
        hold on
        subjnum
        [xcoordC_2{subjnum,h},PLOTC_2{subjnum,h},ERRORC_2{subjnum,h}]  = lineplot(t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,6) , t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,1), 'subset' , ismember(t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,4) , [4 5]));
        [xcoordR_2{subjnum,h},PLOTR_2{subjnum,h},ERRORR_2{subjnum,h}]  = lineplot(t2_Nums(subjnum,h).ReadyToPlot_Rand(:,6) , t2_Nums(subjnum,h).ReadyToPlot_Rand(:,1), 'subset' , ismember(t2_Nums(subjnum,h).ReadyToPlot_Rand(:,4) , [4 5]));
        [xcoordA_2{subjnum,h},PLOTA_2{subjnum,h},ERRORA_2{subjnum,h}]  = lineplot(t2_Nums(subjnum,h).ReadyToPlot_All(:,6) , t2_Nums(subjnum,h).ReadyToPlot_All(:,1), 'subset' , ismember(t2_Nums(subjnum,h).ReadyToPlot_All(:,4) , [4 5]));
        close(h1)
        
        temp = corrcoef(t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,1) , t2_Nums(subjnum,h).ReadyToPlot_Chunked(:,7));
        C2_chunked(subjnum,h) = temp(2);
        temp  = corrcoef(t2_Nums(subjnum,h).ReadyToPlot_Rand(:,1) , t2_Nums(subjnum,h).ReadyToPlot_Rand(:,7));
        C2_random(subjnum,h) = temp(2);
        temp  = corrcoef(t2_Nums(subjnum,h).ReadyToPlot_All(:,1) , t2_Nums(subjnum,h).ReadyToPlot_All(:,7));
        C2_all(subjnum,h) = temp(2);
        
        
        
        
        clear xtick_r xtick_c xticklab_r xticklab_c
        allT3 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-2) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-2) length(ANA.AllPress)*(size(ANA.AllPress , 2)-2)];
        t3_Nums(subjnum,h).Chunked = zeros(length(CMB.comb3) , 1);
        t3_Nums(subjnum,h).Rand    = zeros(length(CMB.comb3) , 1);
        t3_Nums(subjnum,h).All     = zeros(length(CMB.comb3) , 1);
        
        ANA1.t3_Nums = zeros(length(ANA1.AllPress) , size(ANA1.AllPress , 2) -2);
        ANA0.t3_Nums = zeros(length(ANA0.AllPress) , size(ANA0.AllPress , 2) -2);
        ANA.t3_Nums = zeros(length(ANA.AllPress) , size(ANA.AllPress , 2) -2);
        for t3 = 1:length(CMB.comb3)
            t3_Nums(subjnum,h).Chunked_IPI{t3,1} = [];
            t3_Nums(subjnum,h).Rand_IPI{t3,1} = [];
            t3_Nums(subjnum,h).All_IPI{t3,1} = [];
            
            t3_Nums(subjnum,h).TranNumb(t3 , 1) = t3;
            t3_Nums(subjnum,h).Transition(t3 , :) = CMB.comb3(t3,:);
            for p = 1:size(ANA1.AllPress , 2) -2
                t3_Nums(subjnum,h).Chunked(t3,1) =  t3_Nums(subjnum,h).Chunked(t3,1) + sum(ismember(ANA1.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                t3_Nums(subjnum,h).Rand(t3,1) =  t3_Nums(subjnum,h).Rand(t3,1) + sum(ismember(ANA0.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                t3_Nums(subjnum,h).All(t3,1) =  t3_Nums(subjnum,h).All(t3,1) + sum(ismember(ANA.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows'));
                
                ANA1.t3_Nums(ismember(ANA1.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                ANA0.t3_Nums(ismember(ANA0.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                ANA.t3_Nums(ismember(ANA.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') , p) = t3;
                
                CorID = ismember(ANA1.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') & ~sum(ANA1.isWrong(:,p:p+2) , 2);
                if ~LastIPI
                    t3_Nums(subjnum,h).Chunked_IPI{t3} = [t3_Nums(subjnum,h).Chunked_IPI{t3} ; [sum(ANA1.IPI(CorID , p:p+1),2) t3*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                else
                    t3_Nums(subjnum,h).Chunked_IPI{t3} = [t3_Nums(subjnum,h).Chunked_IPI{t3} ; [sum(ANA1.IPI(CorID , p+1),2) t3*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                end
                
                CorID = ismember(ANA0.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') & ~sum(ANA0.isWrong(:,p:p+2) , 2);
                if ~LastIPI
                    t3_Nums(subjnum,h).Rand_IPI{t3} = [t3_Nums(subjnum,h).Rand_IPI{t3} ; [sum(ANA0.IPI(CorID , p:p+1),2) t3*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                else
                    t3_Nums(subjnum,h).Rand_IPI{t3} = [t3_Nums(subjnum,h).Rand_IPI{t3} ; [sum(ANA0.IPI(CorID , p+1),2) t3*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                end
                
                CorID = ismember(ANA.AllPress(:,p:p+2) ,CMB.comb3(t3,:), 'rows') & ~sum(ANA.isWrong(:,p:p+2) , 2);
                if ~LastIPI
                    t3_Nums(subjnum,h).All_IPI{t3} = [t3_Nums(subjnum,h).All_IPI{t3} ; [sum(ANA.IPI(CorID , p:p+1),2) t3*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                else
                    t3_Nums(subjnum,h).All_IPI{t3} = [t3_Nums(subjnum,h).All_IPI{t3} ; [sum(ANA.IPI(CorID , p+1),2) t3*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                end
            end
            t3_Nums(subjnum,h).MeanChunked_IPI(t3,1) = nanmean(t3_Nums(subjnum,h).Chunked_IPI{t3}(:,1));
            t3_Nums(subjnum,h).MeanRand_IPI(t3,1) = nanmean(t3_Nums(subjnum,h).Rand_IPI{t3}(:,1));
            t3_Nums(subjnum,h).MeanAll_IPI(t3,1) = nanmean(t3_Nums(subjnum,h).All_IPI{t3}(:,1));
        end
        t3_Nums(subjnum,h).Chunked = t3_Nums_allh(subjnum).Chunked;%t3_Nums(h).Chunked/allT3(1);
        t3_Nums(subjnum,h).Rand = t3_Nums_allh(subjnum).Rand;%t3_Nums(h).Rand/allT3(2);
        t3_Nums(subjnum,h).All = t3_Nums_allh(subjnum).All;%t3_Nums(h).Rand/allT3(2);
        [~ , t3_Nums(subjnum,h).sort_ID] = sort(t3_Nums(subjnum,h).All , 'descend');
        t3_Nums(subjnum,h).All_normalized = t3_Nums_allh(subjnum).All_normalized;%t3_Nums(h).Rand/allT3(2);
        [~ , t3_Nums(subjnum,h).sort_norm_ID] = sort(t3_Nums(subjnum,h).All_normalized , 'descend');
        
        t3_Nums(subjnum,h).SortedMeanIPI_Chunked = t3_Nums(subjnum,h).MeanChunked_IPI(t3_Nums(subjnum,h).sort_ID);
        t3_Nums(subjnum,h).SortedIPI_Chunked = t3_Nums(subjnum,h).Chunked_IPI(t3_Nums(subjnum,h).sort_ID);
        t3_Nums(subjnum,h).ReadyToPlot_Chunked = cell2mat(t3_Nums(subjnum,h).SortedIPI_Chunked);
        counter = 1;
        counter_n = 1;
        ChunkNum3 = [];
        ChunkNum3_4 = [];
        for i = 1:length(CMB.comb3)
            
            idd = t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,2) == t3_Nums(subjnum,h).sort_ID(i);
            if sum(ismember(CMB.Chunks{3}(:,1:3) , CMB.comb3(i,:) , 'rows'))
                ChunkNum3 = [ChunkNum3 ; [i , find(t3_Nums(subjnum,h).sort_ID == i)]];
            end
            if sum(ismember(CMB.Chunks{4}(:,1:3) , CMB.comb3(i,:) , 'rows')) | sum(ismember(CMB.Chunks{4}(:,2:4) , CMB.comb3(i,:) , 'rows'))
                ChunkNum3_4 = [ChunkNum3_4 ; [i , find(t3_Nums(subjnum,h).sort_ID == i)]];
            end
            if sum(idd)
                T3(subjnum,h).xticklab_c{counter} = num2str(unique(t3_Nums(subjnum,h).ReadyToPlot_Chunked(idd,2)));
                t3_Nums(subjnum,h).ReadyToPlot_Chunked(idd,7) = t3_Nums(subjnum,h).All(t3_Nums(subjnum,h).sort_ID(i));
                T3(subjnum,h).xtick_c(counter) = i;
                counter = counter + 1;
            end
            t3_Nums(subjnum,h).ReadyToPlot_Chunked(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            
            idd_norm = t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,2) == t3_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T3(subjnum,h).xticklab_a{counter_n} = num2str(unique(t3_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,2)));
                t3_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,9) = t3_Nums(subjnum,h).All_normalized(t3_Nums(subjnum,h).sort_norm_ID(i));
                T3(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t3_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
            
        end
        
        
        t3_Nums(subjnum,h).SortedMeanIPI_Rand = t3_Nums(subjnum,h).MeanRand_IPI(t3_Nums(subjnum,h).sort_ID);
        t3_Nums(subjnum,h).SortedIPI_Rand = t3_Nums(subjnum,h).Rand_IPI(t3_Nums(subjnum,h).sort_ID);
        t3_Nums(subjnum,h).ReadyToPlot_Rand = cell2mat(t3_Nums(subjnum,h).SortedIPI_Rand);
        counter = 1;
        counter_n = 1;
        for i = 1:length(CMB.comb3)
            idd = t3_Nums(subjnum,h).ReadyToPlot_Rand(:,2) == t3_Nums(subjnum,h).sort_ID(i);
            if sum(idd)
                T3(subjnum,h).xticklab_r{counter} = num2str(unique(t3_Nums(subjnum,h).ReadyToPlot_Rand(idd,2)));
                t3_Nums(subjnum,h).ReadyToPlot_Rand(idd,7) = t3_Nums(subjnum,h).All(t3_Nums(subjnum,h).sort_ID(i));
                T3(subjnum,h).xtick_r(counter) = i;
                counter = counter + 1;
            end
            t3_Nums(subjnum,h).ReadyToPlot_Rand(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            
            idd_norm = t3_Nums(subjnum,h).ReadyToPlot_Rand(:,2) == t3_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T3(subjnum,h).xticklab_a{counter_n} = num2str(unique(t3_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,2)));
                t3_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,9) = t3_Nums(subjnum,h).All_normalized(t3_Nums(subjnum,h).sort_norm_ID(i));
                T3(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t3_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
            
        end
        
        
        t3_Nums(subjnum,h).SortedMeanIPI_All = t3_Nums(subjnum,h).MeanAll_IPI(t3_Nums(subjnum,h).sort_ID);
        t3_Nums(subjnum,h).SortedIPI_All = t3_Nums(subjnum,h).All_IPI(t3_Nums(subjnum,h).sort_ID);
        t3_Nums(subjnum,h).ReadyToPlot_All = cell2mat(t3_Nums(subjnum,h).SortedIPI_All);
        counter = 1;
        counter_n = 1;
        for i = 1:length(CMB.comb3)
            idd = t3_Nums(subjnum,h).ReadyToPlot_All(:,2) == t3_Nums(subjnum,h).sort_ID(i);
            if sum(idd)
                T3(subjnum,h).xticklab_a{counter} = num2str(unique(t3_Nums(subjnum,h).ReadyToPlot_All(idd,2)));
                t3_Nums(subjnum,h).ReadyToPlot_All(idd,7) = t3_Nums(subjnum,h).All(t3_Nums(subjnum,h).sort_ID(i));
                T3(subjnum,h).xtick_a(counter) = i;
                counter = counter + 1;
            end
            t3_Nums(subjnum,h).ReadyToPlot_All(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            
            idd_norm = t3_Nums(subjnum,h).ReadyToPlot_All(:,2) == t3_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T3(subjnum,h).xticklab_a{counter_n} = num2str(unique(t3_Nums(subjnum,h).ReadyToPlot_All(idd_norm,2)));
                t3_Nums(subjnum,h).ReadyToPlot_All(idd_norm,9) = t3_Nums(subjnum,h).All_normalized(t3_Nums(subjnum,h).sort_norm_ID(i));
                T3(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t3_Nums(subjnum,h).ReadyToPlot_All(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
        end
        
        
        h1 = figure;
        hold on
        [xcoordC_3{subjnum,h},PLOTC_3{subjnum,h},ERRORC_3{subjnum,h}]  = lineplot(t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,6) , t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,1) , 'subset' , ismember(t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,4) , [4 5]));
        [xcoordR_3{subjnum,h},PLOTR_3{subjnum,h},ERRORR_3{subjnum,h}]  = lineplot(t3_Nums(subjnum,h).ReadyToPlot_Rand(:,6) , t3_Nums(subjnum,h).ReadyToPlot_Rand(:,1) , 'subset' , ismember(t3_Nums(subjnum,h).ReadyToPlot_Rand(:,4) , [4 5]));
        [xcoordA_3{subjnum,h},PLOTA_3{subjnum,h},ERRORA_3{subjnum,h}]  = lineplot(t3_Nums(subjnum,h).ReadyToPlot_All(:,6) , t3_Nums(subjnum,h).ReadyToPlot_All(:,1) , 'subset' , ismember(t3_Nums(subjnum,h).ReadyToPlot_All(:,4) , [4 5]));
        close(h1)
        
        %
        %             anovan(t3_Nums(h).ReadyToPlot_Rand(:,1) , t3_Nums(h).ReadyToPlot_Rand(:,2))
        temp = corrcoef(t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,1) , t3_Nums(subjnum,h).ReadyToPlot_Chunked(:,7));
        C3_chunked(subjnum,h) = temp(2);
        temp  = corrcoef(t3_Nums(subjnum,h).ReadyToPlot_Rand(:,1) , t3_Nums(subjnum,h).ReadyToPlot_Rand(:,7));
        C3_random(subjnum,h) = temp(2);
        temp  = corrcoef(t3_Nums(subjnum,h).ReadyToPlot_All(:,1) , t3_Nums(subjnum,h).ReadyToPlot_All(:,7));
        C3_all(subjnum,h) = temp(2);
        
        
        
        clear xtick_r xtick_c xticklab_r xticklab_c
        allt4 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-3) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-3) length(ANA.AllPress)*(size(ANA.AllPress , 2)-3)];
        t4_Nums(subjnum,h).Chunked = zeros(length(CMB.comb4) , 1);
        t4_Nums(subjnum,h).Rand    = zeros(length(CMB.comb4) , 1);
        t4_Nums(subjnum,h).All     = zeros(length(CMB.comb4) , 1);
        
        ANA1.t4_Nums = zeros(length(ANA1.AllPress) , size(ANA1.AllPress , 2) -3);
        ANA0.t4_Nums = zeros(length(ANA0.AllPress) , size(ANA0.AllPress , 2) -3);
        ANA.t4_Nums  = zeros(length(ANA.AllPress) , size(ANA.AllPress , 2) -3);
        for t4 = 1:length(CMB.comb4)
            t4_Nums(subjnum,h).Chunked_IPI{t4,1} = [];
            t4_Nums(subjnum,h).Rand_IPI{t4,1} = [];
            t4_Nums(subjnum,h).All_IPI{t4,1} = [];
            
            t4_Nums(subjnum,h).TranNumb(t4 , 1) = t4;
            t4_Nums(subjnum,h).Transition(t4 , :) = CMB.comb4(t4,:);
            for p = 1:size(ANA1.AllPress , 2) -3
                t4_Nums(subjnum,h).Chunked(t4,1) =  t4_Nums(subjnum,h).Chunked(t4,1) + sum(ismember(ANA1.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                t4_Nums(subjnum,h).Rand(t4,1) =  t4_Nums(subjnum,h).Rand(t4,1) + sum(ismember(ANA0.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                t4_Nums(subjnum,h).All(t4,1) =  t4_Nums(subjnum,h).All(t4,1) + sum(ismember(ANA.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows'));
                
                ANA1.t4_Nums(ismember(ANA1.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                ANA0.t4_Nums(ismember(ANA0.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                ANA.t4_Nums(ismember(ANA.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') , p) = t4;
                
                CorID = ismember(ANA1.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') & ~sum(ANA1.isWrong(:,p:p+3) , 2);
                if ~LastIPI
                    t4_Nums(subjnum,h).Chunked_IPI{t4} = [t4_Nums(subjnum,h).Chunked_IPI{t4} ; [sum(ANA1.IPI(CorID , p:p+2),2) t4*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                else
                    t4_Nums(subjnum,h).Chunked_IPI{t4} = [t4_Nums(subjnum,h).Chunked_IPI{t4} ; [sum(ANA1.IPI(CorID , p+2),2) t4*ones(length(ANA1.IPI(CorID , p)) , 1) ANA1.SN(CorID) ANA1.Day(CorID) ANA1.Horizon(CorID)]];
                end
                
                CorID = ismember(ANA0.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') & ~sum(ANA0.isWrong(:,p:p+3) , 2);
                if ~LastIPI
                    t4_Nums(subjnum,h).Rand_IPI{t4} = [t4_Nums(subjnum,h).Rand_IPI{t4} ; [sum(ANA0.IPI(CorID , p:p+2),2) t4*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                else
                    t4_Nums(subjnum,h).Rand_IPI{t4} = [t4_Nums(subjnum,h).Rand_IPI{t4} ; [sum(ANA0.IPI(CorID , p+2),2) t4*ones(length(ANA0.IPI(CorID , p)) , 1) ANA0.SN(CorID) ANA0.Day(CorID) ANA0.Horizon(CorID)]];
                end
                
                CorID = ismember(ANA.AllPress(:,p:p+3) ,CMB.comb4(t4,:), 'rows') & ~sum(ANA.isWrong(:,p:p+3) , 2);
                if ~LastIPI
                    t4_Nums(subjnum,h).All_IPI{t4} = [t4_Nums(subjnum,h).All_IPI{t4} ; [sum(ANA.IPI(CorID , p:p+2),2) t4*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                else
                    t4_Nums(subjnum,h).All_IPI{t4} = [t4_Nums(subjnum,h).All_IPI{t4} ; [sum(ANA.IPI(CorID , p+2),2) t4*ones(length(ANA.IPI(CorID , p)) , 1) ANA.SN(CorID) ANA.Day(CorID) ANA.Horizon(CorID)]];
                end
            end
            t4_Nums(subjnum,h).MeanChunked_IPI(t4,1) = nanmean(t4_Nums(subjnum,h).Chunked_IPI{t4}(:,1));
            t4_Nums(subjnum,h).MeanRand_IPI(t4,1) = nanmean(t4_Nums(subjnum,h).Rand_IPI{t4}(:,1));
            t4_Nums(subjnum,h).MeanAll_IPI(t4,1) = nanmean(t4_Nums(subjnum,h).All_IPI{t4}(:,1));
        end
        t4_Nums(subjnum,h).Chunked = t4_Nums_allh(subjnum).Chunked;%t4_Nums(h).Chunked/allt4(1);
        t4_Nums(subjnum,h).Rand = t4_Nums_allh(subjnum).Rand;%t4_Nums(h).Rand/allt4(2);
        t4_Nums(subjnum,h).All = t4_Nums_allh(subjnum).All;%t4_Nums(h).Rand/allt4(2);
        t4_Nums(subjnum,h).All_normalized = t4_Nums_allh(subjnum).All_normalized;%t4_Nums(h).Rand/allt4(2);
        [~ , t4_Nums(subjnum,h).sort_ID] = sort(t4_Nums(subjnum,h).All , 'descend');
        [~ , t4_Nums(subjnum,h).sort_norm_ID] = sort(t4_Nums(subjnum,h).All_normalized , 'descend');
        
        t4_Nums(subjnum,h).SortedMeanIPI_Chunked = t4_Nums(subjnum,h).MeanChunked_IPI(t4_Nums(subjnum,h).sort_ID);
        t4_Nums(subjnum,h).SortedIPI_Chunked = t4_Nums(subjnum,h).Chunked_IPI(t4_Nums(subjnum,h).sort_ID);
        t4_Nums(subjnum,h).ReadyToPlot_Chunked = cell2mat(t4_Nums(subjnum,h).SortedIPI_Chunked);
        counter = 1;
        counter_n = 1;
        ChunkNum4 = [];
        for i = 1:length(CMB.comb4)
            idd = t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,2) == t4_Nums(subjnum,h).sort_ID(i);
            
            if sum(ismember(CMB.Chunks{4} , CMB.comb4(i,:)  , 'rows'))
                ChunkNum4 = [ChunkNum4 ; [i , find(t4_Nums(subjnum,h).sort_ID == i)]];
            end
            
            if sum(idd)
                T4(subjnum,h).xticklab_c{counter} = num2str(unique(t4_Nums(subjnum,h).ReadyToPlot_Chunked(idd,2)));
                t4_Nums(subjnum,h).ReadyToPlot_Chunked(idd,7) = t4_Nums(subjnum,h).All(t4_Nums(subjnum,h).sort_ID(i));
                T4(subjnum,h).xtick_c(counter) = i;
                counter = counter + 1;
            end
            idd_norm = t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,2) == t4_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T4(subjnum,h).xticklab_c{counter_n} = num2str(unique(t4_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,2)));
                t4_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,9) = t4_Nums(subjnum,h).All_normalized(t4_Nums(subjnum,h).sort_norm_ID(i));
                T4(subjnum,h).xtick_c(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t4_Nums(subjnum,h).ReadyToPlot_Chunked(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            t4_Nums(subjnum,h).ReadyToPlot_Chunked(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
            
        end
        
        t4_Nums(subjnum,h).SortedMeanIPI_Rand = t4_Nums(subjnum,h).MeanRand_IPI(t4_Nums(subjnum,h).sort_ID);
        t4_Nums(subjnum,h).SortedIPI_Rand = t4_Nums(subjnum,h).Rand_IPI(t4_Nums(subjnum,h).sort_ID);
        t4_Nums(subjnum,h).ReadyToPlot_Rand = cell2mat(t4_Nums(subjnum,h).SortedIPI_Rand);
        counter = 1;
        counter_n = 1;
        for i = 1:length(CMB.comb4)
            idd = t4_Nums(subjnum,h).ReadyToPlot_Rand(:,2) == t4_Nums(subjnum,h).sort_ID(i);
            
            if sum(idd)
                T4(subjnum,h).xticklab_r{counter} = num2str(unique(t4_Nums(subjnum,h).ReadyToPlot_Rand(idd,2)));
                t4_Nums(subjnum,h).ReadyToPlot_Rand(idd,7) = t4_Nums(subjnum,h).All(t4_Nums(subjnum,h).sort_ID(i));
                T4(subjnum,h).xtick_r(counter) = i;
                counter = counter + 1;
            end
            t4_Nums(subjnum,h).ReadyToPlot_Rand(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            idd_norm = t4_Nums(subjnum,h).ReadyToPlot_Rand(:,2) == t4_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd_norm)
                T4(subjnum,h).xticklab_r{counter_n} = num2str(unique(t4_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,2)));
                t4_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,9) = t4_Nums(subjnum,h).All_normalized(t4_Nums(subjnum,h).sort_norm_ID(i));
                T4(subjnum,h).xtick_r(counter_n) = i;
                counter_n = counter_n + 1;
            end
            
            t4_Nums(subjnum,h).ReadyToPlot_Rand(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
        end
        
        t4_Nums(subjnum,h).SortedMeanIPI_All = t4_Nums(subjnum,h).MeanAll_IPI(t4_Nums(subjnum,h).sort_ID);
        t4_Nums(subjnum,h).SortedIPI_All = t4_Nums(subjnum,h).All_IPI(t4_Nums(subjnum,h).sort_ID);
        t4_Nums(subjnum,h).ReadyToPlot_All = cell2mat(t4_Nums(subjnum,h).SortedIPI_All);
        counter = 1;
        counter_n = 1;
        for i = 1:length(CMB.comb4)
            idd = t4_Nums(subjnum,h).ReadyToPlot_All(:,2) == t4_Nums(subjnum,h).sort_ID(i);
            idd_norm = t4_Nums(subjnum,h).ReadyToPlot_All(:,2) == t4_Nums(subjnum,h).sort_norm_ID(i);
            if sum(idd)
                T4(subjnum,h).xticklab_a{counter} = num2str(unique(t4_Nums(subjnum,h).ReadyToPlot_All(idd,2)));
                t4_Nums(subjnum,h).ReadyToPlot_All(idd,7) = t4_Nums(subjnum,h).All(t4_Nums(subjnum,h).sort_ID(i));
                T4(subjnum,h).xtick_a(counter) = i;
                counter = counter + 1;
            end
            if sum(idd_norm)
                T4(subjnum,h).xticklab_a{counter_n} = num2str(unique(t4_Nums(subjnum,h).ReadyToPlot_All(idd,2)));
                t4_Nums(subjnum,h).ReadyToPlot_All(idd_norm,9) = t4_Nums(subjnum,h).All_normalized(t4_Nums(subjnum,h).sort_norm_ID(i));
                T4(subjnum,h).xtick_a(counter_n) = i;
                counter_n = counter_n + 1;
            end
            t4_Nums(subjnum,h).ReadyToPlot_All(idd,6) = i; %rank the triplet numbers based on their frequency of occurence
            t4_Nums(subjnum,h).ReadyToPlot_All(idd_norm,8) = i; %rank the triplet numbers based on their frequency of occurence
            
        end
        
        h1 = figure;
        hold on
        [xcoordC_4{subjnum,h},PLOTC_4{subjnum,h},ERRORC_4{subjnum,h}]  = lineplot(t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,6) , t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,1) , 'subset' , ismember(t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,4) , [4 5]));
        [xcoordR_4{subjnum,h},PLOTR_4{subjnum,h},ERRORR_4{subjnum,h}]  = lineplot(t4_Nums(subjnum,h).ReadyToPlot_Rand(:,6) , t4_Nums(subjnum,h).ReadyToPlot_Rand(:,1) , 'subset' , ismember(t4_Nums(subjnum,h).ReadyToPlot_Rand(:,4) , [4 5]));
        [xcoordA_4{subjnum,h},PLOTA_4{subjnum,h},ERRORA_4{subjnum,h}]  = lineplot(t4_Nums(subjnum,h).ReadyToPlot_All(:,6) , t4_Nums(subjnum,h).ReadyToPlot_All(:,1) , 'subset' , ismember(t4_Nums(subjnum,h).ReadyToPlot_All(:,4) , [4 5]));
        close(h1)
        
        %
        %             anovan(t4_Nums(h).ReadyToPlot_Rand(:,1) , t4_Nums(h).ReadyToPlot_Rand(:,2))
        temp = corrcoef(t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,1) , t4_Nums(subjnum,h).ReadyToPlot_Chunked(:,7));
        C4_chunked(subjnum,h) = temp(2);
        temp  = corrcoef(t4_Nums(subjnum,h).ReadyToPlot_Rand(:,1) , t4_Nums(subjnum,h).ReadyToPlot_Rand(:,7));
        C4_random(subjnum,h) = temp(2);
        temp  = corrcoef(t4_Nums(subjnum,h).ReadyToPlot_All(:,1) , t4_Nums(subjnum,h).ReadyToPlot_All(:,7));
        C4_all(subjnum,h) = temp(2);
    end
end

%%



C.IPI = [];
C.IPI_norm = [];
C.t2 = [];
C.SN =[];
C.BN =[];
C.Day = [];
C.Horizon  = [];
C.IPIarrangement = [];
C.IPIarrangement = [];
C.estIPIarrangement = [];
C.t2Rank = [];
C.t2Prob = [];
C.t2Rank_n = [];
C.t2Prob_n = [];
C.t3 = [];
C.t4 = [];
C.t3Rank = [];
C.t3Prob = [];
C.t4Rank = [];
C.t4Prob = [];
C.t3Rank_n = [];
C.t3Prob_n = [];
C.t4Rank_n = [];
C.t4Prob_n = [];
R = C;
All = C;
for subjnum = 1:length(subj_name)-1
    ANA1 = getrow(Dall ,ismember(Dall.seqNumb , [1:2]) & ismember(Dall.SN , subjnum) & ~Dall.isError);
    for tn = 1:length(ANA1.TN)
        n = (ANA1.AllPressIdx(tn , sum(~isnan(ANA1.AllPressIdx(tn , :))))  - ANA1.AllPressIdx(tn , 1)) / 1000;
        nIdx(tn , :) = (ANA1.AllPressIdx(tn , :) - ANA1.AllPressIdx(tn , 1))/n;
        ANA1.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
    end
    [a,d]=find(isnan(ANA1.IPI) | (ANA1.IPI> 2000));
    tns = ones(length(ANA1.TN),1);
    tns(a) = 0;
    ANA1 = getrow(ANA1 , logical(tns));
    ANA0 = getrow(Dall ,ismember(Dall.seqNumb , 0) & ismember(Dall.SN , subjnum)  & ~Dall.isError);
    for tn = 1:length(ANA0.TN)
        n = (ANA0.AllPressIdx(tn , sum(~isnan(ANA0.AllPressIdx(tn , :))))  - ANA0.AllPressIdx(tn , 1)) / 1000;
        nIdx(tn , :) = (ANA0.AllPressIdx(tn , :) - ANA0.AllPressIdx(tn , 1))/n;
        ANA0.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
    end
    [a,d]=find(isnan(ANA0.IPI) | (ANA0.IPI> 2000));
    tns = ones(length(ANA0.TN),1);
    tns(a) = 0;
    ANA0 = getrow(ANA0 , logical(tns));
    ANA  = addstruct(ANA1 , ANA0);
    
    
    
    
    
    allT2 = [length(ANA1.AllPress)*(size(ANA1.AllPress , 2)-1) length(ANA0.AllPress)*(size(ANA0.AllPress , 2)-1) length(ANA.AllPress)*(size(ANA.AllPress , 2)-1)];
    
    [t2sort_ID_n(:,1) , t2sort_ID_n(:,2)] = sort(t2_Nums_allh(subjnum).All_normalized , 'descend');
    [t3sort_ID_n(:,1) , t3sort_ID_n(:,2)] = sort(t3_Nums_allh(subjnum).All_normalized , 'descend');
    [t4sort_ID_n(:,1) , t4sort_ID_n(:,2)] = sort(t4_Nums_allh(subjnum).All_normalized , 'descend');
    [t2sort_ID(:,1) , t2sort_ID(:,2)] = sort(t2_Nums_allh(subjnum).All , 'descend');
    [t3sort_ID(:,1) , t3sort_ID(:,2)] = sort(t3_Nums_allh(subjnum).All , 'descend');
    [t4sort_ID(:,1) , t4sort_ID(:,2)] = sort(t4_Nums_allh(subjnum).All , 'descend');
    
    for p = 3:11
        for t2 = 1:length(CMB.comb2)
            
            t2id = find(t2sort_ID(:,2) == t2);
            t2id_n = find(t2sort_ID_n(:,2) == t2);
            
            clear t3  t4 t3Rank t4Rank t4Prob t3Prob t4Rank_n t3Rank_n t4Prob_n t3Prob_n
            CorID = ismember(ANA1.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows');
            tempID = find(CorID);
            if ~isempty(tempID)
                C.IPI = [C.IPI ; ANA1.IPI(CorID , p)];
                C.IPI_norm = [C.IPI_norm ; ANA1.IPI_norm(CorID , p)];
                C.t2 = [C.t2 ; t2*ones(length(ANA1.IPI(CorID , p)) , 1)];
                C.SN = [C.SN ; ANA1.SN(CorID)];
                C.BN = [C.BN ; ANA1.BN(CorID)];
                C.Day = [C.Day ;ANA1.Day(CorID)];
                C.Horizon  = [C.Horizon ;ANA1.Horizon(CorID)];
                C.IPIarrangement = [C.IPIarrangement ;ANA1.IPIarrangement(CorID , p)];
                C.estIPIarrangement = [C.estIPIarrangement ;ANA1.estChnkBndry(CorID , p)];
                C.t2Rank = [C.t2Rank ; t2id*ones(length(ANA1.IPI(CorID , p)) , 1)];
                C.t2Prob = [C.t2Prob ; t2sort_ID(t2id,1)*ones(length(ANA1.IPI(CorID , p)) , 1)];
                
                C.t2Rank_n = [C.t2Rank_n ; t2id_n*ones(length(ANA1.IPI(CorID , p)) , 1)];
                C.t2Prob_n = [C.t2Prob_n ; t2sort_ID_n(t2id_n,1)*ones(length(ANA1.IPI(CorID , p)) , 1)];
                
                for jj = 1:length(tempID)
                    t3(jj,:) = [find(ismember(CMB.comb3 , ANA1.AllPress(tempID(jj),p-1:p+1), 'rows'))   ,   find(ismember(CMB.comb3 , ANA1.AllPress(tempID(jj),p:p+2) ,'rows'))];
                    t3Rank(jj , :) = [find(t3sort_ID(:,2) == t3(jj,1)) , find(t3sort_ID(:,2) == t3(jj,2))];
                    t3Prob(jj,:) = [t3sort_ID(t3Rank(jj , 1),1) , t3sort_ID(t3Rank(jj , 2),1)];
                    
                    t3Rank_n(jj , :) = [find(t3sort_ID_n(:,2) == t3(jj,1)) , find(t3sort_ID_n(:,2) == t3(jj,2))];
                    t3Prob_n(jj,:) = [t3sort_ID_n(t3Rank_n(jj , 1),1) , t3sort_ID_n(t3Rank_n(jj , 2),1)];
                    
                    
                    t4(jj,:) = [find(ismember(CMB.comb4 , ANA1.AllPress(tempID(jj),p-2:p+1), 'rows'))   ,   find(ismember(CMB.comb4 , ANA1.AllPress(tempID(jj),p-1:p+2) ,'rows'))   ,   find(ismember(CMB.comb4 , ANA1.AllPress(tempID(jj),p:p+3) ,'rows'))];
                    t4Rank(jj , :) = [find(t4sort_ID(:,2) == t4(jj,1))  , find(t4sort_ID(:,2) == t4(jj,2))  , find(t4sort_ID(:,2) == t4(jj,3))];
                    t4Prob(jj,:) = [t4sort_ID(t4Rank(jj , 1),1) , t4sort_ID(t4Rank(jj , 2),1)  ,  t4sort_ID(t4Rank(jj , 3),1)];
                    
                    t4Rank_n(jj , :) = [find(t4sort_ID_n(:,2) == t4(jj,1))  , find(t4sort_ID_n(:,2) == t4(jj,2))  , find(t4sort_ID_n(:,2) == t4(jj,3))];
                    t4Prob_n(jj,:) = [t4sort_ID_n(t4Rank_n(jj , 1),1) , t4sort_ID_n(t4Rank_n(jj , 2),1)  ,  t4sort_ID_n(t4Rank_n(jj , 3),1)];
                end
                C.t3 = [C.t3 ; t3];
                C.t4 = [C.t4 ; t4];
                C.t3Rank = [C.t3Rank ; t3Rank];
                C.t3Prob = [C.t3Prob ; t3Prob];
                C.t4Rank = [C.t4Rank ; t4Rank];
                C.t4Prob = [C.t4Prob ; t4Prob];
                
                C.t3Rank_n = [C.t3Rank_n ; t3Rank_n];
                C.t3Prob_n = [C.t3Prob_n ; t3Prob_n];
                C.t4Rank_n = [C.t4Rank_n ; t4Rank_n];
                C.t4Prob_n = [C.t4Prob_n ; t4Prob_n];
            end
            
            
            
            
            clear t3  t4 t3Rank t4Rank t4Prob t3Prob t4Rank_n t3Rank_n t4Prob_n t3Prob_n
            CorID = ismember(ANA0.AllPress(:,p:p+1) ,CMB.comb2(t2,:), 'rows');
            tempID = find(CorID);
            if ~isempty(tempID)
                R.IPI = [R.IPI ; ANA0.IPI(CorID , p)];
                R.IPI_norm = [R.IPI_norm ; ANA0.IPI_norm(CorID , p)];
                R.t2 = [R.t2 ; t2*ones(length(ANA0.IPI(CorID , p)) , 1)];
                R.SN = [R.SN ; ANA0.SN(CorID)];
                R.BN  = [R.BN ; ANA0.BN(CorID)];
                R.Day = [R.Day ;ANA0.Day(CorID)];
                R.Horizon  = [R.Horizon; ANA0.Horizon(CorID)];
                R.estIPIarrangement = [R.estIPIarrangement ;ANA0.estChnkBndry(CorID , p)];
                R.t2Rank = [R.t2Rank ; t2id*ones(length(ANA0.IPI(CorID , p)) , 1)];
                R.t2Prob = [R.t2Prob ; t2sort_ID(t2id,1)*ones(length(ANA0.IPI(CorID , p)) , 1)];
                R.t2Rank_n = [R.t2Rank_n ; t2id_n*ones(length(ANA0.IPI(CorID , p)) , 1)];
                R.t2Prob_n = [R.t2Prob_n ; t2sort_ID_n(t2id_n,1)*ones(length(ANA0.IPI(CorID , p)) , 1)];
                for jj = 1:length(tempID)
                    % find chunks in raondom sequences and give them chunk arrangements
                    CorID_findChunks = ismember(ANA0.AllPress(:,p-1:p+1) , CMB.Chunks{3}(:,1:3) , 'rows');
                    ANA0.IPIarrangement(CorID_findChunks,p-1:p) = repmat([2 2] , sum(CorID_findChunks),1);
                    t3(jj,:) = [find(ismember(CMB.comb3 , ANA0.AllPress(tempID(jj),p-1:p+1), 'rows'))   ,   find(ismember(CMB.comb3 , ANA0.AllPress(tempID(jj),p:p+2) ,'rows'))];
                    t3Rank(jj , :) = [find(t3sort_ID(:,2) == t3(jj,1)) , find(t3sort_ID(:,2) == t3(jj,2))];
                    t3Prob(jj,:) = [t3sort_ID(t3Rank(jj , 1),1) , t3sort_ID(t3Rank(jj , 2),1)];
                    
                    t3Rank_n(jj , :) = [find(t3sort_ID_n(:,2) == t3(jj,1)) , find(t3sort_ID_n(:,2) == t3(jj,2))];
                    t3Prob_n(jj,:) = [t3sort_ID_n(t3Rank_n(jj , 1),1) , t3sort_ID_n(t3Rank_n(jj , 2),1)];
                    
                    % find chunks in raondom sequences and give them chunk arrangements
                    CorID_findChunks = ismember(ANA0.AllPress(:,p-2:p+1) , CMB.Chunks{3},  'rows');
                    ANA0.IPIarrangement(CorID_findChunks,p-2:p) = repmat([2 2 2] , sum(CorID_findChunks),1);
                    t4(jj,:) = [find(ismember(CMB.comb4 , ANA0.AllPress(tempID(jj),p-2:p+1), 'rows'))   ,   find(ismember(CMB.comb4 , ANA0.AllPress(tempID(jj),p-1:p+2) ,'rows'))   ,   find(ismember(CMB.comb4 , ANA0.AllPress(tempID(jj),p:p+3) ,'rows'))];
                    t4Rank(jj , :) = [find(t4sort_ID(:,2) == t4(jj,1))  , find(t4sort_ID(:,2) == t4(jj,2))  , find(t4sort_ID(:,2) == t4(jj,3))];
                    t4Prob(jj,:) = [t4sort_ID(t4Rank(jj , 1),1) , t4sort_ID(t4Rank(jj , 2),1)  ,  t4sort_ID(t4Rank(jj , 3),1)];
                    
                    t4Rank_n(jj , :) = [find(t4sort_ID_n(:,2) == t4(jj,1))  , find(t4sort_ID_n(:,2) == t4(jj,2))  , find(t4sort_ID_n(:,2) == t4(jj,3))];
                    t4Prob_n(jj,:) = [t4sort_ID_n(t4Rank_n(jj , 1),1) , t4sort_ID_n(t4Rank_n(jj , 2),1)  ,  t4sort_ID_n(t4Rank_n(jj , 3),1)];
                end
                R.IPIarrangement = [R.IPIarrangement ;ANA0.IPIarrangement(CorID , p)];
                R.t3 = [R.t3 ; t3];
                R.t4 = [R.t4 ; t4];
                R.t3Rank = [R.t3Rank ; t3Rank];
                R.t3Prob = [R.t3Prob ; t3Prob];
                R.t4Rank = [R.t4Rank ; t4Rank];
                R.t4Prob = [R.t4Prob ; t4Prob];
                
                R.t3Rank_n = [R.t3Rank_n ; t3Rank_n];
                R.t3Prob_n = [R.t3Prob_n ; t3Prob_n];
                R.t4Rank_n = [R.t4Rank_n ; t4Rank_n];
                R.t4Prob_n = [R.t4Prob_n ; t4Prob_n];
            end
            
        end
        
    end
end
C.seqNumb(1:length(C.SN) , :) =1;
R.seqNumb(1:length(R.SN) , :) =0;
All = addstruct(C,R);

%% bin the probabilities into 5 classes of probability and test the effect of probability on IPIs in Random sequences
C.t2Rank_n_binned = C.t2Rank_n;
rr = 1;
for j = 0:5:25
    C.t2Rank_n_binned(C.t2Rank_n_binned>=j & C.t2Rank_n_binned<j+5) = rr;
    rr = rr+1;
end
R.t2Rank_n_binned = R.t2Rank_n;
rr = 1;
for j = 0:6:25
    R.t2Rank_n_binned(R.t2Rank_n_binned>=j & R.t2Rank_n_binned<j+6) = rr;
    rr = rr+1;
end

C.t3Rank_n_binned = C.t3Rank_n(:,1);
rr = 1;
for j = 0:20:125
    C.t3Rank_n_binned(C.t3Rank_n_binned>=j & C.t3Rank_n_binned<j+20) = rr;
    rr = rr+1;
end
R.t3Rank_n_binned = R.t3Rank_n(:,1);
rr = 1;
for j = 0:25:125
    R.t3Rank_n_binned(R.t3Rank_n_binned>=j & R.t3Rank_n_binned<j+25) = rr;
    rr = rr+1;
end

C.t4Rank_n_binned = C.t4Rank_n(:,1);
rr = 1;
for j = 0:60:625
    C.t4Rank_n_binned(C.t4Rank_n_binned>=j & C.t4Rank_n_binned<j+60) = rr;
    rr = rr+1;
end
R.t4Rank_n_binned = R.t4Rank_n(:,1);
rr = 1;
for j = 0:80:625
    R.t4Rank_n_binned(R.t4Rank_n_binned>=j & R.t4Rank_n_binned<j+80) = rr;
    rr = rr+1;
end
% map the block to less / half per day / so bin every 4 blocks to 1

All = addstruct(C,R);
dd = unique(All.Day);
for db= 1:length(dd)
    T = getrow(All , All.Day == dd(db));
    bls = unique(T.BN);
    id1 = ismember(All.BN , bls(1:floor(length(bls)/2)));
    id2 = ismember(All.BN , bls(floor(length(bls)/2):end));
    All.BN(id1) = 2*db -1;
    All.BN(id2) = 2*db;
end
CC = getrow(All,All.seqNumb == 1);
RR = getrow(All,All.seqNumb == 0);

C_Summarized  = tapply(CC , {'Horizon' ,'SN' , 'IPIarrangement' , 'BN' ,'t2Rank_n_binned' , 't3Rank_n_binned' , 't4Rank_n_binned'} , {'IPI' , 'nanmedian(x)'} , {'IPI_norm' , 'nanmedian(x)'} , {'Day' , 'nanmedian(x)'});
R_Summarized  = tapply(RR , {'Horizon' ,'SN' , 'IPIarrangement' , 'BN' ,'t2Rank_n_binned' , 't3Rank_n_binned' , 't4Rank_n_binned'} , {'IPI' , 'nanmedian(x)'} , {'IPI_norm' , 'nanmedian(x)'} , {'Day' , 'nanmedian(x)'});



save([baseDir , '/se2_TranProb.mat'])
