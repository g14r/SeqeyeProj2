function out  = se2_visualize(Dall , subjnum, what, distance, calc , day ,rep, GroupCode)
%%  distances:
% 'euclidean'	      Euclidean distance (default).
% 'squaredeuclidean'  Squared Euclidean distance. (This option is provided for efficiency only. It does not satisfy the triangle inequality.)
% 'seuclidean'	      Standardized Euclidean distance. Each coordinate difference between rows in X and Y is scaled by dividing by the corresponding element of the standard deviation computed from X, S=nanstd(X). To specify another value for S, use D = PDIST2(X,Y,'seuclidean',S).
% 'cityblock'	      City block metric.
% 'minkowski'         Minkowski distance. The default exponent is 2. To compute the distance with a different exponent, use D = pdist2(X,Y,'minkowski',P), where the exponent P is a scalar positive value.
% 'Chebychev'	      Chebychev distance (maximum coordinate difference).
% 'mahalanobis'	      Mahalanobis distance, using the sample covariance of X as computed by nancov. To compute the distance with a different covariance, use D = pdist2(X,Y,'mahalanobis',C) where the matrix C is symmetric and positive definite.
% 'cosine'	          One minus the cosine of the included angle between points (treated as vectors).
% 'correlation'	      One minus the sample correlation between points (treated as sequences of values).
% 'spearman'	      One minus the sample Spearman's rank correlation between observations, treated as sequences of values.
% 'hamming'           Hamming distance, the percentage of coordinates that differ.
% 'jaccard'           One minus the Jaccard coefficient, the percentage of nonzero coordinates that differ.%     case 'chunk_est_instance'
%%  Cases
%     case 'IPI_ttest_rand'
%     case 'MT_ttest_rand'
%     case 'MT_Vs_Horizon'
%     case 'IPI_Vs_Horizon'
%     case 'RT_Vs_Horizon'
%     case 'Points'
%     case 'MT_asymptote'
%     case 'test_MT_asymptote'
%     case 'eyepress_pos_traces'
%     case 'eyepress_vel_traces'
%     case 'crossvaldist_pos'
%     case 'crossvaldist_vel'
%     case 'crossval_IPI_dist'
%     case 'Errors'
%     case 'Sacc_Chunked'
%     case 'Sacc_All'
%     case 'transitions_All'
%     case 'glm_IPIs'
%     case 'crossvaldist_chunk'
%%
prefix = 'se1_';
% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye2/analyze';          %iMac
subj_name = {'AT1' , 'CG1' , 'HB1' , 'JT1' , 'CB1' , 'YM1' , 'NL1' , 'SR1' , 'IB1' , 'MZ1' , 'DW1', 'All'};

% subj_name = {'AT1' , 'CG1' , 'HB1' , 'JT1' , 'CB1' , 'YM1' , 'NL1' , 'SR1' , 'All'};
colors = [0.0899988149868883,0.789028923313949,0.710408685278170;0.320941032647761,0.317833053726229,0.311859945147533;0.511408938819178,0.452207453762982,0.291457127647727;0.0606063665682423,0.752227970049942,0.850357337374621;0.725687923545844,0.109861705750686,0.911647424007853;0.556555748561992,0.109742368593904,0.639276147276064;0.529359902481257,0.269883663704401,0.255370297944443;0.829982432033195,0.524637345396311,0.0886658400322831;0.858759034071804,0.972651076977497,0.838255587537226];
colors = [0 0 1;...
    0 1 0;...
    1 0 0;...
    0 1 1;...
    1 0 1;...
    1 0.69 0.39;...
    0.6 0.2 0;...
    0 0.75 0.75;...
    0.22 0.44 0.34;...
    0.32 0.19 0.19];
if subjnum == length(subj_name)
    %     subjnum = [1 3:length(subj_name)-1];
    subjnum = 1:length(subj_name)-1;
end


% load([baseDir , '/CMB.mat'])
%load([baseDir , '/se1_all.mat'])
days  = {1 ,2 ,3 ,4 ,5,[1:5] ,[2:5] [2:3] [4:5],[3:5]};
% subjnum = 2;

tid = {[1:250] [251 :750] [751:1000] [1:1000]};
switch what
    case 'IPI_ttest_rand'
        
        %% within and between IPI ttest with randomization
        %% IPI per day
        nn= input('norm (1) or non-norm (2)?');
        structNumb = [1 2];
        %         plotfcn = input('nanmean or nanmedian?' , 's');
        %% IPIs vs horizon
        % this is the output of the case: 'transitions_All' that is saved to disc
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 , structNumb])  & ~Dall.isError & ismember(Dall.Day , days{day}));
        
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        ANA.seqNumb(ANA.seqNumb>1) = 1;
        for tn  = 1:length(ANA.TN)
            ANA.ChunkBndry(tn , :) = diff(ANA.ChnkArrang(tn,:));
            a = find(ANA.ChunkBndry(tn , :));
            ANA.ChunkBndry(tn , a-1) = 3;
            ANA.ChunkBndry(tn , end) = 3;
            ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
            ANA.ChunkBndry(tn , 1:2) = [-1 -1];  % dont account for the first and last sseqeuce presses
            ANA.ChunkBndry(tn , end-1:end) = [-2 -2];% dont account for the first and last sseqeuce presses
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,13);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,13);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,13);
            ANA.IPI_prsnumb(tn , :) = [1 :13];
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,13);
        end
        
        newANA.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        newANA.IPI_norm = reshape(ANA.IPI_norm , numel(ANA.IPI) , 1);
        newANA.ChunkBndry = reshape(ANA.ChunkBndry , numel(ANA.IPI) , 1);
        newANA.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        newANA.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        newANA.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        newANA.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        newANA.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        newANA.ChunkBndry(newANA.ChunkBndry>2) = 2;
        IPItable  = tapply(newANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb' , 'ChunkBndry' , 'prsnumb'} , {'IPI' , 'nanmedian(x)'},{'IPI_norm' , 'nanmedian(x)'});
        
        figure('color' , 'white')
        Days = unique(IPItable.Day);
        hrz = unique(IPItable.Horizon);
        % gives a scatter plot of the normalized vs non-normalized IPIs in h = 1
        for d = 1:length(Days)
            subplot(2,3,d)
            id = IPItable.seqNumb == 1 & IPItable.ChunkBndry==1 & IPItable.Day==Days(d) & ismember(IPItable.Horizon , 1);
            plot(IPItable.IPI(id) , IPItable.IPI_norm(id) , '*')
            hold on
            id = IPItable.seqNumb == 1 & IPItable.ChunkBndry==2 & IPItable.Day==Days(d) & ismember(IPItable.Horizon , 1);
            plot(IPItable.IPI(id) , IPItable.IPI_norm(id) , '*')
            xlabel('msec')
            ylabel('Norm time')
            legend({'Between Chunk' , 'Within Chunk'})
            title(['Day ' , num2str(Days(d))])
            set(gca , 'FontSize' , 20 , 'XLim' , [50 1200] , 'YLim' , [10 180]);
            grid on
        end
        
        
        if nn == 2
            ylim  = [200 610];
        elseif nn == 1
            IPItable.IPI  = IPItable.IPI_norm;
            ylim  = [70 95];
        end
        
        h1 = figure;
        for h  =[1:8 , 13]
            for d = 1:5
                [cooIPI_sb(:,h) ,pltIPI_sb(:,h) , errIPI_sb(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 1) & ismember(IPItable.Horizon,h) &  ~ismember(IPItable.ChunkBndry, [-1 , -2]) &  ismember(IPItable.ChunkBndry, [1]));
                [cooIPI_sw(:,h) ,pltIPI_sw(:,h) , errIPI_sw(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 1) & ismember(IPItable.Horizon,h) &  ~ismember(IPItable.ChunkBndry, [-1 , -2]) &  ismember(IPItable.ChunkBndry, [2]));
                [cooIPI_r(:,h) ,pltIPI_r(:,h) , errIPI_r(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 0) & ismember(IPItable.Horizon,h) &  ~ismember(IPItable.ChunkBndry, [-1 , -2]));
            end
        end
        close(h1)
        figure('color', 'white')
        hCounter = 1;
        for h  =[1:8 , 13]
            subplot(3,3,hCounter)
            errorbar(cooIPI_sb(:,h) ,pltIPI_sb(:,h) , errIPI_sb(:,h) , 'LineWidth' , 3)
            hold on
            errorbar(cooIPI_sw(:,h) ,pltIPI_sw(:,h) , errIPI_sw(:,h), 'LineWidth' , 3)
            errorbar(cooIPI_r(:,h) ,pltIPI_r(:,h) , errIPI_r(:,h), 'LineWidth' , 3)
            set(gca , 'FontSize' ,20)
            title(['IPI, H = ' , num2str(h)])
            ylabel('msec' )
            xlabel('Days')
            set(gca , 'YLim' , ylim ,'XTick' , [1:5], 'FontSize' , 20)
            grid on
            hCounter = hCounter +  1;
            
        end
        legend({'Between' 'Within' 'Random'})
        
        
        hCounter = 1;
        fInd = {[1] [2] [3] [4:6]}; % figure index
        horz = {[1] [2] [3] [4:13]};
        HTile = {'1' , '2' , '3' , '4 - Full'};
        h1 = figure;
        for h  = 1:length(horz)
            for d = 1:5
                [cooIPI_sb1(:,h) ,pltIPI_sb1(:,h) , errIPI_sb1(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 1) & ismember(IPItable.Horizon,horz{h}) &  ~ismember(IPItable.ChunkBndry, [-1 , -2]) &  ismember(IPItable.ChunkBndry, [1]));
                [cooIPI_sw1(:,h) ,pltIPI_sw1(:,h) , errIPI_sw1(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 1) & ismember(IPItable.Horizon,horz{h}) &  ~ismember(IPItable.ChunkBndry, [-1 , -2]) &  ismember(IPItable.ChunkBndry, [2]));
                [cooIPI_r1(:,h) ,pltIPI_r1(:,h) , errIPI_r1(:,h)] = lineplot(IPItable.Day, IPItable.IPI, 'subset' , ismember(IPItable.seqNumb , 0) & ismember(IPItable.Horizon,horz{h}) &  ~ismember(IPItable.ChunkBndry, [-1 , -2]));
            end
        end
        close(h1)
        figure('color', 'white')
        for h  = 1:length(horz)
            subplot(2,3,h)
            errorbar(cooIPI_sb1(:,h) ,pltIPI_sb1(:,h) , errIPI_sb1(:,h) , 'LineWidth' , 3)
            hold on
            errorbar(cooIPI_sw1(:,h) ,pltIPI_sw1(:,h) , errIPI_sw1(:,h), 'LineWidth' , 3)
            errorbar(cooIPI_r1(:,h) ,pltIPI_r1(:,h) , errIPI_r1(:,h), 'LineWidth' , 3)
            set(gca , 'FontSize' ,20 , 'GridAlpha' , 1 , 'Box' , 'off')
            title(['Iner-Press Intervals , H = ' , HTile{h}])
            ylabel('msec' )
            xlabel('Days')
            set(gca , 'YLim' , ylim ,'XTick' , [1:5], 'FontSize' , 20 , 'Box' , 'off')
            grid on
        end
        legend({'Between' 'Within' 'Random'})
        
        figure('color', 'white')
        hCounter = 1;
        subplot(131)
        for h  =[1:8 , 13]
            errorbar(cooIPI_sb(:,h) ,pltIPI_sb(:,h) , errIPI_sb(:,h) , 'LineWidth' , 3 , 'color', colors(hCounter ,:))
            hold on
            title(['Between chunk IPIs'])
            ylabel('msec' )
            xlabel('Days')
            set(gca , 'YLim' , ylim ,'XTick' , [1:5], 'FontSize' , 20 , 'Box' , 'off')
            grid on
            hCounter = hCounter +  1;
        end
        %         legend({'H=1' 'H=2' 'H=3' 'H=4' 'H=5' 'H=6' 'H=7' 'H=8' 'H=13'})
        
        subplot(132)
        hCounter = 1;
        for h  =[1:8 , 13]
            errorbar(cooIPI_sw(:,h) ,pltIPI_sw(:,h) , errIPI_sw(:,h), 'LineWidth' , 3, 'color', colors(hCounter ,:))
            hold on
            title(['Within Chunk IPIs'])
            ylabel('msec' )
            xlabel('Days')
            set(gca , 'YLim' , ylim ,'XTick' , [1:5], 'FontSize' , 20 , 'Box' , 'off')
            grid on
            hCounter = hCounter +  1;
        end
        %         legend({'H=1' 'H=2' 'H=3' 'H=4' 'H=5' 'H=6' 'H=7' 'H=8' 'H=13'})
        
        subplot(133)
        hCounter = 1;
        for h  =[1:8 , 13]
            errorbar(cooIPI_r(:,h) ,pltIPI_r(:,h) , errIPI_r(:,h), 'LineWidth' , 3, 'color', colors(hCounter ,:))
            hold on
            title(['Random IPIs'])
            ylabel('msec' )
            xlabel('Days')
            set(gca , 'YLim' , ylim ,'XTick' , [1:5], 'FontSize' , 20 , 'Box' , 'off')
            grid on
            hCounter = hCounter +  1;
        end
        legend({'H=1' 'H=2' 'H=3' 'H=4' 'H=5' 'H=6' 'H=7' 'H=8' 'H=13'})
        
        for d = 1:5
            hcount = 1;
            for h  = [1:8 13]
                temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.ChunkBndry,{'within/between'},'intercept',1 ,'subset' , ismember(IPItable.Horizon , [h]) & ismember(IPItable.seqNumb , 1) & ismember(IPItable.Day , d))  ;
                out.chunkEffect_Chunkedseq(hcount , d) = temp.eff(2).p;
                temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.ChunkBndry,{'within/between'},'intercept',1 ,'subset' , ismember(IPItable.Horizon , [h]) & ismember(IPItable.seqNumb , 0) & ismember(IPItable.Day , d))  ;
                out.chunkEffect_Randomseq(hcount , d) = temp.eff(2).p;
                hcount = hcount+1;
            end
            temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',[IPItable.ChunkBndry IPItable.Horizon] ,{'within/between' , 'horizon'},'intercept',1 ,'subset'  & ismember(IPItable.seqNumb , 1) & ismember(IPItable.Day , d))  ;
            out.chunkHorizonInterac_Chunkedseq(d) = temp.eff(4).p;
            temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',[IPItable.ChunkBndry IPItable.Horizon] ,{'within/between' , 'horizon'},'intercept',1 ,'subset'  & ismember(IPItable.seqNumb , 0) & ismember(IPItable.Day , d))  ;
            out.chunkHorizonInterac_Randseq(d) = temp.eff(4).p;
        end
        hcount = 1;
        for h  = [1:8 13]
            temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.Day,{'day'},'intercept',1 ,'subset' , ismember(IPItable.Horizon , [h]) & ismember(IPItable.seqNumb , 1) & ismember(IPItable.ChunkBndry , 1))  ;
            out.dayEffect_bet_Chunkedseq(hcount) = temp.eff(2).p;
            temp = anovaMixed(IPItable.IPI , IPItable.SN,'within',IPItable.Day,{'day'},'intercept',1 ,'subset' , ismember(IPItable.Horizon , [h]) & ismember(IPItable.seqNumb , 1) & ismember(IPItable.ChunkBndry , 2))  ;
            out.dayEffect_wit_Chunkedseq(hcount) = temp.eff(2).p;
            hcount = hcount+1;
        end
    case 'MT_ttest_rand'
        ANA1_all = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [1:2])  & ~Dall.isError & ismember(Dall.Day , days{day}));
        ANA0_all = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0])  & ~Dall.isError & ismember(Dall.Day , days{day}));
        allMT = [];
        figure('color' , 'white')
        h_counter = 1;
        for hrzn = [8:-1:1]
            ANA1 = getrow(ANA1_all , ANA1_all.Horizon == hrzn);
            ANA0 = getrow(ANA0_all , ANA0_all.Horizon == hrzn);
            
            ANA1.MT = ANA1.AllPressTimes(:,end-2) - ANA1.AllPressTimes(:,3);
            ANA1 = getrow(ANA1 , ANA1.MT <= 9000);
            ANA0.MT = ANA0.AllPressTimes(:,end-2) - ANA0.AllPressTimes(:,3);
            ANA0 = getrow(ANA0 , ANA0.MT <= 9000);
            
            out.MT(hrzn) = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
            %             out.MT(hrzn) = anovan([ANA1.MT ; ANA0.MT] , [0*ANA1.MT ; 1+0*ANA0.MT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
            h1 = figure('color' , 'white');
            [xcoord,PLOTs,ERRORs] = lineplot(ANA1.BN,ANA1.MT);
            hold on
            [xcoord,PLOTr,ERRORr] = lineplot(ANA0.BN,ANA0.MT);
            close(h1);
            
            subplot(3,4,h_counter)
            errorbar(PLOTs , ERRORs,'LineWidth' , 3)
            hold on
            errorbar(PLOTr , ERRORr,'LineWidth' , 3)
            hold on
            ax = gca;
            ax.FontSize = 20;
            title(['MT, H = ' , num2str(hrzn) , ' p = ' , num2str(out.MT(hrzn).eff(2).p)])
            ylabel('msec'  ,'FontSize' , 20)
            %             ax.XTick = [1:4];
            %             ax.YLim = [3000 10000];
            legend({'Chunked' , 'Random'})
            grid on
            h_counter = h_counter +  1;
            allMT = [allMT; [ANA1.MT ; ANA0.MT] [0*ANA1.MT ; 1+0*ANA0.MT] [ANA1.SN ;ANA0.SN] [ANA1.Horizon ;ANA0.Horizon]];
        end
        hrzn = 13;
        ANA1 = getrow(ANA1_all , ANA1_all.Horizon == hrzn);
        ANA0 = getrow(ANA0_all , ANA0_all.Horizon == hrzn);
        
        ANA1.MT = ANA1.AllPressTimes(:,end-2) - ANA1.AllPressTimes(:,3);
        ANA1 = getrow(ANA1 , ANA1.MT <= 9000);
        ANA0.MT = ANA0.AllPressTimes(:,end-2) - ANA0.AllPressTimes(:,3);
        ANA0 = getrow(ANA0 , ANA0.MT <= 9000);
        
        allMT = [allMT; [ANA1.MT ; ANA0.MT] [0*ANA1.MT ; 1+0*ANA0.MT] [ANA1.SN ;ANA0.SN] [ANA1.Horizon ;ANA0.Horizon]];
        out.MT(hrzn) = anovaMixed([ANA1.MT ; ANA0.MT] , [ANA1.SN ;ANA0.SN],'within',[0*ANA1.MT ; 1+0*ANA0.MT],{'Random/Chunked'},'intercept',1)  ;
        %         out.MT(hrzn) = anovan([ANA1.MT ; ANA0.MT] , [0*ANA1.MT ; 1+0*ANA0.MT] , 'display' , 'off' , 'varnames' , 'Rand/Chunked');
        %         out.allMT = anovan(allMT(:,1) , [allMT(:,2)  allMT(:,3)] , 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Rand/Chunked' , 'Horizon'}  );
        out.allMT = anovaMixed(allMT(:,1) , allMT(:,3) , 'within' , allMT(:,[2 , 4]) , {'Random/Chunked' , 'Horizon'},'intercept',1);
        disp(['The interaction effect between Chunked/Random , horizon on MTs is ' num2str(out.allMT(end).eff(3).p)])
        h1 = figure('color' , 'white');
        [xcoords,PLOTs,ERRORs] = lineplot(ANA1.BN,ANA1.MT);
        hold on
        [xcoordr,PLOTr,ERRORr] = lineplot(ANA0.BN,ANA0.MT);
        close(h1);
        
        subplot(3,4,[9:12])
        errorbar(xcoords,PLOTs , ERRORs,'LineWidth' , 3)
        hold on
        errorbar(xcoordr,PLOTr , ERRORr,'LineWidth' , 3)
        hold on
        ax = gca;
        ax.FontSize = 20;
        title(['MT, H = ' , num2str(hrzn) , ' p = ' , num2str(out.MT(hrzn).eff(2).p)])
        ylabel('msec'  ,'FontSize' , 20)
        %         ax.YLim = [3000 10000];
        legend({'Chunked' , 'Random'})
        grid on
    case 'MT_Vs_Horizon'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError);
        ANA.seqNumb(ANA.seqNumb == 2) = 1;
        MT  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb','BN'} , {'MT' , 'nanmedian(x)'});
        
        MT = getrow(MT , MT.MT <= 9000 );
        dayz = {1 [2 3] [4 5]};
        for d=  1:length(dayz)
            hc = 1;
            for h = [1:8 13]
                temp = anovaMixed(MT.MT  , MT.SN ,'within',MT.seqNumb ,{'Random/Chunked'},'intercept',1 ,...
                    'subset' , ismember(MT.Horizon , h) & ismember(MT.Day , d))  ;
                out.ChVsRand(hc,d)  = temp.eff(2).p;
                hc = hc+1;
            end
            h1 = figure('color' , 'white');
            [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(MT.Horizon,MT.MT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 1 & ismember(MT.Day , dayz{d}));
            hold on
            [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(MT.Horizon,MT.MT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 0 & ismember(MT.Day , dayz{d}));
            close(h1);
        end
        
        
        figure('color' , 'white');
        for d=  1:length(dayz)
            temp = anovaMixed(MT.MT  , MT.SN ,'within',[MT.Horizon MT.seqNumb] ,{'Horizon' 'Random/Chunked'},'intercept',1 ,...
                'subset' , ismember(MT.Day , d))  ;
            subplot(1,3,d)
            h1 = plotshade(xcoords{d}',PLOTs{d} , ERRORs{d},'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
            hold on
            h2 = plotshade(xcoordr{d}',PLOTr{d} , ERRORr{d},'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':');
            set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                'GridAlpha' , 1 , 'Box' , 'off');
            title(['Execution time - Day(s) ' , num2str(dayz{d})])
            ylabel('msec' )
            xlabel('Horizon' )
            legend([h1 h2] ,{'Chunked' , 'Random'})
            grid on
        end
        
        
        
        %%
        
        h1 = figure('color' , 'white');
        hold on
        
        hc = 1;
        for h  = [1:7]
            [xcoords_med{h},PLOTs_med{h},ERRORs_med{h}] = lineplot([MT.seqNumb MT.Day],MT.MT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(MT.Horizon , [h]) & ismember(MT.seqNumb , 1));
            [xcoordr_med{h},PLOTr_med{h},ERRORr_med{h}] = lineplot([MT.seqNumb MT.Day],MT.MT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(MT.Horizon , [h]) & ismember(MT.seqNumb , 0));
            for d=  1:length(dayz)
                [h d]
                temp = anovaMixed(MT.MT , MT.SN,'within',MT.Horizon,{'Horizon'},'intercept',1 ,'subset' ,...
                    ~ismember(MT.Horizon , [1:h]) & ismember(MT.seqNumb , 1) & ismember(MT.Day , dayz{d}))  ;
                out.HorizonEffect_Chunkedseq_not(h,d) = temp.eff(2).p;
                temp  = anovaMixed(MT.MT , MT.SN,'within',MT.Horizon,{'Horizon'},'intercept',1 ,'subset' ,...
                    ~ismember(MT.Horizon , [1:h]) & ismember(MT.seqNumb , 0) & ismember(MT.Day , dayz{d}))  ;
                out.HorizonEffect_Randomseq_not(h,d) = temp.eff(2).p;
            end
        end
        h = h +1;
        for d=  1:length(dayz)
            temp= anovaMixed(MT.MT , MT.SN,'within',MT.Horizon,{'Horizon'},'intercept',1,'subset' , ismember(MT.seqNumb , 1) & ismember(MT.Day , dayz{d}))  ;
            out.HorizonEffect_Chunkedseq_not(h,d) = temp.eff(2).p;
            temp = anovaMixed(MT.MT , MT.SN,'within',MT.Horizon,{'Horizon'},'intercept',1,'subset' , ismember(MT.seqNumb , 0) & ismember(MT.Day , dayz{d}))  ;
            out.HorizonEffect_Randomseq_not(h,d) = temp.eff(2).p;
        end
        h = 1;
        for hc  = [1:8 13]
            [xcoords_med{h},PLOTs_med{h},ERRORs_med{h}] = lineplot([MT.seqNumb MT.Day],MT.MT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(MT.Horizon , [hc]) & ismember(MT.seqNumb , 1));
            [xcoordr_med{h},PLOTr_med{h},ERRORr_med{h}] = lineplot([MT.seqNumb MT.Day],MT.MT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(MT.Horizon , [hc]) & ismember(MT.seqNumb , 0));
            h  = h + 1;
        end
        
        figure('color' , 'white');
        subplot(221)
        imagesc(cell2mat(PLOTs_med') , [2500 7000])
        set(gca , 'XTick' , [1:5] , 'YTick', [1:9] , 'YTickLabels' , [1 :8 , 13],'FontSize' , 20)
        title('Execution Time in  Chunked Sequences')
        ylabel('Horizon Size')
        xlabel('Days')
        
        subplot(222)
        imagesc(cell2mat(PLOTr_med'), [2500 7000])
        colorbar
        set(gca , 'XTick' , [1:5] , 'YTick', [1:9] , 'YTickLabels' , [1 :8 , 13],'FontSize' , 20)
        title('Execution Time in  Random Sequences')
        ylabel('Horizon Size')
        xlabel('Days')
        
        subplot(223)
        A = cell2mat(PLOTs_med')';
        B = cell2mat(ERRORs_med')';
        for i = 1:9
            errorbar(A(:,i) , B(:,i) , 'lineWidth' , 3,'color' , colors(i,:));
            hold on
        end
        
        grid on
        set(gca , 'XLim' , [0 6] ,'XTick' , [1:5] , 'YLim', [2500 7000] ,'FontSize' , 20, 'Box' , 'off',...
            'GridAlpha' , 1)
        ylabel('msec')
        xlabel('Days')
        title('Execution Time in  Chunked Sequences')
        %         legend(fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' }));
        
        subplot(224)
        A = cell2mat(PLOTr_med')';
        B = cell2mat(ERRORr_med')';
        for i = 1:9
            errorbar(A(:,i) , B(:,i) , 'lineWidth' , 3,'color' , colors(i,:));
            hold on
        end
        grid on
        set(gca , 'XLim' , [0 6] , 'XTick' , [1:5] , 'YLim', [2500 7000] ,'FontSize' , 20 ,...
            'GridAlpha' , 1, 'Box' , 'off')
        ylabel('msec')
        xlabel('Days')
        title('Execution Time in  Random Sequences')
        legend(fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' }) , 'Box' , 'off');
        
        
        
        
        figure('color' , 'white')
        As = cell2mat(PLOTs_med')';
        Bs = cell2mat(ERRORs_med')';
        Ar = cell2mat(PLOTr_med')';
        Br = cell2mat(ERRORr_med')';
        horz = {[1] [2] [3] [4:9]};
        hlab = {'1' , '2' , '3' , '4 - 13'};
        for i = 1:length(horz)
            subplot(2,2,i)
            errorbar(mean(As(:,horz{i}) , 2) , mean(Bs(:,horz{i}) , 2) , 'lineWidth' , 3);
            hold on
            errorbar(mean(Ar(:,horz{i}) , 2) , mean(Br(:,horz{i}) , 2) , 'lineWidth' , 3);
            grid on
            set(gca , 'XLim' , [0 6] ,'XTick' , [1:5] , 'YLim', [2500 7000] ,'FontSize' , 20, 'Box' , 'off',...
                'GridAlpha' , 1)
            ylabel('msec')
            xlabel('Days')
            title(['Execution Time in horizon(s) ' , hlab{i}])
        end
        legend('Chunked' , 'Random' , 'Box' , 'off');
        
        %
        %         for h  = [1:7]
        %             out.chunked_exclude_h_p(h,1) = out.HorizonEffect_Chunkedseq_not(h).eff(2).p;
        %             disp(['Effect of Horizon on Chunked MTs Excluding Horizons ',num2str([1:h]),' p val = ' , num2str(out.HorizonEffect_Chunkedseq_not(h).eff(2).p)])
        %             out.random_exclude_h_p(h,1) = out.HorizonEffect_Randomseq_not(h).eff(2).p;
        %             disp(['Effect of Horizon on Random MTs Excluding Horizons ',num2str([1:h]),' p val = ' , num2str(out.HorizonEffect_Randomseq_not(h).eff(2).p)])
        %         end
        %         out.chunked_exclude_h_p(h+1,1) = out.HorizonEffect_Chunkedseq_allh.eff(2).p;
        %         out.random_exclude_h_p(h+1,1) = out.HorizonEffect_Randomseq_allh.eff(2).p;
    case 'IPI_Vs_Horizon'
        structNumb = input('Which structure? (1/2)');
        %         plotfcn = input('nanmean or nanmedian?' , 's');
        %% IPIs vs horizon
        % this is the output of the case: 'transitions_All' that is saved to disc
        if ~ calc
            load([baseDir , '/se2_TranProb.mat'] , 'All');
            newANA = All;
            newANA = getrow(newANA , ismember(newANA.Day , days{day}));
        else
            ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 , structNumb])  & ~Dall.isError & ismember(Dall.Day , days{day}));
            
            for tn = 1:length(ANA.TN)
                n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
                nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
                ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
            end
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = diff(ANA.ChnkArrang(tn,:));
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a-1) = 3;
                ANA.ChunkBndry(tn , end) = 3;
                ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                ANA.ChunkBndry(tn , 1:2) = [-1 -1];  % dont account for the first and last sseqeuce presses
                ANA.ChunkBndry(tn , end-1:end) = [-2 -2];% dont account for the first and last sseqeuce presses
                ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,13);
                ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,13);
                ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,13);
                ANA.IPI_prsnumb(tn , :) = [1 :13];
                ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,13);
            end
            
            newANA.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
            newANA.ChunkBndry = reshape(ANA.ChunkBndry , numel(ANA.IPI) , 1);
            newANA.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
            newANA.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
            newANA.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
            newANA.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
            newANA.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
            newANA.ChunkBndry(newANA.ChunkBndry>2) = 2;
        end
        IPItable  = tapply(newANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb' , 'ChunkBndry' , 'prsnumb'} , {'IPI' , 'nanmedian(x)'});
        
        lineplot([IPItable.ChunkBndry IPItable.Day], IPItable.IPI, 'subset' , ismember(IPItable.Horizon,[3:13]) & ismember(IPItable.seqNumb, structNumb) & ~ismember(IPItable.ChunkBndry, [-1 , -2]))
        
        h0 = figure;
        hC = 1;
        for h = [1:8 , 13]
            [IPI0_x(hC,:),IPI0_plot(hC,:) , IPI0_plot_error(hC,:)] = lineplot(IPItable.prsnumb ,IPItable.IPI , 'subset' , IPItable.Horizon == h & ismember(IPItable.seqNumb , 0),'plotfcn' , 'nanmean');
            hold on
            [IPI1_x(hC,:),IPI1_plot(hC,: ), IPI1_plot_error(hC,:)] = lineplot(IPItable.prsnumb ,IPItable.IPI  , 'subset' , IPItable.Horizon == h & ismember(IPItable.seqNumb , structNumb) ,'plotfcn' , 'nanmean');
            hC = hC+1;
        end
        close(h0)
        
        
        h0 = figure;
        [x1 , plot1 , error1] = lineplot(IPItable.Horizon , IPItable.IPI , 'subset' , IPItable.ChunkBndry == 1 ,'plotfcn' , 'nanmean');
        hold on
        [x2 , plot2 , error2] = lineplot(IPItable.Horizon , IPItable.IPI , 'subset' , IPItable.ChunkBndry == 2 ,'plotfcn' , 'nanmean');
        
        close(h0)
        ChunkedIPI  = pivottable(newANA.Horizon, newANA.ChunkBndry, newANA.IPI ,'nanmedian(x)' , 'subset' , ~ismember(newANA.ChunkBndry ,[-1 -2]));
        
        
        out.IPI_allh = anovaMixed(IPItable.IPI , IPItable.SN , 'within' , [IPItable.ChunkBndry , IPItable.Horizon] , {'within/between' , 'horizon'} , 'subset' ,...
            ismember(IPItable.Horizon,[1:13]) & ismember(IPItable.seqNumb, structNumb) & ~ismember(IPItable.ChunkBndry, [-1 , -2]));
        for h = [1:8]
            temp = anovaMixed(IPItable.IPI , IPItable.SN , 'within' , [IPItable.ChunkBndry , IPItable.Horizon] , {'within/between' , 'horizon'} , 'subset' ,...
                ~ismember(IPItable.Horizon,[1:h]) & ismember(IPItable.seqNumb, structNumb) & ~ismember(IPItable.ChunkBndry, [-1 , -2]));
            out.W_B_IPI_not_h(h,1) = temp.eff(2).p;
        end
        out.W_B_IPI_not_h(h+1,1) = out.IPI_allh.eff(2).p;
        
        disp(['Effect of Between/within chunk on Chunked IPIs including all Horizons p val = ' , num2str(out.IPI_allh.eff(2).p)])
        for h  = [1:8]
            disp(['Effect of Between/within chunk on Chunked IPIs Excluding Horizons ',num2str([1:h]),' p val = ' , num2str(out.W_B_IPI_not_h(h))])
        end
        
        figure('color' , 'white');
        subplot(2,2,1)
        imagesc(IPI1_plot, [100 600]);
        colorbar
        hold on
        title(['Structure ' , num2str(structNumb) ,' - IPI vs Horizons'])
        xlabel('Press Number')
        set(gca,'FontSize' , 16, 'XTick' , [1:13] ,'YTick' , [1:9] , 'YTickLabel' ,...
            fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' }));
        axis square
        
        
        subplot(2,2,2)
        imagesc(IPI0_plot , [100 600]);
        colorbar
        title(['Random - IPI vs Horizons'])
        xlabel('Press Number'  ,'FontSize' , 10)
        ax.XTick = [1:13];
        set(gca,'FontSize' , 16, 'XTick' , [1:13] ,'YTick' , [1:9], ...
            'YTickLabel', fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' }));
        axis square
        
        
        
        
        subplot(2,2,3)
        for horzz = 1:9
            errorbar(IPI1_x(horzz,:)',IPI1_plot(horzz,:) , IPI1_plot_error(horzz,:) , 'LineWidth' , 3 , 'color' , colors(horzz  , :))
            hold on
        end
        
        title(['Structure ' , num2str(structNumb) , ' - Day ' , num2str(days{day})])
        grid on
        xlabel('Presses')
        ylabel('msec')
        set(gca,'FontSize' , 16, 'XTick' , [1:13] , 'YLim' , [100 600] , 'Box' , 'off' , 'GridAlpha' , 1);
        
        subplot(2,2,4)
        for horzz = 1:9
            errorbar(IPI0_x(horzz,:)',IPI0_plot(horzz,:) , IPI0_plot_error(horzz,:) , 'LineWidth' , 3 , 'color' , colors(horzz  , :))
            hold on
        end
        legend(fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' }))
        title(['Random - Day ' , num2str(days{day})])
        grid on
        xlabel('Presses')
        ylabel('msec')
        set(gca,'FontSize' , 16, 'XTick' , [1:13], 'YLim' , [100 600] , 'Box' , 'off' , 'GridAlpha' , 1);
        %%
        figure('color' , 'white');
        subplot(1,2,1)
        imagesc(ChunkedIPI);
        colorbar
        hold on
        set(gca,'FontSize' , 16 , 'XTick' , [1 :2] , 'XTickLabel' , {'First' , 'Middle' } , ...
            'YTickLabel' , fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' }))
        title('Median IPI vs Horizons in the Chunked sequences')
        
        
        subplot (1,2,2)
        h1 = plotshade(x1',plot1,error1,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':')
        hold on
        h2 = plotshade(x2',plot2,error2,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':')
        %errorbar(xcoord1,PLOT1,ERROR1,'LineWidth' , 3)
        %hold on
        %errorbar(xcoord2,PLOT2,ERROR2,'LineWidth' , 3)
        %errorbar(xcoord3,PLOT3,ERROR3,'LineWidth' , 3)
        hold on
        ax = gca;
        title('Chunked sequence IPIs vs horizon')
        xlabel('Horizon')
        ylabel('msec')
        legend([h1 h2] , {'First Chunk Press' , 'Middle Chunk Press' })
        grid on
        ax.XLim = [0 9];
        ax.FontSize = 16;
    case 'RT_Vs_Horizon'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & ismember(Dall.Day , days{day}));
        RT = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb','BN'} , {'pressTime0' , 'nanmedian(x)'});
        RT.pressTime0 = RT.pressTime0 - 1500;
        
        RT.seqNumb(RT.seqNumb == 2) = 1;
        Hcounter = 1;
        
        
        
        for h = [1:8 , 13]
            temp = anovaMixed(RT.pressTime0 , RT.SN,'within',RT.seqNumb,{'randon/chunked'},'intercept',1,'subset' , ismember(RT.Horizon , [h]) & ismember(RT.seqNumb , [0 1]));
            out.R_C_RT_not_h(h,1) = temp.eff(2).p;
        end
        temp = anovaMixed(RT.pressTime0 , RT.SN,'within',RT.seqNumb,{'randon/chunked'},'intercept',1,'subset', ismember(RT.seqNumb , [0 1]));
        out.R_C_RT_not_h(h+1,1) = temp.eff(2).p;
        
        disp(['Effect of Random/Chunked on RT including all Horizons p val = ' , num2str(out.R_C_RT_not_h(14,1))])
        for h  = [1:8 , 13]
            disp(['Effect of Random/Chunked on RT in Horizon ',num2str([h]),' p val = ' , num2str(out.R_C_RT_not_h(h))])
        end
        
        
        
        h=figure('color' , 'white');
        [xcoord1,ePLOT1,ERROR1] = lineplot(RT.Horizon , RT.pressTime0 ,'subset', ismember(RT.seqNumb , [1]));
        hold on
        [xcoord0,ePLOT0,ERROR0] = lineplot(RT.Horizon , RT.pressTime0 ,'subset', ismember(RT.seqNumb , [0]));
        close(h)
        figure('color' , 'white')
        h1 = plotshade(xcoord1',ePLOT1,ERROR1,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
        hold on
        h2 = plotshade(xcoord0',ePLOT0,ERROR0,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':');
        legend([h1 , h2] , {'Chunked' , 'Random'})
        grid on
        xlabel('Horizon')
        ylabel('msec')
        title(['Reaction time vs. horizon Days = ' , num2str(days{day})])
        set(gca ,'XTick' ,[1:8 , 13],'FontSize' , 20);
    case 'Points'
        for tn = 1:length(Dall.TN)
            Dall.MT(tn , 1) = Dall.AllPressTimes(tn , Dall.seqlength(tn)) - Dall.AllPressTimes(tn , 1);
        end
        pionts_rate = [];
        MT = [];
        MTall = [];
        for h = [1:8 , 13]
            for sub = subjnum
                ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [1:2]) & Dall.isgood & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
                pionts_rate = [pionts_rate ; [100 * sum(ANA.points == 3)/length(ANA.TN) ,...
                    100 * sum(ANA.points == 1)/length(ANA.TN) ,...
                    100 * sum(ANA.points == 0 & ~ANA.isError)/length(ANA.TN),...
                    100 * sum(ANA.isError)/length(ANA.TN), h , sub ,1]];
                MTall = [MTall ; [ANA.MT ANA.Horizon]];
                MT = [MT ; [nanmean(ANA.MT(ANA.points == 3)) ,...
                    nanmean(ANA.MT(ANA.points == 1)) ,...
                    nanmean(ANA.MT(ANA.points == 0 & ~ANA.isError)),...
                    nanmean(ANA.MT(ANA.isError)), h , sub ,1]];
                ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [0]) & Dall.isgood & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
                pionts_rate = [pionts_rate ; [100 * sum(ANA.points == 3)/length(ANA.TN) ,...
                    100 * sum(ANA.points == 1)/length(ANA.TN) ,...
                    100 * sum(ANA.points == 0)/length(ANA.TN) ,...
                    100 * sum(ANA.isError)/length(ANA.TN) ,h , sub ,0]];
                MT = [MT ; [nanmean(ANA.MT(ANA.points == 3)) ,...
                    nanmean(ANA.MT(ANA.points == 1)) ,...
                    nanmean(ANA.MT(ANA.points == 0 & ~ANA.isError)) ,...
                    nanmean(ANA.MT(~ANA.isError)) , h , sub ,0 ]];
                MTall = [MTall ; [ANA.MT ANA.Horizon]];
            end
        end
        MTall = MTall(MTall(:,1) <= 9000 , :);
        
        P3_temp = pivottable(pionts_rate(:,5) ,pionts_rate(:,7), pionts_rate(:,1) , 'nanmedian');
        P2_temp = pivottable(pionts_rate(:,5) ,pionts_rate(:,7), pionts_rate(:,2) , 'nanmedian');
        P1_temp = pivottable(pionts_rate(:,5) ,pionts_rate(:,7), pionts_rate(:,3) , 'nanmedian');
        E_temp  = pivottable(pionts_rate(:,5) ,pionts_rate(:,7), pionts_rate(:,4) , 'nanmedian');
        
        MT3_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,1) , 'nanmedian');
        MT2_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,2) , 'nanmedian');
        MT1_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,3) , 'nanmedian');
        MTe_temp  = pivottable(MT(:,5) ,MT(:,7), MT(:,4) , 'nanmedian');
        
        out.Err = anovan(pionts_rate(:,4) , pionts_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        out.Po3 = anovan(pionts_rate(:,1) , pionts_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        out.Po2 = anovan(pionts_rate(:,2) , pionts_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        out.Po1 = anovan(pionts_rate(:,3) , pionts_rate(:,[5,7]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunked'});
        
        
        figure('color' , 'white')
        subplot(1,2,1)
        imagesc(E_temp)
        colorbar
        hold on
        ax = gca;
        %         ax.FontSize = 20;
        title(['Error rate vs Horizons -interaction = ' ,num2str(out.Err(3))] )
        ax.XTick = [1 :3];
        ax.XTickLabel = {'Random' , 'Chunked'};
        ax.YTickLabel  = fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' });
        
        subplot(1,2,2)
        imagesc(P3_temp)
        colorbar
        hold on
        ax = gca;
        %         ax.FontSize = 20;
        title(['3-pionts rate vs Horizons -interaction = ' ,num2str(out.Po3(3))] )
        ax.XTick = [1 :3];
        ax.XTickLabel = {'Random' , 'Chunked'};
        ax.YTickLabel  = fliplr({'H = 13' 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' });
        
        figure('color' , 'white')
        subplot(4,1,1)
        barplot(pionts_rate(:,[5,7]) , pionts_rate(:,1));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('% of all the trials')
        title('P = 3')
        
        
        subplot(4,1,2)
        barplot(pionts_rate(:,[5,7]) , pionts_rate(:,2));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('% of all the trials')
        title('P = 1')
        
        subplot(4,1,3)
        barplot(pionts_rate(:,[5,7]) , pionts_rate(:,3));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('% of all the trials')
        title('P = 0')
        
        subplot(4,1,4)
        barplot(pionts_rate(:,[5,7]) , pionts_rate(:,4));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('% of all the trials')
        title('Error')
        
        
        figure('color' , 'white')
        subplot(4,1,1)
        barplot(MT(:,[5,7]) , MT(:,1));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('msec')
        title('MT in P = 3')
        
        subplot(4,1,2)
        barplot(MT(:,[5,7]) , MT(:,2));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('msec')
        title('MT in P = 1')
        
        subplot(4,1,3)
        barplot(MT(:,[5,7]) , MT(:,3));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('msec')
        title('MT in P = 0')
        
        subplot(4,1,4)
        barplot(MT(:,[5,7]) , MT(:,4));
        xlabel('Horizon/0(Random) OR 1(Chunked)')
        ylabel('msec')
        title('MT in Error')
    case 'MT_asymptote'
        Hex = input('What horizons to exclude? (0 = include all)');
        ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError &~ismember(Dall.Horizon , Hex));
        ANA.seqNumb(ANA.seqNumb == 2) = 1;
        MT  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmedian(x)'});
        MT.MT_pred = zeros(size(MT.MT));
        MT.b1 = zeros(size(MT.MT));
        MT.b2 = zeros(size(MT.MT));
        MT.b3 = zeros(size(MT.MT));
        
        for d = 1:5
            for subjnum = 1:length(subj_name)-1
                [d subjnum 1]
                id = ismember(MT.SN , subjnum) & ismember(MT.Day , d) & ismember(MT.seqNumb , [1]);
                MTsn = getrow(MT , id);
                exp_model1 = @(b,x) b(1) + (b(2) - b(1))*exp(-(x-1)/b(3)); % Model Function
                %                 exp_model1 = @(b,x) b(1) + b(2)*exp(b(3)*x); % Model Function
                x = [1:length(unique(MTsn.Horizon))];
                yx = MTsn.MT';                                    % this would be a typical MT vs Horizon vector: [5422 3548 2704 2581 2446 2592 2418 2528 2500]
                OLS = @(b) sum((exp_model1(b,x) - yx).^2);                % Ordinary Least Squares cost function
                opts = optimset('MaxIter', 300,'TolFun',1e-5);
                [B1 Fval] = fminsearch(OLS,[3500 7500  1], opts);        % Use ?fminsearch? to minimise the ?OLS? function
                MT.MT_pred(id) = exp_model1(B1,x);
                MT.b1(id) = B1(1);
                MT.b2(id) = B1(2);
                MT.b3(id) = B1(3);
                
                
                [d subjnum 2]
                id = ismember(MT.SN , subjnum) & ismember(MT.Day , d) & ismember(MT.seqNumb , [0]);
                MTsn = getrow(MT , id);
                exp_model0 = @(b,x) b(1) + (b(2) - b(1))*exp(-(x-1)/b(3)); % Model Function
                yx = MTsn.MT';                                    % this would be a typical MT vs Horizon vector: [5422 3548 2704 2581 2446 2592 2418 2528 2500]
                OLS = @(b) sum((exp_model0(b,x) - yx).^2);                % Ordinary Least Squares cost function
                opts = optimset('MaxIter', 300,'TolFun',1e-5);
                B0 = fminsearch(OLS,[3500 7500 1], opts);        % Use ?fminsearch? to minimise the ?OLS? function
                MT.MT_pred(id) = exp_model0(B0,x);
                MT.b1(id) = B0(1);
                MT.b2(id) = B0(2);
                MT.b3(id) = B0(3);
            end
            [coo1(:,d),plot1(:,d),err1(:,d)] = lineplot([MT.Horizon] , MT.MT , 'subset' , ismember(MT.seqNumb , [1])  & ismember(MT.Day , d));
            [coo0(:,d),plot0(:,d),err0(:,d)] = lineplot([MT.Horizon] , MT.MT , 'subset' , ismember(MT.seqNumb , [0]) & ismember(MT.Day , d));
            
            [coo1_pred(:,d),plot1_pred(:,d),err1_pred(:,d)] = lineplot([MT.Horizon] , MT.MT_pred , 'subset' , ismember(MT.seqNumb , [1]) & ismember(MT.Day , d));
            [coo0_pred(:,d),plot0_pred(:,d),err0_pred(:,d)] = lineplot([MT.Horizon] , MT.MT_pred , 'subset' , ismember(MT.seqNumb , [0]) & ismember(MT.Day , d));
        end
        %% if you wanted to d othe parameter estimation in NonLinearModel:
        %         % New Model
        %         exp_model1 = @(b,x) b(1) + (b(2) - b(1))*exp((x-1)/b(3)); % Model Function
        %         opts1 = statset('Display','final','TolFun',1e-5, 'MaxIter', 100);
        %         % regressing out the horizon effect, i.e the exponential decrease in the MT
        %         % Horizon number is the input and MT is the output
        %         nlmf1 = NonLinearModel.fit([1:8 13], MTsn.MT, exp_model1, initial_values, 'Options', opts1);
        %         MT.MT_pred(id) = nlmf1.Fitted;
        %         MT.b1(id) = nlmf1.Coefficients.Estimate(1)*ones(sum(id),1);
        %         MT.b2(id) = nlmf1.Coefficients.Estimate(2)*ones(sum(id),1);
        %         MT.b3(id) = nlmf1.Coefficients.Estimate(3)*ones(sum(id),1);
        %
        
        
        %%
        h0 = figure;
        [coo1_b1,plot1_b1,err1_b1] = lineplot([MT.Day] , MT.b1 , 'subset' , ismember(MT.seqNumb , [1]) );
        [coo0_b1,plot0_b1,err0_b1] = lineplot([MT.Day] , MT.b1 , 'subset' , ismember(MT.seqNumb , [0]) );
        
        [coo1_b2,plot1_b2,err1_b2] = lineplot([MT.Day] , MT.b2 , 'subset' , ismember(MT.seqNumb , [1]) );
        [coo0_b2,plot0_b2,err0_b2] = lineplot([MT.Day] , MT.b2 , 'subset' , ismember(MT.seqNumb , [0]) );
        
        [coo1_b3,plot1_b3,err1_b3] = lineplot([MT.Day] , MT.b3 , 'subset' , ismember(MT.seqNumb , [1]) );
        [coo0_b3,plot0_b3,err0_b3] = lineplot([MT.Day] , MT.b3 , 'subset' , ismember(MT.seqNumb , [0]) );
        
        [coo1_invb3,plot1_invb3,err1_invb3] = lineplot([MT.Day] , (MT.b3).^-1 , 'subset' , ismember(MT.seqNumb , [1]) );
        [coo0_invb3,plot0_invb3,err0_invb3] = lineplot([MT.Day] , (MT.b3).^-1 , 'subset' , ismember(MT.seqNumb , [0]) );
        close(h0)
        for d = 1:5
            id1  = MT.Day == d & MT.seqNumb ==1;
            id2  = MT.Day == d & MT.seqNumb ==0;
            [~,out.b1_significance(d)] = ttest2(MT.b1(id1) , MT.b1(id2));
            [~,out.b2_significance(d)] = ttest2(MT.b2(id1) , MT.b2(id2));
            [~,out.b3_significance(d)] = ttest2(MT.b3(id1) , MT.b3(id2));
        end
        
        figure('color' , 'white')
        for d = 1:5
            subplot(2,3,d)
            errorbar(coo1(:,d)',plot1(:,d)',err1(:,d)' , 'LineWidth' , 3);
            hold on
            errorbar(coo0(:,d)',plot0(:,d)',err0(:,d)' , 'LineWidth' , 3);
            set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1)
            hold on
            xlabel('Horizon')
            ylabel('msec')
            title(['Chunked MT vs. Random on Day ' , num2str(d)])
            grid on
        end
        legend({'Chunked' , 'Random'}, 'Box' , 'off')
        
        figure('color' , 'white')
        for d = 1:5
            subplot(2,3,d)
            errorbar(coo1_pred(:,d)',plot1_pred(:,d)',err1_pred(:,d)' , 'LineWidth' , 3);
            hold on
            errorbar(coo0_pred(:,d)',plot0_pred(:,d)',err0_pred(:,d)' , 'LineWidth' , 3);
            set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1)
            hold on
            xlabel('Horizon')
            ylabel('msec')
            title(['fitted Chunked vs. fitted Random- Day ' , num2str(d)])
            grid on
        end
        legend({'Chunked' , 'Random'}, 'Box' , 'off')
        
        figure('color' , 'white')
        for d = 1:5
            subplot(2,3,d)
            errorbar(coo1(:,d)',plot1(:,d)',err1(:,d)' , 'LineWidth' , 3,'color','r');
            hold on
            errorbar(coo1_pred(:,d)',plot1_pred(:,d)',err1_pred(:,d)' , 'LineWidth' , 3 ,'color','c','LineStyle',':');
            set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1 , 'Box' , 'off')
            hold on
            xlabel('Horizon')
            ylabel('msec')
            title(['Chunked vs. fitted - Day ' , num2str(d)])
            grid on
        end
        legend({'Chunked' , 'Fitted Chunked'})
        
        figure('color' , 'white')
        for d = 1:5
            subplot(2,3,d)
            errorbar(coo0(:,d)',plot0(:,d)',err0(:,d)' , 'LineWidth' , 3,'color','r');
            hold on
            errorbar(coo0_pred(:,d)',plot0_pred(:,d)',err0_pred(:,d)' , 'LineWidth' ,3,'color','c','LineStyle',':');
            set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
            hold on
            xlabel('Horizon')
            ylabel('msec')
            title(['Random vs. fitted Random- Day ' , num2str(d)])
            grid on
        end
        legend({'Random' , 'Fitted Random'}, 'Box' , 'off')
        
        
        figure('color' , 'white')
        subplot(221)
        for d = 1:5
            errorbar(coo1(:,d)',plot1(:,d)',err1(:,d)' , 'LineWidth' , 3);
            hold on
        end
        set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
        hold on
        xlabel('Horizon')
        ylabel('msec')
        title('Chunked MT vs. Horizon')
        grid on
        
        subplot(222)
        for d = 1:5
            errorbar(coo0(:,d)',plot0(:,d)',err0(:,d)' , 'LineWidth' , 3);
            hold on
        end
        set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
        xlabel('Horizon')
        ylabel('msec')
        title('Random MT vs. Horizon')
        legend({'Day1' , 'Day2' ,'Day3','Day4','Day5'}, 'Box' , 'off')
        grid on
        
        subplot(223)
        for d = 1:5
            errorbar(coo1_pred(:,d)',plot1_pred(:,d)',err1_pred(:,d)' , 'LineWidth' , 3);
            hold on
        end
        set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
        hold on
        xlabel('Horizon')
        ylabel('msec')
        title('Fitted Chunked MT vs. Horizon')
        grid on
        
        subplot(224)
        for d = 1:5
            errorbar(coo0_pred(:,d)',plot0_pred(:,d)',err0_pred(:,d)' , 'LineWidth' , 3);
            hold on
        end
        set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20 , 'GridAlpha' , 1, 'Box' , 'off')
        xlabel('Horizon')
        title('Fitted Random MT vs. Horizon')
        grid on
        
        figure('color' , 'white')
        subplot(141)
        errorbar(coo1_b1,plot1_b1,err1_b1 , 'LineWidth' , 3);
        hold on
        errorbar(coo0_b1,plot0_b1,err0_b1 , 'LineWidth' , 3);
        set(gca, 'XTick' , [1:5],'FontSize' , 20, 'Box' , 'off','GridAlpha' , 1)
        hold on
        xlabel('Day')
        title('b1 Coefficient')
        grid on
        
        subplot(142)
        errorbar(coo1_b2,plot1_b2,err1_b2 , 'LineWidth' , 3);
        hold on
        errorbar(coo0_b2,plot0_b2,err0_b2 , 'LineWidth' , 3);
        set(gca, 'XTick' , [1:5],'FontSize' , 20, 'Box' , 'off','GridAlpha' , 1)
        hold on
        xlabel('Day')
        title('b2 Coefficient')
        grid on
        
        subplot(143)
        errorbar(coo1_b3,plot1_b3,err1_b3 , 'LineWidth' , 3);
        hold on
        errorbar(coo0_b3,plot0_b3,err0_b3 , 'LineWidth' , 3);
        set(gca, 'XTick' , [1:5],'FontSize' , 20, 'Box' , 'off','GridAlpha' , 1)
        hold on
        xlabel('Day')
        ylabel('msec')
        title('b3 Coefficient')
        grid on
        legend({'Chunked' , 'Random'}, 'Box' , 'off')
        
        subplot(144)
        errorbar(coo1_invb3,plot1_invb3,err1_invb3 , 'LineWidth' , 3);
        hold on
        errorbar(coo0_invb3,plot0_invb3,err0_invb3 , 'LineWidth' , 3);
        set(gca, 'XTick' , [1:5],'FontSize' , 20, 'Box' , 'off','GridAlpha' , 1)
        hold on
        xlabel('Day')
        ylabel('msec')
        title('1/b3 (Decay constant)')
        grid on
        legend({'Chunked' , 'Random'}, 'Box' , 'off')
        
        x = [1:10];
        b1 = 6;
        b2 = 10;
        figure('color' , 'white')
        %         subplot(121)
        for b3 = .2:.1:.7
            plot(x , b1 + (b2 - b1)*exp(-(x-1)/b3) ,'LineWidth' , 3)
            hold on
        end
        legend({'b3  = 0.2' , 'b3  = 0.3' ,'b3  = 0.4' ,'b3  = 0.5' ,'b3  = 0.6' , 'b3 = 0.7'}, 'Box' , 'off')
        title(['b1 + (b2 - b1)*exp(-(x-1)/b3)         for          b1 = ' ,num2str(b1)  , ',     b2 = ' ,num2str(b2)])
        set(gca , 'FontSize' , 20, 'Box' , 'off', 'GridAlpha' , 1, 'Box' , 'off')
        grid on
        x = [1:35];
        b1 = 6000;
        b3 = 3;
        %         figure('color' , 'white')
        subplot(122)
        for b2 = 0:20000:10^5
            plot(x , b1 + (b2 - b1)*exp(-(x-1)/b3) ,'LineWidth' , 3)
            hold on
        end
        legend({'b2  = 0' , 'b2  = 20000 ' ,'b2  = 40000' ,'b2  = 60000' ,'b2  = 80000' ,'b2  = 100000'}, 'Box' , 'off')
        title(['b1 = ' ,num2str(b1)  , ',     b3 = ' ,num2str(b3)])
        set(gca , 'FontSize' , 20 , 'Box' , 'off', 'GridAlpha' , 1)
        grid on
    case 'test_MT_asymptote'
        Hex = input('What horizons to exclude? (0 = include all)');
        ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError &~ismember(Dall.Horizon , Hex) & ismember(Dall.Day , [3:5]));
        ANA.seqNumb(ANA.seqNumb == 2) = 1;
        MT  = tapply(ANA , {'Horizon' , 'seqNumb'} , {'MT' , 'nanmedian(x)'});
        MT.MT_pred = zeros(size(MT.MT));
        MT.b1 = zeros(size(MT.MT));
        MT.b2 = zeros(size(MT.MT));
        MT.b3 = zeros(size(MT.MT));
        subjnum = 1 : length(subj_name) - 1;
        init_val = [3500 7500];   % [b1 b2]
        purt = 50*rand(1,100);
        itermax = [20 50 80 110 200 300];
        
        for i = 1:length(purt)
            d = 1
            id = ismember(MT.seqNumb , [1]);
            MTsn = getrow(MT , id);
            exp_model1 = @(b,x) b(1) + (b(2) - b(1))*exp(-(x-1)/b(3)); % Model Function
            %                 exp_model1 = @(b,x) b(1) + b(2)*exp(b(3)*x); % Model Function
            x = [1:length(unique(MTsn.Horizon))];
            yx = MTsn.MT';                                    % this would be a typical MT vs Horizon vector: [5422 3548 2704 2581 2446 2592 2418 2528 2500]
            OLS = @(b) sum((exp_model1(b,x) - yx).^2);                % Ordinary Least Squares cost function
            for j = 1:length(itermax)
                opts = optimset( 'MaxIter', itermax(j),'TolFun',1e-5);
                [B1{j}(:,i) Fval1{j}(i)] = fminsearch(OLS,[init_val+purt(i)*init_val 1], opts);        % Use ?fminsearch? to minimise the ?OLS? function
                MT_pred1{j}(i,:) = exp_model1(B1{j}(:,i),x);
            end
            
            id = ismember(MT.seqNumb , [0]);
            MTsn = getrow(MT , id);
            exp_model0 = @(b,x) b(1) + (b(2) - b(1))*exp(-(x-1)/b(3)); % Model Function
            yx = MTsn.MT';                                    % this would be a typical MT vs Horizon vector: [5422 3548 2704 2581 2446 2592 2418 2528 2500]
            OLS = @(b) sum((exp_model0(b,x) - yx).^2);                % Ordinary Least Squares cost function
            for j = 1:length(itermax)
                opts = optimset('MaxIter', itermax(j),'TolFun',1e-5);
                [B0{j}(:,i) fval0{j}(i)]= fminsearch(OLS,[init_val+purt(i)*init_val 1], opts);        % Use ?fminsearch? to minimise the ?OLS? function
                MT_pred0{j}(i,:) = exp_model0(B0{j}(:,i),x);
            end
        end
        figure('color' , 'white')
        for j = 1:length(itermax)
            subplot(2,3,j)
            plot(fval0{j} , purt , '*' , 'MarkerSize' , 5)
            set(gca, 'FontSize' , 20)
            xlabel('MSE')
            ylabel('percentage')
            
            
            hold on
            plot(Fval1{j} , purt , '*' , 'MarkerSize' , 5)
            set(gca, 'FontSize' , 20)
            xlabel('MSE')
            ylabel('percentage')
            grid on
            legend({'Random' , 'Chunked'})
            title(['SSE with iteration of ' ,num2str(itermax(j))])
        end
    case 'eyepress_pos_traces'
        clear len EyePatterne pressPattern eye press chunk chunkPattern estchunkPattern
        if isequal(subjnum , [1:length(subj_name)-1])
            subjnum = length(subj_name);
        end
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:2
                out.EyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            if ~isempty(rep)
                switch rep
                    case 1
                        load([baseDir , '/se2_tnorm_r1.mat'])
                    case 2
                        load([baseDir , '/se2_tnorm_r2.mat'])
                    otherwise
                        load([baseDir , '/se2_tnorm.mat'])
                end
            else
                load([baseDir , '/se2_tnorm.mat'])
            end
        end
        
        for h = [1:8 , 13]
            out(h).eye = N(h).norm(day,subjnum).eye;
            out(h).press= N(h).norm(day,subjnum).press;
            out(h).eyeveloc= N(h).norm(day,subjnum).eyeveloc;
            out(h).prsveloc= N(h).norm(day,subjnum).prsveloc;
            
            out(h).EyePattern  = N(h).norm(day,subjnum).eyePattern;
            out(h).pressPattern = N(h).norm(day,subjnum).pressPattern;
            out(h).eyevelocPattern = N(h).norm(day,subjnum).eyevelocPattern;
            out(h).prsvelocPattern = N(h).norm(day,subjnum).prsvelocPattern;
        end
        %% eye vs press
        figure('color' , [1 1 1])
        h_counter = 1;
        for h = 1:8
            subplot(3,4,h_counter)
            for s = 2:3
                plot(out(h).pressPattern(s,:) , out(h).EyePattern(s,:) , 'LineWidth' , 3);
                hold on
            end
            plot(out(h).pressPattern(1,:) , out(h).EyePattern(1,:)  , 'LineWidth' , 3 , 'color' , 'black');
            grid on
            xlabel('Finger Press' , 'FontSize' , 10)
            ylabel('Eye Position' , 'FontSize' , 10)
            title(['H =  ' , num2str(h)] )
            %             legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
            hold on
            line([1,14] , [1 14] ,'LineStyle' , ':','color' , 'r' , 'LineWidth' , 3  )
            ax = gca;
            ax.YTick = [1:14];
            ax.YLim  = [1 14];
            ax.XTick = [1:14];
            ax.XLim  = [1 14];
            %             ax.FontSize = 20;
            ax.Box = 'off';
            axis square
            h_counter = h_counter +1;
        end
        
        %         figure('color' , [1 1 1])
        subplot(3,4,[9:12])
        for h = 13
            for s = 2:3
                plot(out(h).pressPattern(s,:) , out(h).EyePattern(s,:) , 'LineWidth' , 3);
                hold on
            end
            plot(out(h).pressPattern(1,:) , out(h).EyePattern(1,:)  , 'LineWidth' , 3 , 'color' , 'black');
            grid on
            xlabel('Finger Press' , 'FontSize' , 10)
            ylabel('Eye Position' , 'FontSize' , 10)
            title(['H =  ' , num2str(h)] )
            %             legend({ 'Structure 1','Structure 2','Structure 3','Structure 4','Structure 5','Structure 6' 'Random'} , 'FontSize' , 20)
            hold on
            line([1,14] , [1 14] ,'LineStyle' , ':','color' , 'r' , 'LineWidth' , 3  )
            ax = gca;
            ax.YTick = [1:14];
            ax.YLim = [1 14];
            ax.XTick = [1:14];
            ax.XLim = [1 14];
            %             ax.FontSize = 20;
            ax.Box = 'off';
            axis square
            h_counter = h_counter +1;
        end
        
        legend({ 'Structure 1','Structure 2', 'Random'} , 'FontSize' , 20)
        %% eye vs press
        %         figure('color' , [1 1 1])
        %         h_counter = 1;
        %         for h = 1:8
        %             subplot(2,4,h_counter)
        %             plot(out(h).EyePattern(2:end , :)' , 'LineWidth' , 3)
        %             hold on
        %             plot(out(h).EyePattern(1 , :)' , 'LineWidth' , 5 , 'color' , 'black')
        %             grid on
        %             xlabel('Normalized time' , 'FontSize' , 20)
        %             ylabel('Eye Position' , 'FontSize' , 20)
        %             title(['H =  ' , num2str(h)] )
        %             ax = gca;
        %             ax.YTick = [1:14];
        %             ax.YLim = [1 max(max(out(h).EyePattern))];
        %             ax.FontSize = 20;
        %             ax.Box = 'off';
        %             h_counter = h_counter+1;
        %         end
        
        %% last 10% of eye
        LTP = 0;
        if LTP
            figure('color' , [1 1 1])
            h_counter = 1;
            for h = 1:8
                subplot(3,4,h_counter)
                plot(out(h).EyePattern(2:end , 901:end)' , 'LineWidth' , 3)
                hold on
                plot(out(h).EyePattern(1 , 901:end)' , 'LineWidth' , 3 , 'color' , 'black')
                grid on
                xlabel('Normalized time' , 'FontSize' , 10)
                ylabel('Eye Position' , 'FontSize' , 10)
                title(['last 10% H =  ' , num2str(h)] )
                ax = gca;
                ax.FontSize = 10;
                ax.Box = 'off';
                ax.XLim = [1 100];
                %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
                ax.YLim = [12 15];
                h_counter = h_counter+1;
            end
            
            %         figure('color' , [1 1 1])
            subplot(3,4,[9:12])
            for h = 13
                plot(out(h).EyePattern(2:end , 901:end)' , 'LineWidth' , 3)
                hold on
                plot(out(h).EyePattern(1 , 901:end)' , 'LineWidth' , 3 , 'color' , 'black')
                grid on
                xlabel('Normalized time' , 'FontSize' , 10)
                ylabel('Eye Position' , 'FontSize' , 10)
                title(['last 10% H =  ' , num2str(h)] )
                ax = gca;
                ax.FontSize = 10;
                ax.Box = 'off';
                ax.XLim = [1 100];
                %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
                ax.YLim = [12 15];
            end
            legend({ 'Structure 1','Structure 2','Random'} , 'FontSize' , 20)
        end
        %% first 10% of eye
        FTP = 0
        if FTP
            figure('color' , [1 1 1])
            h_counter = 1;
            for h = 1:8
                subplot(3,4,h_counter)
                plot(out(h).EyePattern(2:end , 1:100)' , 'LineWidth' , 3)
                hold on
                plot(out(h).EyePattern(1 , 1:100)' , 'LineWidth' , 3 , 'color' , 'black')
                grid on
                xlabel('Normalized time' , 'FontSize' , 10)
                ylabel('Eye Position' , 'FontSize' , 10)
                title(['First 10% H =  ' , num2str(h)] )
                ax = gca;
                ax.FontSize = 10;
                ax.Box = 'off';
                ax.XLim = [1 100];
                %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
                ax.YLim = [1 4];
                h_counter = h_counter+1;
            end
            
            %         figure('color' , [1 1 1])
            subplot(3,4,[9:12])
            for h = 13
                plot(out(h).EyePattern(2:end , 1:100)' , 'LineWidth' , 3)
                hold on
                plot(out(h).EyePattern(1 , 1:100)' , 'LineWidth' , 3 , 'color' , 'black')
                grid on
                xlabel('Normalized time' , 'FontSize' , 10)
                ylabel('Eye Position' , 'FontSize' , 10)
                title(['First 10% H =  ' , num2str(h)] )
                ax = gca;
                ax.FontSize = 10;
                ax.Box = 'off';
                ax.XLim = [1 100];
                %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
                ax.YLim = [1 4];
            end
            legend({ 'Structure 1 ()','Structure 2', 'Random'} , 'FontSize' , 20)
        end
        out = [];
    case 'eyepress_vel_traces'
        if isequal(subjnum , [1:length(subj_name)-1])
            subjnum = length(subj_name);
        end
        %         calc = 0;
        if calc
            N = se1_timeNorm('single_sub' ,subjnum,day);
            out.eye = N.norm.eye;
            out.press= N.norm.press;
            out.eyeveloc= N.norm.eyeveloc;
            out.prsveloc= N.norm.prsveloc;
            for s = 0:6
                out.EyePattern(s+1, :)  = nanmean(out.eye{s+1});
                out.pressPattern(s+1,:) = nanmean(out.press{s+1});
                out.eyevelocPattern(s+1,:) = nanmean(out.eyeveloc{s+1});
                out.prsvelocPattern(s+1,:) = nanmean(out.prsveloc{s+1});
            end
        else
            if ~isempty(rep)
                switch rep
                    case 1
                        load([baseDir , '/se2_tnorm_r1.mat'])
                    case 2
                        load([baseDir , '/se2_tnorm_r2.mat'])
                    otherwise
                        load([baseDir , '/se2_tnorm.mat'])
                end
            else
                load([baseDir , '/se2_tnorm.mat'])
            end
        end
        
        for h = [1:8 , 13]
            out(h).eye = N(h).norm(day,subjnum).eye;
            out(h).press= N(h).norm(day,subjnum).press;
            out(h).eyeveloc= N(h).norm(day,subjnum).eyeveloc;
            out(h).prsveloc= N(h).norm(day,subjnum).prsveloc;
            
            out(h).EyePattern  = N(h).norm(day,subjnum).eyePattern;
            out(h).pressPattern = N(h).norm(day,subjnum).pressPattern;
            out(h).eyevelocPattern = N(h).norm(day,subjnum).eyevelocPattern;
            out(h).prsvelocPattern = N(h).norm(day,subjnum).prsvelocPattern;
        end
        
        
        
        
        %% Eye Vel
        figure('color' , [1 1 1])
        h_counter = 1;
        for h = 1:8
            subplot(3,4,h_counter)
            plot(out(h).eyevelocPattern(2:end,:)' , 'LineWidth' , 3)
            hold on
            plot(out(h).eyevelocPattern(1,:)' , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Normalized time' , 'FontSize' , 10)
            ylabel('deg/sec' , 'FontSize' , 10)
            title(['Eye Velocity H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            ax.YLim = [-3 5];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
            h_counter = h_counter+1;
        end
        
        subplot(3,4,[9:12])
        for h = 13
            plot(out(h).eyevelocPattern(2:end , :)' , 'LineWidth' , 3)
            hold on
            plot(out(h).eyevelocPattern(1 , :)' , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Normalized time' , 'FontSize' , 10)
            ylabel('deg/sec' , 'FontSize' , 10)
            title(['Eye Velocity H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            ax.YLim = [-3 5];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
        end
        legend({ 'Structure 1','Structure 2','Random'} , 'FontSize' , 20)
        
        %% Press VEl
        
        figure('color' , [1 1 1])
        h_counter = 1;
        for h = 1:8
            subplot(3,4,h_counter)
            plot(out(h).prsvelocPattern(2:end,:)' , 'LineWidth' , 3)
            hold on
            plot(out(h).prsvelocPattern(1,:)' , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Normalized time' , 'FontSize' , 10)
            ylabel('press/sec' , 'FontSize' , 10)
            title(['Press Velocity, H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            ax.YLim = [0 8];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
            h_counter = h_counter+1;
        end
        
        subplot(3,4,[9:12])
        for h = 13
            plot(out(h).prsvelocPattern(2:end , :)' , 'LineWidth' , 3)
            hold on
            plot(out(h).prsvelocPattern(1 , :)' , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Normalized time' , 'FontSize' , 10)
            ylabel('press/sec' , 'FontSize' , 10)
            title(['Press Velocity, H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            ax.YLim = [0 8];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
        end
        legend({ 'Structure 1','Structure 2' 'Random'} , 'FontSize' , 20)
        
        figure('color' , [1 1 1])
        h_counter = 1;
        for h = 1:8
            subplot(3,4,h_counter)
            plot(out(h).prsvelocPattern(2:end,:)' , 'LineWidth' , 3)
            hold on
            plot(out(h).prsvelocPattern(1,:)' , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Normalized time' , 'FontSize' , 10)
            ylabel('press/sec' , 'FontSize' , 10)
            title(['Press Velocity, H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            ax.YLim = [0 8];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
            h_counter = h_counter+1;
        end
        
        subplot(3,4,[9:12])
        for h = 13
            plot(out(h).prsvelocPattern(2:end , :)' , 'LineWidth' , 3)
            hold on
            plot(out(h).prsvelocPattern(1 , :)' , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Normalized time' , 'FontSize' , 10)
            ylabel('press/sec' , 'FontSize' , 10)
            title(['Press Velocity, H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            ax.YLim = [0 8];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
        end
        legend({ 'Structure 1','Structure 2' 'Random'} , 'FontSize' , 20)
        
        
        %% press velocity vs position
        
        figure('color' , [1 1 1])
        h_counter = 1;
        for h = 1:8
            subplot(3,4,h_counter)
            for s = 2:3
                plot(out(h).pressPattern(s,:) ,out(h).prsvelocPattern(s,:), 'LineWidth' , 3)
                hold on
            end
            hold on
            plot(out(h).EyePattern(1,:),out(h).prsvelocPattern(1,:) , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Press Position' , 'FontSize' , 10)
            ylabel('press/sec' , 'FontSize' , 10)
            title(['Press Vel vs. Pos, H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            %             ax.YLim = [0 8];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
            h_counter = h_counter+1;
        end
        
        h = 13;
        subplot(3,4,[9:12])
        for s = 2:3
            plot(out(h).pressPattern(s,:) ,out(h).prsvelocPattern(s,:), 'LineWidth' , 3)
            hold on
        end
        hold on
        plot(out(h).EyePattern(1,:),out(h).prsvelocPattern(1,:) , 'LineWidth' , 3 , 'color' , 'black')
        grid on
        xlabel('Press Position' , 'FontSize' , 10)
        ylabel('press/sec' , 'FontSize' , 10)
        title(['Press Vel vs. Pos, H =  ' , num2str(h)] )
        ax = gca;
        ax.FontSize = 10;
        ax.Box = 'off';
        
        %% Eye velocity vs position
        
        figure('color' , [1 1 1])
        h_counter = 1;
        for h = 1:8
            subplot(3,4,h_counter)
            for s = 2:3
                plot(out(h).EyePattern(s,:) ,out(h).eyevelocPattern(s,:), 'LineWidth' , 3)
                hold on
            end
            hold on
            plot(out(h).EyePattern(1,:),out(h).eyevelocPattern(1,:) , 'LineWidth' , 3 , 'color' , 'black')
            grid on
            xlabel('Press Position' , 'FontSize' , 10)
            ylabel('press/sec' , 'FontSize' , 10)
            title(['Press Vel vs. Pos, H =  ' , num2str(h)] )
            ax = gca;
            ax.FontSize = 10;
            ax.Box = 'off';
            %             ax.YLim = [0 8];
            %             ax.YLim = [min(min(out(h).EyePattern(: , 901:end))) 1.5+min(min(out(h).EyePattern(: , 901:end))) ]
            h_counter = h_counter+1;
        end
        
        h = 13;
        subplot(3,4,[9:12])
        for s = 2:3
            plot(out(h).EyePattern(s,:) ,out(h).eyevelocPattern(s,:), 'LineWidth' , 3)
            hold on
        end
        hold on
        plot(out(h).EyePattern(1,:),out(h).eyevelocPattern(1,:) , 'LineWidth' , 3 , 'color' , 'black')
        grid on
        xlabel('Press Position' , 'FontSize' , 10)
        ylabel('press/sec' , 'FontSize' , 10)
        title(['Press Vel vs. Pos, H =  ' , num2str(h)] )
        ax = gca;
        ax.FontSize = 10;
        ax.Box = 'off';
        
        
        %
    case 'crossvaldist_pos'
        %         if isequal(subjnum , 1:length(subj_name)-1)
        %             subjnum = 1:length(subj_name);
        %         end
        calcDist = 1;
        if calcDist
            if ~isempty(rep)
                switch rep
                    case 1
                        load([baseDir , '/se2_tnorm_r1.mat'])
                    case 2
                        load([baseDir , '/se2_tnorm_r2.mat'])
                    otherwise
                        load([baseDir , '/se2_tnorm.mat'])
                end
            else
                load([baseDir , '/se2_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [2:3];
            for hrzn = [1:8 ,13]
                for d  = [1 8 9]
                    dwithclass{d,hrzn} = [];
                    dbetclass{d,hrzn} = [];
                    for SubN = 1:length(subjnum)
                        out(hrzn, SubN, d).eye= N(hrzn).norm(d,subjnum(SubN)).eye;
                        ll = length(out(hrzn, SubN, d).eye);
                        if ll~=0
                            if ~isempty(out(hrzn,SubN, d).eye{1})
                                X = [];clear dis
                                lab = [];
                                for s = 1:length(seqs)
                                    X = [X ; out(hrzn,SubN, d).eye{seqs(s)}];
                                    lab = [lab ; s*ones(size(out(hrzn,SubN, d).eye{seqs(s)} , 1) ,1)];
                                end
                                for i = 1:size(X,1)
                                    id = ones(size(X , 1) , 1);
                                    id(i) = 0;
                                    Y = inpaint_nans(X(~id , :)); % The left out trial
                                    X1 = X(id==1 , :);     % all the other trials
                                    lab1 = lab(id==1);
                                    s1 = lab(id ==0);
                                    SS = seqs(seqs~= s1+1)-1;
                                    for s = SS
                                        if sum(lab1==s >=1)
                                            m = nanmean(X1(lab1==s , :),1);
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;[pdist([Y;m], distance) s1 s SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                        else
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;NaN s1 s SubN d ];
                                        end
                                    end
                                    if sum(lab1==s1 >=1)
                                        m = nanmean(X1(lab1==s1 , :),1);
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;[pdist([Y;m], distance) s1 s1 SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                    else
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;NaN s1 s1 SubN d ];
                                    end
                                end
                                alldist = [dwithclass{d,hrzn}(dwithclass{d,hrzn}(:,4) == SubN , :) ;dbetclass{d,hrzn}(dbetclass{d,hrzn}(:,4) == SubN , :)];
                                for l =  1:length(seqs)
                                    for l1 =  1:length(seqs)
                                        out(hrzn, SubN, d).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l)-1 & alldist(:,3)==seqs(l1)-1) , 1);
                                    end
                                end
                            end
                        end
                    end
                    if ~isempty(dwithclass{d,hrzn})
                        alldist = [dwithclass{d,hrzn} ;dbetclass{d,hrzn}];
                        alldist = [alldist [ones(length(dwithclass{d,hrzn}),1);zeros(length(dbetclass{d,hrzn}),1)]];
                        for l =  1:length(seqs)
                            for l1 =  1:length(seqs)
                                out(hrzn ,SubN+1, d).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                            end
                        end
                        eyesig(d,hrzn).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
                    end
                end
                
                allday_dwithclass{hrzn} = [];
                allday_dbetclass{hrzn} = [];
                for d = [1 8 9]
                    allday_dwithclass{hrzn} = [allday_dwithclass{hrzn} ; dwithclass{d,hrzn}];
                    allday_dbetclass{hrzn} = [allday_dbetclass{hrzn} ; dbetclass{d,hrzn}];
                end
                
                h1 = figure('color' , 'white');
                
                for SubN = 1:length(subjnum)
                    idx1 = allday_dwithclass{hrzn}(:,4) == SubN & ~isnan(allday_dwithclass{hrzn}(: , 1));
                    [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN , :),eERRORw{hrzn}(SubN , :)] = lineplot(allday_dwithclass{hrzn}(idx1, 5) , allday_dwithclass{hrzn}(idx1, 1));
                    hold on
                    idx1 = allday_dbetclass{hrzn}(:,4) == SubN & ~isnan(allday_dbetclass{hrzn}(: , 1));
                    [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN , :),eERRORb{hrzn}(SubN , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1, 1));
                end
                idx1 = ~isnan(allday_dwithclass{hrzn}(: , 1));
                [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN+1 , :),eERRORw{hrzn}(SubN+1 , :)] = lineplot(allday_dwithclass{hrzn}(idx1 , 5) , allday_dwithclass{hrzn}(idx1 , 1));
                hold on
                idx1 = ~isnan(allday_dbetclass{hrzn}(: , 1));
                [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN+1 , :),eERRORb{hrzn}(SubN+1 , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1 , 1));
                close(h1);
                
            end
            %% visualize EYE within class
            figure('color' , 'white')
            SubN = 1: length(subj_name);
            LineWidth = [.5*ones(1,length(subj_name)-1) , 3];
            colors = [repmat({'k'} , 1 , length(subj_name)-1) , {'r'}];
            counter = 1;
            
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Pos within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Pos within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            %% visualize EYE Between class
            figure('color' , 'white')
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Pos between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Pos between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                ax = gca;
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            clear dwithclass dbetclass eERRORb ePLOTb ePLOTw eERRORw allday_dwithclass allday_dbetclass
            seqs = [2:3];
            for hrzn = [1:8 ,13]
                for d = [1 8 9]
                    dwithclass{d,hrzn} = [];
                    dbetclass{d,hrzn} = [];
                    for SubN = 1:length(subjnum)
                        out(hrzn, SubN, d).press= N(hrzn).norm(d,subjnum(SubN)).press;
                        ll = length(out(hrzn, SubN, d).press);
                        if ll~=0
                            if ~isempty(out(hrzn,SubN, d).press{1})
                                X = [];clear dis
                                lab = [];
                                for s = 1:length(seqs)
                                    X = [X ; out(hrzn,SubN, d).press{seqs(s)}];
                                    lab = [lab ; s*ones(size(out(hrzn,SubN, d).press{seqs(s)} , 1) ,1)];
                                end
                                for i = 1:size(X,1)
                                    id = ones(size(X , 1) , 1);
                                    id(i) = 0;
                                    Y = inpaint_nans(X(~id , :)); % The left out trial
                                    X1 = X(id==1 , :);     % all the other trials
                                    lab1 = lab(id==1);
                                    s1 = lab(id ==0);
                                    SS = seqs(seqs~= s1+1)-1;
                                    for s = SS
                                        if sum(lab1==s >=1)
                                            m = nanmean(X1(lab1==s , :),1);
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;[pdist([Y;m], distance) s1 s SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                        else
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;NaN s1 s SubN d ];
                                        end
                                    end
                                    if sum(lab1==s1 >=1)
                                        m = nanmean(X1(lab1==s1 , :),1);
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;[pdist([Y;m], distance) s1 s1 SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                    else
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;NaN s1 s1 SubN d ];
                                    end
                                end
                                alldist = [dwithclass{d,hrzn}(dwithclass{d,hrzn}(:,4) == SubN , :) ;dbetclass{d,hrzn}(dbetclass{d,hrzn}(:,4) == SubN , :)];
                                for l =  1:length(seqs)
                                    for l1 =  1:length(seqs)
                                        out(hrzn, SubN, d).D_pres(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l)-1 & alldist(:,3)==seqs(l1)-1) , 1);
                                    end
                                end
                            end
                        end
                    end
                    if ~isempty(dwithclass{d,hrzn})
                        alldist = [dwithclass{d,hrzn} ;dbetclass{d,hrzn}];
                        alldist = [alldist [ones(length(dwithclass{d,hrzn}),1);zeros(length(dbetclass{d,hrzn}),1)]];
                        for l =  1:length(seqs)
                            for l1 =  1:length(seqs)
                                out(hrzn ,SubN+1, d).press(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                            end
                        end
                        prssig(d,hrzn).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
                    end
                end
                
                allday_dwithclass{hrzn} = [];
                allday_dbetclass{hrzn} = [];
                for d = [1 8 9]
                    allday_dwithclass{hrzn} = [allday_dwithclass{hrzn} ; dwithclass{d,hrzn}];
                    allday_dbetclass{hrzn} = [allday_dbetclass{hrzn} ; dbetclass{d,hrzn}];
                end
                
                h1 = figure('color' , 'white');
                
                for SubN = 1:length(subjnum)
                    idx1 = allday_dwithclass{hrzn}(:,4) == SubN & ~isnan(allday_dwithclass{hrzn}(: , 1));
                    [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN , :),eERRORw{hrzn}(SubN , :)] = lineplot(allday_dwithclass{hrzn}(idx1, 5) , allday_dwithclass{hrzn}(idx1, 1));
                    hold on
                    idx1 = allday_dbetclass{hrzn}(:,4) == SubN & ~isnan(allday_dbetclass{hrzn}(: , 1));
                    [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN , :),eERRORb{hrzn}(SubN , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1, 1));
                end
                idx1 = ~isnan(allday_dwithclass{hrzn}(: , 1));
                [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN+1 , :),eERRORw{hrzn}(SubN+1 , :)] = lineplot(allday_dwithclass{hrzn}(idx1 , 5) , allday_dwithclass{hrzn}(idx1 , 1));
                hold on
                idx1 = ~isnan(allday_dbetclass{hrzn}(: , 1));
                [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN+1 , :),eERRORb{hrzn}(SubN+1 , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1 , 1));
                close(h1);
            end
            
            %% visualize Press within class
            figure('color' , 'white')
            SubN = 1: length(subj_name);
            LineWidth = [.5*ones(1,length(subj_name)-1) , 3];
            colors = [repmat({'k'} , 1 , length(subj_name)-1) , {'r'}];
            counter = 1;
            
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Pos within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                legend({'W' , 'B'})
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Pos within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                legend({'W' , 'B'})
                hold on
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            %% visualize EYE Between class
            figure('color' , 'white')
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Pos between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Pos between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                ax = gca;
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            
        else
            seqs = [2:7];
            load([baseDir , '/eyepos_dbetclass.mat'])
            load([baseDir , '/eyepos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_eyepos_dbetclass.mat'])
                load([baseDir , '/EUC_eyepos_dwithclass.mat'])
            end
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prspos_dbetclass.mat'])
            load([baseDir , '/prspos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_prspos_dbetclass.mat'])
                load([baseDir , '/EUC_prspos_dwithclass.mat'])
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        end
    case 'crossvaldist_vel'
        
        %         if isequal(subjnum , 1:length(subj_name)-1)
        %             subjnum = 1:length(subj_name);
        %         end
        calcDist = 1;
        if calcDist
            if ~isempty(rep)
                switch rep
                    case 1
                        load([baseDir , '/se2_tnorm_r1.mat'])
                    case 2
                        load([baseDir , '/se2_tnorm_r2.mat'])
                    otherwise
                        load([baseDir , '/se2_tnorm.mat'])
                end
            else
                load([baseDir , '/se2_tnorm.mat'])
            end
            
            %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EYE
            seqs = [2:3];
            for hrzn = [1:8 ,13]
                for d  = [1 8 9]
                    dwithclass{d,hrzn} = [];
                    dbetclass{d,hrzn} = [];
                    for SubN = 1:length(subjnum)
                        out(hrzn, SubN, d).eyeveloc= N(hrzn).norm(d,subjnum(SubN)).eyeveloc;
                        ll = length(out(hrzn, SubN, d).eyeveloc);
                        if ll~=0
                            if ~isempty(out(hrzn,SubN, d).eyeveloc{1})
                                X = [];clear dis
                                lab = [];
                                for s = 1:length(seqs)
                                    X = [X ; out(hrzn,SubN, d).eyeveloc{seqs(s)}];
                                    lab = [lab ; s*ones(size(out(hrzn,SubN, d).eyeveloc{seqs(s)} , 1) ,1)];
                                end
                                for i = 1:size(X,1)
                                    id = ones(size(X , 1) , 1);
                                    id(i) = 0;
                                    Y = inpaint_nans(X(~id , :)); % The left out trial
                                    X1 = X(id==1 , :);     % all the other trials
                                    lab1 = lab(id==1);
                                    s1 = lab(id ==0);
                                    SS = seqs(seqs~= s1+1)-1;
                                    for s = SS
                                        if sum(lab1==s >=1)
                                            m = nanmean(X1(lab1==s , :),1);
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;[pdist([Y;m], distance) s1 s SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                        else
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;NaN s1 s SubN d ];
                                        end
                                    end
                                    if sum(lab1==s1 >=1)
                                        m = nanmean(X1(lab1==s1 , :),1);
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;[pdist([Y;m], distance) s1 s1 SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                    else
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;NaN s1 s1 SubN d ];
                                    end
                                end
                                alldist = [dwithclass{d,hrzn}(dwithclass{d,hrzn}(:,4) == SubN , :) ;dbetclass{d,hrzn}(dbetclass{d,hrzn}(:,4) == SubN , :)];
                                for l =  1:length(seqs)
                                    for l1 =  1:length(seqs)
                                        out(hrzn, SubN, d).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l)-1 & alldist(:,3)==seqs(l1)-1) , 1);
                                    end
                                end
                            end
                        end
                    end
                    %                     if ~isempty(dwithclass{d,hrzn})
                    %                         alldist = [dwithclass{d,hrzn} ;dbetclass{d,hrzn}];
                    %                         alldist = [alldist [ones(length(dwithclass{d,hrzn}),1);zeros(length(dbetclass{d,hrzn}),1)]];
                    %                         for l =  1:length(seqs)
                    %                             for l1 =  1:length(seqs)
                    %                                 out(hrzn ,SubN+1, d).D_eye(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                    %                             end
                    %                         end
                    %                         eyesig(d,hrzn).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
                    %                     end
                end
                
                allday_dwithclass{hrzn} = [];
                allday_dbetclass{hrzn} = [];
                for d = [1 8 9]
                    allday_dwithclass{hrzn} = [allday_dwithclass{hrzn} ; dwithclass{d,hrzn}];
                    allday_dbetclass{hrzn} = [allday_dbetclass{hrzn} ; dbetclass{d,hrzn}];
                end
                
                h1 = figure('color' , 'white');
                
                for SubN = 1:length(subjnum)
                    idx1 = allday_dwithclass{hrzn}(:,4) == SubN & ~isnan(allday_dwithclass{hrzn}(: , 1));
                    [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN , :),eERRORw{hrzn}(SubN , :)] = lineplot(allday_dwithclass{hrzn}(idx1, 5) , allday_dwithclass{hrzn}(idx1, 1));
                    hold on
                    idx1 = allday_dbetclass{hrzn}(:,4) == SubN & ~isnan(allday_dbetclass{hrzn}(: , 1));
                    [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN , :),eERRORb{hrzn}(SubN , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1, 1));
                end
                idx1 = ~isnan(allday_dwithclass{hrzn}(: , 1));
                [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN+1 , :),eERRORw{hrzn}(SubN+1 , :)] = lineplot(allday_dwithclass{hrzn}(idx1 , 5) , allday_dwithclass{hrzn}(idx1 , 1));
                hold on
                idx1 = ~isnan(allday_dbetclass{hrzn}(: , 1));
                [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN+1 , :),eERRORb{hrzn}(SubN+1 , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1 , 1));
                close(h1);
                
            end
            %% visualize EYE within class
            figure('color' , 'white')
            SubN = 1: length(subj_name);
            LineWidth = [.5*ones(1,length(subj_name)-1) , 3];
            colors = [repmat({'k'} , 1 , length(subj_name)-1) , {'r'}];
            counter = 1;
            
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Vel within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Vel within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            %% visualize EYE Between class
            figure('color' , 'white')
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Vel between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Eye Vel between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                ax = gca;
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRESS
            clear dwithclass dbetclass eERRORb ePLOTb ePLOTw eERRORw allday_dwithclass allday_dbetclass
            seqs = [2:3];
            for hrzn = [1:8 ,13]
                for d = [1 8 9]
                    dwithclass{d,hrzn} = [];
                    dbetclass{d,hrzn} = [];
                    for SubN = 1:length(subjnum)
                        out(hrzn, SubN, d).prsveloc= N(hrzn).norm(d,subjnum(SubN)).prsveloc;
                        ll = length(out(hrzn, SubN, d).prsveloc);
                        if ll~=0
                            if ~isempty(out(hrzn,SubN, d).prsveloc{1})
                                X = [];clear dis
                                lab = [];
                                for s = 1:length(seqs)
                                    X = [X ; out(hrzn,SubN, d).prsveloc{seqs(s)}];
                                    lab = [lab ; s*ones(size(out(hrzn,SubN, d).prsveloc{seqs(s)} , 1) ,1)];
                                end
                                for i = 1:size(X,1)
                                    id = ones(size(X , 1) , 1);
                                    id(i) = 0;
                                    Y = inpaint_nans(X(~id , :)); % The left out trial
                                    X1 = X(id==1 , :);     % all the other trials
                                    lab1 = lab(id==1);
                                    s1 = lab(id ==0);
                                    SS = seqs(seqs~= s1+1)-1;
                                    for s = SS
                                        if sum(lab1==s >=1)
                                            m = nanmean(X1(lab1==s , :),1);
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;[pdist([Y;m], distance) s1 s SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                        else
                                            dbetclass{d,hrzn} = [dbetclass{d,hrzn} ;NaN s1 s SubN d ];
                                        end
                                    end
                                    if sum(lab1==s1 >=1)
                                        m = nanmean(X1(lab1==s1 , :),1);
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;[pdist([Y;m], distance) s1 s1 SubN d]];%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                                    else
                                        dwithclass{d,hrzn} = [dwithclass{d,hrzn} ;NaN s1 s1 SubN d ];
                                    end
                                end
                                alldist = [dwithclass{d,hrzn}(dwithclass{d,hrzn}(:,4) == SubN , :) ;dbetclass{d,hrzn}(dbetclass{d,hrzn}(:,4) == SubN , :)];
                                for l =  1:length(seqs)
                                    for l1 =  1:length(seqs)
                                        out(hrzn, SubN, d).D_pres(l,l1) = nanmean(alldist(alldist(:,2) == seqs(l)-1 & alldist(:,3)==seqs(l1)-1) , 1);
                                    end
                                end
                            end
                        end
                    end
                    if ~isempty(dwithclass{d,hrzn})
                        alldist = [dwithclass{d,hrzn} ;dbetclass{d,hrzn}];
                        alldist = [alldist [ones(length(dwithclass{d,hrzn}),1);zeros(length(dbetclass{d,hrzn}),1)]];
                        for l =  1:length(seqs)
                            for l1 =  1:length(seqs)
                                out(hrzn ,SubN+1, d).prsveloc(l,l1) = nanmean(alldist(alldist(:,2) == l & alldist(:,3)==l1,1));
                            end
                        end
                        prssig(d,hrzn).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
                    end
                end
                
                allday_dwithclass{hrzn} = [];
                allday_dbetclass{hrzn} = [];
                for d = [1 8 9]
                    allday_dwithclass{hrzn} = [allday_dwithclass{hrzn} ; dwithclass{d,hrzn}];
                    allday_dbetclass{hrzn} = [allday_dbetclass{hrzn} ; dbetclass{d,hrzn}];
                end
                
                h1 = figure('color' , 'white');
                
                for SubN = 1:length(subjnum)
                    idx1 = allday_dwithclass{hrzn}(:,4) == SubN & ~isnan(allday_dwithclass{hrzn}(: , 1));
                    [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN , :),eERRORw{hrzn}(SubN , :)] = lineplot(allday_dwithclass{hrzn}(idx1, 5) , allday_dwithclass{hrzn}(idx1, 1));
                    hold on
                    idx1 = allday_dbetclass{hrzn}(:,4) == SubN & ~isnan(allday_dbetclass{hrzn}(: , 1));
                    [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN , :),eERRORb{hrzn}(SubN , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1, 1));
                end
                idx1 = ~isnan(allday_dwithclass{hrzn}(: , 1));
                [xcoordw{hrzn}(SubN , :),ePLOTw{hrzn}(SubN+1 , :),eERRORw{hrzn}(SubN+1 , :)] = lineplot(allday_dwithclass{hrzn}(idx1 , 5) , allday_dwithclass{hrzn}(idx1 , 1));
                hold on
                idx1 = ~isnan(allday_dbetclass{hrzn}(: , 1));
                [xcoordb{hrzn}(SubN , :),ePLOTb{hrzn}(SubN+1 , :),eERRORb{hrzn}(SubN+1 , :)] = lineplot(allday_dbetclass{hrzn}(idx1 , 5) , allday_dbetclass{hrzn}(idx1 , 1));
                close(h1);
            end
            
            %% visualize Press within class
            figure('color' , 'white')
            SubN = 1: length(subj_name);
            LineWidth = [.5*ones(1,length(subj_name)-1) , 3];
            colors = [repmat({'k'} , 1 , length(subj_name)-1) , {'r'}];
            counter = 1;
            
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Vel within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTw{hrzn}(sn , :),eERRORw{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Vel within class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            %% visualize Press Between class
            figure('color' , 'white')
            for hrzn  = 1:8
                subplot(3,4,hrzn)
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Vel between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            for hrzn  = 13
                subplot(3,4,[9:12])
                for sn = 1:length(SubN)
                    errorbar(ePLOTb{hrzn}(sn , :),eERRORb{hrzn}(sn , :) , 'LineWidth' , LineWidth(sn) , 'color' ,  colors{sn})
                    hold on
                end
                grid on
                title(['Press Vel between class dist, h = ' , num2str(hrzn) , ' day ' num2str(days{d})])
                hold on
                ax = gca;
                set(gca , 'XTick' ,[1:3] , 'XTickLabel' , {'Day 1' , 'Day 2,3' , 'Day 4,5'})
                ylabel('Average distance');
            end
            
            
        else
            seqs = [2:7];
            load([baseDir , '/eyepos_dbetclass.mat'])
            load([baseDir , '/eyepos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_eyepos_dbetclass.mat'])
                load([baseDir , '/EUC_eyepos_dwithclass.mat'])
            end
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                eyesig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,ePLOTw(SubN , :),eERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,ePLOTb(SubN , :),eERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,ePLOTw(SubN+1 , :),eERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,ePLOTb(SubN+1 , :),eERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
            
            for d = 1:4
                alldist = [dwithclass{d} ;dbetclass{d}];
                alldist = [alldist [ones(length(dwithclass{d}),1);zeros(length(dbetclass{d}),1)]];
                
                prssig(d).veldist = anovaMixed(alldist(:,1) ,alldist(:,4) ,'within',alldist(:,6),{'Within/Between'},'intercept',1)  ;
            end
            load([baseDir , '/prspos_dbetclass.mat'])
            load([baseDir , '/prspos_dwithclass.mat'])
            if isequall(distance , 'euclidean')
                load([baseDir , '/EUC_prspos_dbetclass.mat'])
                load([baseDir , '/EUC_prspos_dwithclass.mat'])
            end
            allday_dwithclass = [];
            allday_dbetclass = [];
            for d = 1:4
                allday_dwithclass = [allday_dwithclass ; dwithclass{d}];
                allday_dbetclass = [allday_dbetclass ; dbetclass{d}];
            end
            h1 = figure('color' , 'white');
            for SubN = 1:length(subjnum)
                [xcoord,pPLOTw(SubN , :),pERRORw(SubN , :)] = lineplot(allday_dwithclass(allday_dwithclass(:,4) == SubN , 5) , allday_dwithclass(allday_dwithclass(:,4) == SubN , 1));
                hold on
                [xcoord,pPLOTb(SubN , :),pERRORb(SubN , :)] = lineplot(allday_dbetclass(allday_dbetclass(:,4) == SubN , 5) , allday_dbetclass(allday_dbetclass(:,4) == SubN , 1));
            end
            [xcoord,pPLOTw(SubN+1 , :),pERRORw(SubN+1 , :)] = lineplot(allday_dwithclass(: , 5) , allday_dwithclass(: , 1));
            hold on
            [xcoord,pPLOTb(SubN+1 , :),pERRORb(SubN+1 , :)] = lineplot(allday_dbetclass(: , 5) , allday_dbetclass(: , 1));
            close(h1);
            
        end
    case 'crossval_IPI_dist'
%         h= input('Which Horizons?');
        horz = {[1] [2] [3] [4:13]};
        count = 1;
        for hr = 1:length(horz)
            clear X Y X1 d ANA ANAs lab nIdx 
            h = horz{hr};
            ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0:2]) & ~Dall.isError & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
            % time normalizing the IPIs (sp the press indecies) to 1 : 1000 normalized samples
            for tn = 1:length(ANA.TN)
                n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
                nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
                ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
            end
            X = ANA.IPI_norm;clear d
            lab = ANA.seqNumb +1;
            
            for i = 1:size(X,1)
                id = ones(size(X , 1) , 1);
                id(i) = 0;
                Y = inpaint_nans(X(~id , :));
                X1 = X(id==1 , :);
                lab1 = lab(id==1);
                for s = 1:3
                    m = nanmean(X1(lab1==s , :));
                    d{hr}(i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                end
                
            end
            seqs =  unique(lab);
            for l =  1:length(seqs)
                for l1 =  1:length(seqs)
                    id = lab == l;s
                    out.D_IPI{hr}(l,l1) = nanmean(d{hr}(id , l1));
                end
            end
            
            
            
            sn = unique(ANA.SN);
            for sub = 1:length(sn)
                ANAs = getrow(ANA , ANA.SN == sn(sub));
                X = ANAs.IPI;clear d
                
                lab = ANAs.seqNumb +1;
                
                for i = 1:size(X,1)
                    id = ones(size(X , 1) , 1);
                    id(i) = 0;
                    Y = inpaint_nans(X(~id , :));
                    X1 = X(id==1 , :);
                    lab1 = lab(id==1);
                    for s = 1:3
                        m = nanmean(X1(lab1==s , :));
                        ds(sub , i , s) = pdist([Y;m], distance);%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                    end
                    
                end
                seqs =  unique(lab);
                for l =  1:length(seqs)
                    for l1 =  1:length(seqs)
                        id = lab == l;
                        out.D_IPI_sub.sub(count,1) = sub;
                        out.D_IPI_sub.d(count,1) = nanmean(ds(sub , id , l1));
                        out.D_IPI_sub.s1(count,1) = l1;
                        out.D_IPI_sub.s2(count,1) = l;
                        out.D_IPI_sub.h(count,1) = hr;
                        out.D_IPI_sub.isWithin(count , 1) = (l == l1);
                        count = count + 1;
                    end
                end
            end
        end
        figure('color' , [1 1 1])
        fIND = {[1] [2] [3] [4:6]};
        for hr = 1:length(horz)
            h = horz{hr};
            subplot(2,3,fIND{hr})
            imagesc(out.D_IPI{hr} , [0.5 1]);
            if length(h)>1
                title(['IPI Dissimilarity Matrix - Horizon(s) ' , num2str(h(1)) , ' - ' ,num2str(h(end))])
            else
                title(['IPI Dissimilarity Matrix - Horizon(s) ' , num2str(h(1))])
            end
            hold on
            line([1.5 7.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            line([1.5 1.5] ,[1.5 7.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            ax = gca;
            ax.XTick = [1:3];
            ax.YTick = [1:3];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' , 'Structure 3' ,'Structure 4' , 'Structure 5' ,'Structure 6'}
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            axis square
            colorbar
            hold on
            
            axis square
            colorbar
        end
        
        
        h1 = figure
        [xw , pw , ew] = lineplot(out.D_IPI_sub.h, out.D_IPI_sub.d,'style_thickline' , ...
            'subset' ,out.D_IPI_sub.isWithin);
        [xb , pb , eb] = lineplot(out.D_IPI_sub.h, out.D_IPI_sub.d,'style_thickline' , ...
            'subset' ,~out.D_IPI_sub.isWithin);
        close(h1);
        figure('color' , 'white')
        errorbar(xw , pw , ew , 'LineWidth' , 3);
        hold on 
        errorbar(xb , pb , eb , 'LineWidth' , 3);
        set(gca , 'FontSize' , 20 , 'XTick' , [1: 4] , 'XTickLabels' , {'H = 1' , 'H = 2' , 'H = 3' , 'H = 4 - Full'},...
            'Box' , 'off' , 'XLim' , [0 5]);
        legend({'Within Structure' , 'Between Structres'})
        ylabel('Correlation Distance')
        title('Within Vs. Between Structure Crossvalidated Dissimilarity')
        grid on
        
    case 'Errors'
        for tn = 1:length(Dall.TN)
            Dall.MT(tn , 1) = Dall.AllPressTimes(tn , Dall.seqlength(tn)) - Dall.AllPressTimes(tn , 1);
        end
        err = [];
        err0 = [];
        proberr = [];
        proberr0 = [];
        figure('color' , 'white')
        for h = [1:8]
            errh = [];
            errh0 = [];
            
            for sub = subjnum
                errhs = [];
                errhs0 = [];
                ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [1:2]) & Dall.isError  & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
                ANA0 = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , 0) & Dall.isError  & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
                
                for tn  = 1:length(ANA.TN)
                    ANA.ChunkBndry(tn , :) = diff(ANA.ChnkArrang(tn,:));
                    a = find(ANA.ChunkBndry(tn , :));
                    ANA.ChunkBndry(tn , a-1) = 3;
                    ANA.ChunkBndry(tn , end) = 3;
                    ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                    ANA.ChunkBndry(tn , 1:2) = [-1 -1];  % dont account for the first and last sseqeuce presses
                    ANA.ChunkBndry(tn , end-1:end) = [-1 -1];% dont account for the first and last sseqeuce presses
                    ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,13);
                    a = find(ANA.AllPress(tn,:) ~= ANA.AllResponse(tn,:));
                    err = [err ; [a(1) ANA.ChnkPlcmnt(tn,a(1)), ANA.Horizon(tn) , ANA.SN(tn) ]];
                    errh = [errh ; [a(1) ANA.ChnkPlcmnt(tn,a(1)), ANA.Horizon(tn) , ANA.SN(tn) ]];
                    errhs = [errhs ; [a(1) ANA.ChnkPlcmnt(tn,a(1)), ANA.Horizon(tn) , ANA.SN(tn) ]];
                end
                
                for tn  = 1:length(ANA0.TN)
                    a = find(ANA0.AllPress(tn,:) ~= ANA0.AllResponse(tn,:));
                    err0 = [err0 ; [a(1) 0, ANA0.Horizon(tn) , ANA0.SN(tn) ]];
                    errh0 = [errh0 ; [a(1) 0, ANA0.Horizon(tn) , ANA0.SN(tn) ]];
                    errhs0 = [errhs0 ; [a(1) 0, ANA0.Horizon(tn) , ANA0.SN(tn) ]];
                end
                for p = 1:14
                    proberr  = [proberr  ; 100*sum(errhs(:,1) == p)/size(errhs,1) sub h p];
                    proberr0 = [proberr0 ; 100*sum(errhs0(:,1) == p)/size(errhs0,1) sub h p];
                end
            end
            subplot(2,4,h)
            ax = gca;
            hold on
            histogram(errh(:,2) , 'Normalization' , 'probability');
            ylabel('Probability');
            xlabel('Chunk placement')
            title(['Errors in different chunk placements, H = ' , num2str(h)])
            ax.YLim =[0: 1];
            ax.XTick = [1:3];
            ax.XTickLabel = {'First' , 'Middle' , 'Last'};
        end
        
        out.Po1 = anovan(err(:,2) , err(:,3), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon'});
        
        temp = pivottable(proberr(:,3) , proberr(:,4) , proberr(:,1) , 'nanmean');
        figure('color' , 'white')
        subplot(1,2,1)
        imagesc(temp , [0 50])
        colorbar
        hold on
        ax = gca;
        ax.XTick = [1:14];
        xlabel('Presses')
        ax.YTickLabel  = fliplr({ 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' });
        title('%Errors per press in chunked sequences')
        
        temp0 = pivottable(proberr0(:,3) , proberr0(:,4) , proberr0(:,1) , 'nanmean');
        subplot(1,2,2)
        imagesc(temp0 , [0 50])
        colorbar
        hold on
        ax = gca;
        ax.XTick = [1:14];
        xlabel('Presses')
        ax.YTickLabel  = fliplr({ 'H = 8' 'H = 7' 'H = 6' 'H = 5' 'H = 4' 'H = 3' 'H = 2' 'H = 1' });
        title('%Errors per press in Random sequences')
        
        %% error percetanges
        errp = [];
        
        for h = [1:8]
            for sub = subjnum
                ANA = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , [1:6])& ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
                ANA0 = getrow(Dall , ismember(Dall.SN , sub) & ismember(Dall.seqNumb , 0) & ismember(Dall.Rep , rep) & ismember(Dall.Day , days{day}) & ismember(Dall.Horizon , h));
                errp = [errp ; 100*sum(ANA.isError)/length(ANA.TN) h sub 1];
                errp = [errp ; 100*sum(ANA0.isError)/length(ANA0.TN) h sub 0];
            end
        end
        h1 = figure;
        [xcoor1,PLOT1,ERROR1] = lineplot(errp(:,2) , errp(:,1) , 'subset' , errp(:,4) ==1);
        hold on
        [xcoor0,PLOT0,ERROR0] = lineplot(errp(:,2) , errp(:,1) , 'subset' , errp(:,4) ==0);
        close(h1)
        
        out.EP = anovan(errp(:,1) , errp(:,[2,4]), 'display' , 'off' , 'model' , 'full' , 'varnames' , {'Horizon' , 'Rand/Chunk'});
        figure('color' , 'white')
        h1 = plotshade(xcoor1',PLOT1,ERROR1,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
        hold on
        h0 = plotshade(xcoor0',PLOT0,ERROR0,'transp' , .2 , 'patchcolor' , 'm' , 'linecolor' , 'm' , 'linewidth' , 3 , 'linestyle' , ':');
        
        title(['Error rate,  p(r/c) = ' , num2str(out.EP(2)) , '       No horizon effect, no interaction'])
        legend([h1,h0] , {'Chunked' , 'Randon'})
        hold on
        ax = gca;
        ax.XTick = [1:8];
        xlabel('Horizon')
        ylabel('Percent of Error trial');
        grid on
    case 'Sacc_Chunked'
        % The alalysis routine allocated the space to digits asymmetrically,
        % so that all the space to the left of a digit is allocated to that digit
        % p-1<x<p ---> allocated to digit p
        % if you want symmetrical digit allocation
        %         prompt = 'Which Horizon?';
        %         h = input(prompt);
        isSymmetric = 1;
        %         subjnum = subjnum(~ismember(subjnum , [2 4 6])); % Chao's eye data is fucked up!!!!
        
        days = {1 [2 3] [4 5]};
        for h  =[1:8 , 13]
            for d = 1:length(days)
                DigFix{h,d}  = [];
                PervBen{h,d} = [];
                isSacc{h,d}  = [];
                ANA = getrow(Dall ,ismember(Dall.seqNumb , [1:2]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & ismember(Dall.Day , days{d}) & ismember(Dall.Horizon , h));
                for tn  = 1:length(ANA.TN)
                    ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                    a = find(ANA.ChunkBndry(tn , :));
                    ANA.ChunkBndry(tn , a(2:end)-1) = 3;
                    ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                    ANA.ChunkBndry(tn , 1:3) = [-1 -1 -1];  % dont account for the first and last sseqeuce presses
                    ANA.ChunkBndry(tn , end-2:end) = [-1 -1 -1];% dont account for the first and last sseqeuce presses
                    ANA.DigFixDur(tn , :) = zeros(1 ,14);
                    
                    
                    if isSymmetric
                        window = 12;
                        ANA.EyeFixDigit{tn , 1} = ANA.xEyePosDigit{tn}((ANA.AllPressIdx(tn , 1)+2)-window :ANA.AllPressIdx(tn , ANA.seqlength(tn)) + 2 + window) .* ANA.SaccFlag{tn};
                        for p = 1:14
                            id = ANA.EyeFixDigit{tn , 1}<=p+.5 & ANA.EyeFixDigit{tn , 1}>p-.5;
                            %                         ANA.EyeFixDigit{tn , 1}(id) = p;
                            if sum(id)
                                ANA.DigFixWeight(tn , p) = (sum(abs(1- (ANA.EyeFixDigit{tn , 1}(id) - p)))/sum(id))*(sum(id)/500);
                            else
                                ANA.DigFixWeight(tn , p) = 0;
                            end
                        end
                    end
                    perv_Ben        = [1:14] - ANA.EyePressTimePos(tn , :);
                    perv_Ben(perv_Ben<=-3.5 | perv_Ben>-.9) = NaN;
                    PervBen{h,d} = [PervBen{h,d} ;[perv_Ben(ANA.ChunkBndry(tn , :) == 1)' ...
                        (d)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                        ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                        ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))...
                        ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 1)'))]];
                    
                    PervBen{h,d} = [PervBen{h,d} ;[perv_Ben(ANA.ChunkBndry(tn , :) == 2)' ...
                        (d)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                        2*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                        ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))...
                        ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 2)'))]];
                    
                    PervBen{h,d} = [PervBen{h,d} ;[perv_Ben(ANA.ChunkBndry(tn , :) == 3)' ...
                        (d)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))...
                        3*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))...
                        ANA.SN(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))...
                        ANA.seqNumb(tn)*ones(size(perv_Ben(ANA.ChunkBndry(tn , :) == 3)'))]];
                    
                    
                    
                    
                    
                    for p = 1:14
                        ANA.DigFixDur(tn , p) = 2*nansum(ANA.EyeFixDigit{tn} == p);
                    end
                    ANA.DigFixDur(tn , perv_Ben<=-3.5 | perv_Ben>-.9) = NaN;
                    DigFix{h,d} = [DigFix{h,d} ;[ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)' ...
                        (d)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                        ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                        ANA.SN(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))...
                        ANA.seqNumb(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 1)'))]];
                    
                    DigFix{h,d} = [DigFix{h,d} ;[ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)' ...
                        (d)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                        2* ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                        ANA.SN(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))...
                        ANA.seqNumb(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 2)'))]];
                    
                    DigFix{h,d} = [DigFix{h,d} ;[ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)' ...
                        (d)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))...
                        3*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))...
                        ANA.SN(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))...
                        ANA.seqNumb(tn)*ones(size(ANA.DigFixWeight(tn , ANA.ChunkBndry(tn , :) == 3)'))]];
                    
                    
                    
                    
                    
                    
                    isSacc{h,d} = [isSacc{h,d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 1))/sum(ANA.ChunkBndry(tn , :) == 1) ...
                        d 1 ANA.SN(tn)  ANA.seqNumb(tn)]];
                    
                    isSacc{h,d} = [isSacc{h,d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 2))/sum(ANA.ChunkBndry(tn , :) == 2) ...
                        d 2 ANA.SN(tn)  ANA.seqNumb(tn)]];
                    
                    isSacc{h,d} = [isSacc{h,d} ;[sum(ANA.isSaccWhilePress(tn , ANA.ChunkBndry(tn , :) == 3))/sum(ANA.ChunkBndry(tn , :) == 3) ...
                        d 3 ANA.SN(tn)  ANA.seqNumb(tn)]];
                    
                end
                
                out.FD(d) = anovaMixed(DigFix{h,d}(:,1) , DigFix{h,d}(:,4) ,'within',[DigFix{h,d}(:,3)],{'Chunk place'},'subset' , DigFix{h,d}(:,2) == d & ~isnan(DigFix{h,d}(:,1)), 'intercept',1);
                out.PB(d) = anovaMixed(PervBen{h,d}(:,1) , PervBen{h,d}(:,4) ,'within',[PervBen{h,d}(:,3)],{'Chunk place'},'subset' , PervBen{h,d}(:,2) == d & ~isnan(PervBen{h,d}(:,1)), 'intercept',1);
                out.iS(d) = anovaMixed(isSacc{h,d}(:,1) , isSacc{h,d}(:,4) ,'within',[isSacc{h,d}(:,3)],{'Chunk place'},'subset' , isSacc{h,d}(:,2) == d &~isnan(isSacc{h,d}(:,1)), 'intercept',1);
                if d == 1
                    all_ANA = ANA;
                else
                    all_ANA = addstruct(all_ANA , ANA);
                end
            end
            FD{h} = [];
            cp1{h} = [];
            for chnkprs = 1:3
                idoi = DigFix{h,end}(:,3) == chnkprs;
                leng(chnkprs) = length(DigFix{h,end}(idoi , 1));
                TEMP{chnkprs} = NaN * zeros(leng(chnkprs) ,d);
                TEMP{chnkprs}(:,d) =  DigFix{h,end}(idoi , 1);
                for dd =  1:length(days)
                    idtemp = DigFix{h,dd}(:,3) == chnkprs;
                    lengtemp(chnkprs) = length(DigFix{h,dd}(idtemp , 1));
                    TEMP{chnkprs}(1:lengtemp(chnkprs),dd) = DigFix{h,dd}(idtemp,1);
                end
                FD{h} = [FD{h};TEMP{chnkprs}];
                cp1{h} = [cp1{h};chnkprs*ones(length(TEMP{chnkprs}) , 1)];
            end
            
            PB{h} = [];
            cp2{h} = [];
            for chnkprs = 1:3
                idoi = PervBen{h,end}(:,3) == chnkprs;
                leng(chnkprs) = length(PervBen{h,end}(idoi , 1));
                TEMP{chnkprs} = NaN * zeros(leng(chnkprs) ,d);
                TEMP{chnkprs}(:,d) =  PervBen{h,end}(idoi , 1);
                for dd =  1:length(days)
                    idtemp = PervBen{h,dd}(:,3) == chnkprs;
                    lengtemp(chnkprs) = length(PervBen{h,dd}(idtemp , 1));
                    TEMP{chnkprs}(1:lengtemp(chnkprs),dd) = PervBen{h,dd}(idtemp,1);
                end
                PB{h} = [PB{h};TEMP{chnkprs}];
                cp2{h} = [cp2{h};chnkprs*ones(length(TEMP{chnkprs}) , 1)];
            end
            
            iS{h} = [];
            cp3{h} = [];
            for chnkprs = 1:3
                idoi = isSacc{h,end}(:,3) == chnkprs;
                leng(chnkprs) = length(isSacc{h,end}(idoi , 1));
                TEMP{chnkprs} = NaN * zeros(leng(chnkprs) ,d);
                TEMP{chnkprs}(:,d) =  isSacc{h,end}(idoi , 1);
                for dd =  1:length(days)
                    idtemp = isSacc{h,dd}(:,3) == chnkprs;
                    lengtemp(chnkprs) = length(isSacc{h,dd}(idtemp , 1));
                    TEMP{chnkprs}(1:lengtemp(chnkprs),dd) = isSacc{h,dd}(idtemp,1);
                end
                iS{h} = [iS{h};TEMP{chnkprs}];
                cp3{h} = [cp3{h};chnkprs*ones(length(TEMP{chnkprs}) , 1)];
            end
        end
        justFix = 0;
        cols = {'blue' , 'red' , 'c'};
        if justFix
            figcount = 1;
            figure('color' , 'white');
            for h = [1:8 , 13]
                htemp = figure;
                [xcoords,PLOT_df,ERROR_df] = lineplot(cp1{h} , FD{h},'plotfcn','nanmean','leg' ,{'Day1' , 'Day2,3' , 'Day4,5'},'markersize',3,'linewidth',3,'markertype' , {'o' , '*' , '+'},'markercolor' , {'blue' , 'red' , 'c'},...
                    'errorcolor' ,{'blue' , 'red' , 'c'} , 'linecolor' ,{'blue' , 'red' , 'c'} );
                close(htemp);
                subplot(3,3,figcount)
                
                chnkprs = 1;
                h1 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                hold on
                chnkprs = 2;
                h2 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                hold on
                chnkprs = 3;
                h3 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                
                legend([h1 h2 h3] ,{'Day1' , 'Day2,3' , 'Day4,5'})
                
                hold on
                set(gca ,'XTick', [1:3] ,'XTickLabel' ,{'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'} ,'XTickLabelRotation' , 45 ,'FontSize' , 16)
                grid on
                ylabel('Sec')
                title(['Digit Fix Duration - Horizon ',  num2str(h)])
                figcount = figcount +1;
            end
        else
            figcount = 1;
            figure('color' , 'white');
            for h = [1:8 , 13]
                htemp = figure;
                [xcoords,PLOT_df,ERROR_df] = lineplot(cp1{h} , FD{h},'plotfcn','nanmean','leg' ,{'Day1' , 'Day2,3' , 'Day4,5'},'markersize',3,'linewidth',3,'markertype' , {'o' , '*' , '+'},'markercolor' , {'blue' , 'red' , 'c'},...
                    'errorcolor' ,{'blue' , 'red' , 'c'} , 'linecolor' ,{'blue' , 'red' , 'c'} );
                close(htemp);
                subplot(3,3,figcount)
                
                chnkprs = 1;
                h1 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                hold on
                chnkprs = 2;
                h2 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                hold on
                chnkprs = 3;
                h3 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                
                legend([h1 h2 h3] ,{'Day1' , 'Day2,3' , 'Day4,5'})
                
                hold on
                set(gca ,'XTick', [1:3] , 'XTickLabelRotation' , 45 ,'FontSize' , 16)
                if ismember(figcount , [7:9])
                    set(gca ,'XTickLabel' ,{'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'} )
                end
                grid on
                ylabel('Sec')
                title(['Digit Fix Duration - Horizon ',  num2str(h)])
                figcount = figcount +1;
            end
            
            figcount = 1;
            figure('color' , 'white');
            for h = [1:8 , 13]
                htemp = figure;
                [xcoords,PLOT_df,ERROR_df] = lineplot(cp2{h} , -PB{h},'plotfcn','nanmean','leg' ,{'Day1' , 'Day2,3' , 'Day4,5' },'markersize',3,'linewidth',3,'markertype' , {'o' , '*' , '+'},'markercolor' , {'blue' , 'red' , 'c'},...
                    'errorcolor' ,{'blue' , 'red' , 'c'} , 'linecolor' ,{'blue' , 'red' , 'c'});
                close(htemp);
                subplot(3,3,figcount)
                
                chnkprs = 1;
                h1 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                hold on
                chnkprs = 2;
                h2 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                hold on
                chnkprs = 3;
                h3 = plotshade(xcoords',PLOT_df(chnkprs,:) , ERROR_df(chnkprs,:),'transp' , .2 , 'patchcolor' , cols{chnkprs} , 'linecolor' , cols{chnkprs} , 'linewidth' , 3 , 'linestyle' , ':');
                
                legend([h1 h2 h3] ,{'Day1' , 'Day2,3' , 'Day4,5'})
                
                hold on
                set(gca ,'XTick', [1:3] , 'XTickLabelRotation' , 45 ,'FontSize' , 16)
                if ismember(figcount , [7:9])
                    set(gca ,'XTickLabel' ,{'First Chunk Presses' , 'Middle Chunk Presses' , 'Last Chunk Presses'} )
                end
                grid on
                ylabel('Sec')
                title(['Fixation/Press distance - Horizon ',  num2str(h)])
                figcount = figcount +1;
            end
            
            
            
            
            
            % prepare the data for lineplot
            
            %             figCount = 1;
            %             a = tapply(all_ANA , {'seqNumb' , 'SN' , 'Day'},{'isSaccWhilePress','nanmean','name','meanisSaccWhilePress'} , {'ChunkBndry','nanmean','name','ChunkBndry'});
            %             figure('color' , 'white')
            %             for d = 1:3
            %                 for seqnum = 1:length(unique(a.seqNumb))
            %                     subplot(3,2,figCount)
            %                     imagesc([a.meanisSaccWhilePress(a.seqNumb == seqnum & a.Day == d ,3:11) ; NaN*ones(3,9) ;...
            %                         repmat(mean(a.meanisSaccWhilePress(a.seqNumb == seqnum & a.Day == d , 3:11))/max(mean(a.meanisSaccWhilePress(a.seqNumb == seqnum & a.Day == d , 3:11))) , 5, 1)]);
            %                     title(['Is making Saccade While Pressing - Structure ' , num2str(seqnum) , ', Day ' , num2str(d)])
            %                     hold on
            %                     ax = gca;
            %                     ax.YTick = [1:3,7];
            %                     ax.YTickLabel = {'Sub 1' 'Sub 2' 'Sub 3' , 'Sub Avg'};
            %                     ax.XTick = [1:8];
            %                     ax.XTickLabel = mean(a.ChunkBndry(a.seqNumb == seqnum , 4:11));
            %                     figCount = figCount+1;
            %                 end
            %             end
            %
            %             a = tapply(all_ANA , {'seqNumb' , 'SN' , 'Day'},{'DigFixWeight','nanmedian','name','meanDigFixDur'} , {'ChunkBndry','nanmean','name','ChunkBndry'});
            %             a.meanDigFixDur = a.meanDigFixDur ./ repmat(max(a.meanDigFixDur') , 14,1)';
            %             figCount = 1;
            %             figure('color' , 'white')
            %             for d = 1:3
            %                 for seqnum = 1:length(unique(a.seqNumb))
            %                     subplot(3,2,figCount)
            %                     imagesc([a.meanDigFixDur(a.seqNumb == seqnum & a.Day == d ,3:11) ; NaN*ones(3,9) ;...
            %                         repmat(nanmean(a.meanDigFixDur(a.seqNumb == seqnum & a.Day == d , 3:11))/max(nanmean(a.meanDigFixDur(a.seqNumb == seqnum & a.Day == d , 3:11))) , 5, 1)]);
            %                     title(['Normalized Fixation Duration - Structure ' , num2str(seqnum) , ', Day ' , num2str(d)])
            %                     hold on
            %                     ax = gca;
            %                     ax.YTick = [1:3,7];
            %                     ax.YTickLabel = {'Sub 1' 'Sub 2' 'Sub 3'  'Sub Avg'};
            %                     ax.XTick = [1:8];
            %                     ax.XTickLabel = mean(a.ChunkBndry(a.seqNumb == seqnum , 4:11));
            %                     figCount = figCount+1;
            %                 end
            %             end
        end
        
        disp(['Day 1 FixDuration p-val = ' , num2str(out.FD(1).eff(2).p)])
        disp(['Day 2 FixDuration p-val = ' , num2str(out.FD(2).eff(2).p)])
        disp(['Day 3 FixDuration p-val = ' , num2str(out.FD(3).eff(2).p)])
        %         disp(['All Day FixDuration p-val = ' , num2str(out.FD(4).eff(2).p)])
        
        
        
        
        
        disp(['Day 1 isSacc p-val = ' , num2str(out.iS(1).eff(2).p)])
        disp(['Day 2 isSacc p-val = ' , num2str(out.iS(2).eff(2).p)])
        disp(['Day 3 isSacc p-val = ' , num2str(out.iS(3).eff(2).p)])
        %         disp(['All Day isSacc p-val = ' , num2str(out.iS(4).eff(2).p)])
        
        
        
        
        
        disp(['Day 1 PervBenefit p-val = ' , num2str(out.PB(1).eff(2).p)])
        disp(['Day 2 PervBenefit p-val = ' , num2str(out.PB(2).eff(2).p)])
        disp(['Day 3 PervBenefit p-val = ' , num2str(out.PB(3).eff(2).p)])
        %         disp(['All Day PervBenefit p-val = ' , num2str(out.PB(4).eff(2).p)])
    case 'Sacc_All'
        horz = input('Which horizon(s)?');
        if calc
            isSymmetric = 1;
%             subjnum = subjnum(~ismember(subjnum , 2)); % Chao's eye data is fucked up!!!!
            eyeinfo.PB         = [];
            eyeinfo.CB         = [];
            eyeinfo.Hor        = [];
            eyeinfo.sn         = [];
            eyeinfo.day        = [];
            eyeinfo.BN         = [];
            eyeinfo.sacPerSec  = [];
            eyeinfo.sacDur     = [];
            eyeinfo.sacPeakVel = [];
            eyeinfo.sacAmp     = [];
            eyeinfo.seqNumb    = [];
            eyeinfo.DigFixDur  = [];
            eyeinfo.prsnumb    = [];
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:2]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & cellfun(@length , Dall.xEyePosDigit)>1);
            ignorDig = 2;
            for tn  = 1:length(ANA.TN)
                ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a(ignorDig:end)-1) = 3;
                ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                ANA.ChunkBndry(tn , 1:ignorDig+1) = [-1 -1 -1];  % dont account for the first and last sseqeuce presses
                ANA.ChunkBndry(tn , end-ignorDig:end) = [-1 -1 -1];% dont account for the first and last sseqeuce presses
                ANA.DigFixWeighted(tn , :) = zeros(1 ,14);
                
                
                if isSymmetric
                    window = 50; %samples = 100ms
                    for p = 1:14
                        id = ANA.xEyePosDigit{tn , 1}<=p+.5 & ANA.xEyePosDigit{tn , 1}>p-.5;
                        if sum(id)
                            ANA.DigFixWeighted(tn , p) = mean(abs(ANA.xEyePosDigit{tn , 1}(id) - p))*(sum(id)/500)*1000;
                        else
                            ANA.DigFixWeighted(tn , p) = 0;
                        end
                        id = [ANA.AllPressIdx(tn , p) - window :ANA.AllPressIdx(tn , p) + window];
                        if id(1) > length(ANA.xEyePosDigit{tn}) | sum(id<0)>0
                             ANA.EyePressTimePos(tn , p) = NaN;
                        elseif id(end)>length(ANA.xEyePosDigit{tn})
                            ANA.EyePressTimePos(tn , p) = nanmedian(ANA.xEyePosDigit{tn}(id(1):end));
                        else
                            ANA.EyePressTimePos(tn , p) = nanmedian(ANA.xEyePosDigit{tn}(id));
                        end
                    end
                end
                perv_Ben           = [1:14] - ANA.EyePressTimePos(tn , :);
                goodid             = (perv_Ben>=-3.5 & perv_Ben<=-0.5);
                prsnumb            = find(goodid);
                count              = sum(goodid);
                eyeinfo.PB         = [eyeinfo.PB ;perv_Ben(goodid)'];
                eyeinfo.prsnumb    = [eyeinfo.prsnumb ;find(goodid')];
                eyeinfo.CB         = [eyeinfo.CB ;ANA.ChunkBndry(tn ,goodid)'];
                eyeinfo.Hor        = [eyeinfo.Hor ; ANA.Horizon(tn)*ones(count , 1)];
                eyeinfo.sn         = [eyeinfo.sn ; ANA.SN(tn)*ones(count , 1)];
                eyeinfo.day        = [eyeinfo.day ; ANA.Day(tn)*ones(count , 1)];
                eyeinfo.BN         = [eyeinfo.BN ; ANA.BN(tn)*ones(count , 1)];
                eyeinfo.sacPerSec  = [eyeinfo.sacPerSec ; ANA.SaccPerSec(tn)*ones(count , 1)];
                eyeinfo.sacDur     = [eyeinfo.sacDur ; mean(ANA.SaccDuration{tn})*ones(count , 1)];
                eyeinfo.sacPeakVel = [eyeinfo.sacPeakVel ; mean(ANA.SaccPeakVel{tn})*ones(count , 1)];
                eyeinfo.sacAmp     = [eyeinfo.sacAmp ; mean(ANA.SaccAmplitude{tn})*ones(count , 1)];
                eyeinfo.seqNumb    = [eyeinfo.seqNumb ; ANA.seqNumb(tn)*ones(count , 1)];
                eyeinfo.DigFixDur  = [eyeinfo.DigFixDur ;ANA.DigFixWeighted(tn ,goodid)'];
            end
            K_withinChunk = tapply(eyeinfo , {'day' , 'Hor' , 'sn','CB'} , {'sacDur' , 'nanmedian'}, {'DigFixDur' , 'nanmedian'},...
                {'sacAmp' , 'nanmedian'} , {'sacPerSec' , 'nanmedian'}, {'PB' , 'nanmedian'} , 'subset' ,...
                eyeinfo.seqNumb ~= 0 & eyeinfo.CB~=-1);
            
            K_sqnum = tapply(eyeinfo , {'day' , 'Hor' , 'sn','seqNumb' , 'prsnumb'} , {'sacDur' , 'nanmedian'}, {'sacPeakVel' , 'nanmedian'},...
                {'sacAmp' , 'nanmedian'} , {'sacPerSec' , 'nanmedian'}, {'PB' , 'nanmedian'},{'DigFixDur' , 'nanmedian'});
            
            K_rand_chnk = tapply(eyeinfo , {'day' , 'Hor' , 'sn','seqNumb'} , {'sacDur' , 'nanmedian'}, {'sacPeakVel' , 'nanmedian'},...
                {'sacAmp' , 'nanmedian'} , {'sacPerSec' , 'nanmedian'}, {'PB' , 'nanmedian'},{'DigFixDur' , 'nanmedian'},...
                'subset' , eyeinfo.CB~=-1);
            
            save([baseDir , '/se2_eyeInfo.mat'] , 'eyeinfo' ,'K_withinChunk' , 'K_sqnum' , 'K_rand_chnk' , '-v7.3')
        else
            load([baseDir , '/se2_eyeInfo.mat'])
        end
        ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:2]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & cellfun(@length , Dall.xEyePosDigit)>1);
        h = figure;
        dayz = {[1] , [2 3] , [4 5]};
        for sqn = 0:2
            for dd = 1:length(dayz)
                [xDur{sqn+1,dd},pDur{sqn+1,dd},eDur{sqn+1,dd}] = lineplot( K_sqnum.prsnumb ,  K_sqnum.DigFixDur , 'plotfcn' , 'nanmean','subset' ,...
                    ismember(K_sqnum.Hor , horz) & ismember(K_sqnum.seqNumb , sqn) & ismember(K_sqnum.day , dayz{dd}) , 'style_shade');
            end
            
        end
        close(h)
        
        h = figure('color' , 'white');
        
        for sqn = 0:2
            chnkpl = unique(ANA.ChnkPlcmnt(ANA.seqNumb == sqn , :) , 'rows');
            subplot(length(dayz) , 1 , sqn+1);
            for d = 1:length(dayz)
                errorbar(xDur{sqn+1,d}',pDur{sqn+1,d},eDur{sqn+1,d}, 'LineWidth' , 3 , 'color' , colors(d, :));
                hold on
                
            end
            firsts = find(chnkpl == 1);
            for j = 1:length(firsts)
                line([firsts(j) firsts(j)] , [0 250] , 'LineWidth'  ,2 , 'color' , 'black')
            end
            
            grid on
            
            set(gca , 'YLim' , [0 250] , 'XTick' , [2:14] , 'XLim' , [1 15])
            legend({'day 1' , 'days 2 3' , 'days 4 5'})
            title(['Fixation Duration on Digits on stucure ' , num2str(sqn) , ', in Horizon(s) '  , num2str(horz)])
        end
        
       

        
        K1 = K_rand_chnk;
        h = figure('color' , 'white');
        subplot(411)
        id = K1.day == 1 & ismember(K1.Hor , horz);
        K0.day = K1.day(id) ;
        K0.DigFixDur = K1.DigFixDur(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'DigFixDur' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) & ismember(K1.Hor , horz));
            K00.day = d*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.DigFixDur , K00.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  K0.DigFixDur , 'plotfcn' , 'nanmean','style_thickline');
        
        grid on
        set(gca , 'FontSize' , 20)
        title(['Weighted Average Fixation Duration in Horizon(s) '  , num2str(horz) , ' - Days 1, [2 3] , [4 5]'])
        ylabel('Wieghted ms')
        xlabel('Structure number')
        legend({'Pre-train' , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        %========================================= Preview
        subplot(412)
        id = K1.day == 1 & ismember(K1.Hor , horz);
        K0.day = K1.day(id) ;
        K0.PB = K1.PB(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'PB' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) & ismember(K1.Hor , horz));
            K00.day = d*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
             templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.PB , K00.sn,'within',templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  -K0.PB , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Looking Ahead in Horizon(s) '  , num2str(horz) , ' - Days 1, [2 3] , [4 5]'])
        ylabel('Digits')
        xlabel('Structure number')
        legend({'Pre-train' , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        %========================================= Saccade per second
        subplot(413)
        id = K1.day == 1 & ismember(K1.Hor , horz);
        K0.day = K1.day(id) ;
        K0.sacPerSec = K1.sacPerSec(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'sacPerSec' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) & ismember(K1.Hor , horz));
            K00.day = d*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.sacPerSec , K00.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  K0.sacPerSec , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Saccade rate in Horizon(s) '  , num2str(horz) , ' - Days 1, [2 3] , [4 5]'])
        xlabel('Structure number')
        ylabel('Saccade per sec')
        legend({'Pre-train' , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        %========================================= Saccade amplitude
        subplot(414)
        id = K1.day == 1 & ismember(K1.Hor , horz);
        K0.day = K1.day(id) ;
        K0.sacAmp = K1.sacAmp(id);
        K0.seqNumb  = K1.seqNumb(id);
        K0.sn = K1.sn(id);
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'seqNumb' , 'sn'} , {'sacAmp' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) & ismember(K1.Hor , horz));
            K00.day = d*ones(size(K00.seqNumb));
            K0 = addstruct(K0 , K00);
            templab = K00.seqNumb >0; % just test Random vs all structured
            temp = anovaMixed(K00.sacAmp , K00.sn,'within', templab,{'seqNumb'},'intercept',1)  ;
            out.rand_chnk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.seqNumb] ,  K0.sacAmp , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Sccade Amplitude in Horizon(s) '  , num2str(horz) , ' - Days 1, [2 3] , [4 5]'])
        xlabel('Structure number')
        legend({'Pre-train' , ['Day 2 P = ' , num2str(out.rand_chnk(2))], ['Day 3 P = ' , num2str(out.rand_chnk(3))]})
        ylabel('Digits')
       
        
        
        K1 = K_withinChunk;
        % within chunked
        h = figure('color' , 'white');
        subplot(211)
        id = K1.day == 1 & ismember(K1.Hor , horz);
        K0.day = K1.day(id) ;
        K0.DigFixDur = K1.DigFixDur(id);
        K0.CB  = K1.CB(id);
        K0.sn = K1.sn(id);
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'CB' , 'sn'} , {'DigFixDur' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) & ismember(K1.Hor , horz));
            K00.day = d*ones(size(K00.CB));
            K0 = addstruct(K0 , K00);
            templab = K00.CB ~= 2; % just test middle vs rest
            temp = anovaMixed(K00.DigFixDur , K00.sn,'within', templab ,{'seqNumb'},'intercept',1)  ;
            out.withinChunk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.CB] ,  K0.DigFixDur , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Average Fixation Duration in Horizon(s) '  , num2str(horz) , ' - Days 1, [2 3] , [4 5]'])
        xlabel('Chunk Placement')
        legend({'-' , num2str(out.withinChunk(2)), num2str(out.withinChunk(3))})
        out.withinChunk
        %========================================= Preview
         subplot(212)
        id = K1.day == 1 & ismember(K1.Hor , horz);
        K0.day = K1.day(id) ;
        K0.PB = K1.PB(id);
        K0.CB  = K1.CB(id);
        K0.sn = K1.sn(id);
        for d = 2:length(dayz)
            K00 = tapply(K1 , {'CB' , 'sn'} , {'PB' , 'nanmean'} , 'subset' , ismember(K1.day , dayz{d}) & ismember(K1.Hor , horz));
            K00.day = d*ones(size(K00.CB));
            K0 = addstruct(K0 , K00);
            templab = K00.CB;% ~= 3; % just test middle vs rest
            temp = anovaMixed(K00.PB , K00.sn,'within', templab ,{'seqNumb'},'intercept',1)  ;
            out.withinChunk(d) = temp.eff(2).p;
        end
        % between random and chunked
        [xDur,pDur,eDur] = lineplot([K0.day , K0.CB] ,  -K0.PB , 'plotfcn' , 'nanmean','style_thickline');
        grid on
        set(gca , 'FontSize' , 20)
        title(['Preview in Horizon(s) '  , num2str(horz) , ' - Days 1, [2 3] , [4 5]'])
        xlabel('Chunk Placement')
        legend({'-' , num2str(out.withinChunk(2)), num2str(out.withinChunk(3))})


    case 'glm_IPIs'
        N1 = input('Use Conditional Transition Probabilities? (y/n)' , 's');
        switch N1
            case 'y'
                cond = 1;
            otherwise
                cond = 0;
        end
        N2 = input('Use the last IPI for 2nd/3rd order transitions? (y/n)' , 's');
        switch N2
            case 'y'
                LastIPI = 1;
            otherwise
                LastIPI = 0;
        end
        
        if GroupCode == 1
            load([baseDir , '/CMB_34_1.mat'])
            CMB = CMB_34_1;
        elseif GroupCode == 1
            load([baseDir , '/CMB_34_2.mat'])
            CMB = CMB_34_2;
        end
        N3 = input('look into Chunked/Random/All sequences? (c/r/a)'  , 's');
        N4 = input('Use normalized IPIs? (Y/N)'  , 's');
        N5 = input('Analyse Medians, or raw IPIs? (M/R)' , 's');
        
        %%
        
        switch N3
            case 'c'
                M = C;
                titleSuffix = 'Chunked';
                hor = {1,2,3,4,5,6,7,8,13,[5:13]};
                legHor = {'H = 1', 'H = 2','H = 3','H = 4','H = 5','H = 6','H = 7','H = 8','H = 13' , 'H = 5-13'};
            case 'r'
                M = R;
                titleSuffix = 'Random';
                hor = {1,2,3,4,5,6,7,8,13,[5:13]};
                legHor = {'H = 1', 'H = 2','H = 3','H = 4','H = 5','H = 6','H = 7','H = 8','H = 13' , 'H = 5-13'};
            case 'a'
                M = All;
                titleSuffix = 'All Sequences';
        end
        
        
        for dd = 1:5
            
            for sn = 1:length(subj_name) - 1
                for h = 1:length(hor)
                    M.IPIarrangement(M.IPIarrangement == 0) = 1;
                    clear T TSS RSS0 RSS1 FSS1 FSS0 L X X1 X2 X3 X4 X5 X6
                    T = getrow(M , ismember(M.SN , sn) & ismember(M.Horizon , hor{h}) & ismember(M.Day , dd));
                    
                    
                    L = length(T.IPI);
                    X1 = ones(L , 1); % intercept
                    X2 = T.IPIarrangement;
                    %             X2 = T.estIPIarrangement;
                    
                    X2(X2==1) = -1; % between
                    X2(X2==0) = -1; % Random = between
                    X2(X2==2) = 1;  % within
                    
                    
                    
                    switch cond
                        case 1
                            X3 = T.t2Rank_n_binned;
                            switch LastIPI
                                case 1
                                    %                                 X4 = T.t3Prob_n(:,1);
                                    %                                 X5 = T.t4Prob_n(:,1);
                                    X4 = T.t3Rank_n_binned;
                                    X5 = T.t4Rank_n_binned;
                                otherwise
                                    X4 = mean(T.t3Prob_n , 2);
                                    X5 = mean(T.t4Prob_n , 2);
                            end
                        case 0
                            X3 = T.t2Prob;
                            switch LastIPI
                                case 1
                                    X4 = T.t3Prob(:,1);
                                    X5 = T.t4Prob(:,1);
                                otherwise
                                    X4 = mean(T.t3Prob , 2);
                                    X5 = mean(T.t4Prob , 2);
                            end
                    end
                    
                    X6 = max(T.BN) - [T.BN] +1;
                    X7 = ismember(T.t2 , [21:25])  +1;
                    switch N4
                        case {'n' 'N'}
                            X = [X1 X2 X3 X4 X5 X6 X7];
                            xx    = {[1 6:7] , [1 2 6:7] ,  [1 3 6:7]  , [1,3:4,6:7]           ,[1 3:7] ,      [1:3,6:7]  ,        [1:4,6:7]          [1:7]};
                            Y = T.IPI;
                            label = {'I+L' '  I+C+L', 'I+1st+L' ,'I+1st+2nd+L'    'I+1st+2nd+3rd+L' , 'I+C+1st+L'  ,  'I+C+1st+2nd+L'  ,  'Full'};
                        otherwise
                            %                         X = [X1 X2 X3 X4 X5];
                            %                         xx    = {[1 2] ,  [1 3]  , [1,3:4]           ,[1 3:5] ,      [1:3]  ,        [1:4]          [1:5]};
                            %                         Y = T.IPI_norm;
                            %                         label = {'I+C', 'I+1st' ,'I+1st+2nd'    'I+1st+2nd+3rd' , 'I+C+1st'  ,  'I+C+1st+2nd'  ,  'Full'};
                            X = [X1 X2 X3 X4 X5 X6 X7];
                            xx    = {[1 6:7] , [1 2 6:7] ,  [1 3 6:7]  , [1,3:4,6:7]           ,[1 3:7] ,      [1:3,6,7]  ,        [1:4,6,7]          [1:7]};
                            Y = T.IPI_norm;
                            label = {'I+L' '  I+C+L', 'I+1st+L' ,'I+1st+2nd+L'    'I+1st+2nd+3rd+L' , 'I+C+1st+L'  ,  'I+C+1st+2nd+L'  ,  'Full'};
                    end
                    Xnew = [X1];
                    B = pinv(Xnew'*Xnew)*Xnew'*Y;
                    Ypred = Xnew*B;
                    Res  = Y-Ypred;
                    TSS = sum(Y.^2); % Total Variance
                    FSS1 = sum(Ypred.^2); % Fitted Variance of the Null Model (just the intercept)
                    SSR1 = sum((Y-Ypred).^2); % Residual Variance of the Null Model
                    for k = 1:length(xx)
                        Xnew = X(:,xx{k});
                        Bnew = pinv(Xnew'*Xnew)*Xnew'*Y;
                        Ypred_new = Xnew*Bnew;
                        FSS0(k) = sum(Ypred_new.^2);  % Fitted Variance of the partial model
                        SSR0(k) = sum((Y-Ypred_new).^2); % Residual Variance of the partial Model
                    end
                    % ___________________________      R_squared = 1 - (Residual variance of the partial model/Residual variance of the null model)
                    R2{h,dd}(sn , :) = 1 - (SSR0./SSR1);
                    xlab{h,dd}(sn , :) = 1:length(xx);
                    % ___________________________
                end
            end
        end
        clear ePLOT xcoord ERROR
        figure('color' , 'white')
        for dd = 1:5
            subplot(2,5,dd)
            cCount = 1;
            for h = 1:length(hor)-1
                f = figure;
                [xcoordh,ePLOTh,ERRORh] = lineplot(reshape(xlab{h,dd} , numel(xlab{h,dd}) , 1) , reshape(R2{h,dd} , numel(R2{h,dd}) , 1) , 'plotfcn' , 'nanmean',...
                    'linecolor' , colors(cCount,:) , 'markercolor' , colors(cCount,:) , 'errorcolor' , colors(cCount,:) , 'linewidth' , 3);
                close(f)
                errorbar(xcoordh,ePLOTh,ERRORh , 'color' , colors(cCount,:) , 'LineWidth' , 3)
                hold on
                %             plotshade(xcoord',ePLOT , ERROR,'transp' , .2 , 'patchcolor' , colors(cCount,:) , 'linecolor' , colors(cCount,:) , 'linewidth' , 3 , 'linestyle' , ':');
                cCount = cCount+1;
            end
            ylabel('R^2')
            set(gca , 'XLim' , [0 length(xx)+1],'XTick' , [1: length(xx)] , 'XTickLabels' ,[],'FontSize' , 20,...
                'YLim' , [0 .4],'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' , Days ' , num2str(dd)])
            grid on
        end
        legend(legHor(1:end-1), 'Box' , 'off')
        for dd = 1:5
            subplot(2,5,5+dd)
            cCount = 1;
            for h = length(hor)
                f = figure;
                [xcoord(dd,:),ePLOT(dd,:),ERROR(dd,:)] = lineplot(reshape(xlab{h,dd} , numel(xlab{h,dd}) , 1) , reshape(R2{h,dd} , numel(R2{h,dd}) , 1) , 'plotfcn' , 'nanmean',...
                    'linecolor' , colors(cCount,:) , 'markercolor' , colors(cCount,:) , 'errorcolor' , colors(cCount,:) , 'linewidth' , 3);
                close(f)
                errorbar(xcoord(dd,:)',ePLOT(dd,:),ERROR(dd,:) , 'color' , colors(cCount,:) , 'LineWidth' , 3)
                hold on
                plotshade(xcoord(dd,:),ePLOT(dd,:) , ERROR(dd,:),'transp' , .2 , 'patchcolor' , colors(cCount,:) , 'linecolor' , colors(cCount,:) , 'linewidth' , 3 , 'linestyle' , ':');
                cCount = cCount+1;
            end
            ylabel('R^2')
            title([titleSuffix, ' , Days ' , num2str(dd)])
            set(gca , 'XLim' , [0 length(xx)+1] , 'XTick' , [1: length(xx)] , 'XTickLabels' , label , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'YLim' , [0 .3],'Box' , 'off' , 'GridAlpha' , 1)
            grid on
        end
        legend(legHor(end) , 'Box' , 'off')
        
        figure('color' , 'white')
        imagesc(ePLOT)
        colorbar
        set(gca , 'XTick' , [1 :length(xx)] , 'XTickLabels' , label , 'FontSize' , 20 ,...
            'XTickLabelRotation',45,'YTick' , [1 :5])
    case 'crossvaldist_chunk'
        %% chunk distances
        h = input('Which horizon?');
        ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:2]) & ismember(Dall.SN , subjnum) & Dall.isgood  & ~Dall.isError & ismember(Dall.Horizon ,h));
        % time normalizing the IPIs (sp the press indecies) to 1 : 1000 normalized samples
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , ANA.seqlength(tn))  - ANA.AllPressIdx(tn , 1)) / 500;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            IPI(tn , :) = detrend(diff(nIdx(tn ,:) , 1 , 2) , 'linear' , 2);
        end
        
        for tn = 1:length(ANA.TN)
            thresh = .3 * std(IPI(tn , :));
            [dum , estChnkBndry] = findpeaks(IPI(tn , :));% ,'MinPeakProminence', thresh); % slope of IPIs
            
            if ~isempty(estChnkBndry)
                goodpeak = ones(1, length(estChnkBndry));
                for cb = 1:length(estChnkBndry)
                    if IPI(tn , estChnkBndry(cb)) < nanmean(IPI(tn  ,:)) + thresh
                        goodpeak(cb) = 0;
                    end
                end
                if sum(goodpeak)
                    estChnkBndry = estChnkBndry(logical(goodpeak));
                else
                    estChnkBndry = [];
                end
            end
            
            
            ANA.estChnkBndry(tn , :) = zeros(1, 14);  % first presses of chunks will be 1
            ANA.estChnkBndry(tn , estChnkBndry+1) = 1;
            ANA.estChnkBndry(tn , 1) = 1;
        end
        
        ChnkArrang = [4 4 3 3 ; 3 4 3 4];
        for j = 1:size(ChnkArrang , 1)
            temp = [];
            temp1 = [];
            for k = 1:length(find(ChnkArrang(j , :)))
                temp = [temp k*ones(1,ChnkArrang(j,k))];
                temp1 = [temp1 1:ChnkArrang(j,k)];
            end
            ChnkArrng(j , :) = temp;
            ChnkPlcmnt(j,:)  = temp1;
        end
        IPIarrangement = diff(ChnkArrng , 1 , 2);   % between IPIs will have the index of 1
        IPIarrangement(~IPIarrangement) = 2;         % within IPIs will have the index of 2
        
        
        %% the boundry distance
        temp = diff(ChnkPlcmnt,1,2);
        temp(temp<0) = 0;
        chbndry = [ones(2,1) ~temp]; % all the first presses are one
        
        
        
        for s = 0:2
            A = getrow(ANA , ANA.seqNumb == s);
            CBD{s+1} = A.estChnkBndry;
        end
        diag_CnkBndry_es_es = [];
        offdiag_CnkBndry_es_es = [];
        
        diag_CnkBndry_es_ac = [];
        offdiag_CnkBndry_es_ac = [];
        figure('color' , 'white');
        
        for d = 1:5
            
            clear meanbnd allmeanbnd CB est_dist lab lab1
            ANA1 = getrow(ANA , ismember(ANA.Day , days{d}));
            
            X{d} = ANA1.estChnkBndry;
            lab = ANA1.seqNumb +1;
            
            for i = 1:size(X{d},1)
                id = ones(size(X{d} , 1) , 1);
                id(i) = 0;
                Y = X{d}(~id , :);
                X1 = X{d}(id==1 , :);
                lab1 = lab(id==1);
                for s = 1:3
                    m{d}(s,:) = nanmean(X1(lab1==s , :));
                    dis{d}(i , s) = pdist([Y;m{d}(s,:)], 'cityblock');%(m-Y)*pinv((m - nanmean(m))'*(Y - nanmean(Y)))*(m-Y)';%pdist2(m , Y , distance);
                end
            end
            
            seqs =  unique(lab);
            for l =  1:length(seqs)
                for l1 =  1:length(seqs)
                    id = lab == l;
                    D_CnkBndry{d}(l,l1) = nanmean(dis{d}(id , l1));
                    Dist_est_est(d).A{l,l1} = dis{d}(id , l1);
                    if l1~=l
                        offdiag_CnkBndry_es_es = [offdiag_CnkBndry_es_es  ; [ Dist_est_est(d).A{l,l1}, d*ones(length(Dist_est_est(d).A{l,l1}) , 1)]];
                    else
                        diag_CnkBndry_es_es = [diag_CnkBndry_es_es  ; [ Dist_est_est(d).A{l,l1}, d*ones(length(Dist_est_est(d).A{l,l1}) , 1)]];
                    end
                    
                end
            end
            % between estimated and actual
            for s = 1:3
                for s2 = 1:2
                    act_est_dist{d}(s,s2) = nanmean(nanmean(pdist2(CBD{s} , chbndry(s2 , :) , 'cityblock')));
                    Dist_est_act(d).A{s,s2} = pdist2(CBD{s} , chbndry(s2 , :) , 'cityblock');
                    if s~=s2+1
                        offdiag_CnkBndry_es_ac = [offdiag_CnkBndry_es_ac  ; [Dist_est_act(d).A{s,s2}, d*ones(length(Dist_est_act(d).A{s,s2}) , 1)]];
                    else
                        diag_CnkBndry_es_ac = [diag_CnkBndry_es_ac  ; [Dist_est_act(d).A{s,s2}, d*ones(length(Dist_est_act(d).A{s,s2}) , 1)]];
                    end
                end
                
            end
            
            %             offdiag_CnkBndry = [offdiag_CnkBndry  ; [reshape(D_CnkBndry{d}(2:end ,2:end)-diag(NaN*ones(length(D_CnkBndry{d}(2:end ,2:end)),1)) , (length(seqs)-1)^2 , 1) d*ones((length(seqs)-1)^2 , 1)]];
            %             diag_CnkBndry = [diag_CnkBndry  ; [diag(D_CnkBndry{d}(2:end ,2:end))  d*ones((length(seqs)-1) , 1)]];
            %
            
            
            figure('color' , [1 1 1])
            
            imagesc(D_CnkBndry{d});
            title(['Crossvalidated RDM for Est. Chunk boundries - day ' , num2str(days{d})] , 'FontSize' , 20)
            hold on
            ax = gca;
            ax.XTick = [1:3];
            ax.YTick = [1:3];
            ax.XTickLabel = {'Random' 'Structure 1' 'Structure 2'};
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2'};
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            hold on
            line([1.5 3.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            line([1.5 1.5] ,[1.5 3.5] ,  'LineWidth' , 3 , 'color' , [0 0 0])
            axis square
            colorbar
            
            figure('color' , [1 1 1])
            
            imagesc(act_est_dist{d});
            title(['Est. vs. Expected Chunking RDM - day ' , num2str(days{d})]  , 'FontSize' , 20)
            hold on
            ax = gca;
            line([.5 3.5] , [1.5 1.5] ,  'LineWidth' , 3 , 'color' , [0 0 0]);
            ax.XTick = [1:2];
            ax.YTick = [1:3];
            ax.XTickLabel = {'Structure 1' 'Structure 2' };
            ax.YTickLabel = {'Random' 'Structure 1' 'Structure 2' };
            xlabel('Expected');
            ylabel('Estimated');
            ax.XTickLabelRotation = 45;
            ax.FontSize = 20;
            axis square
            colorbar
        end
        
        
        figure('color' , 'white')
        imcounter = {[3:4] 1 2 };
        titl = {'Random'};
        for s = 1:size(ChnkPlcmnt , 1)
            titl = [titl num2str(ChnkPlcmnt(s , :))];
        end
        cnt = 1;
        for s = 1:3
            subplot(2,2,imcounter{cnt})
            for d = 1:5
                plot(m{d}(s,:), 'LineWidth' , 3)
                hold on
            end
            grid on
            set(gca , 'XLim' , [1 14] )
            if s>=2
                ids  = find(ChnkPlcmnt(s-1 , :) == 1);
                for j = 1:length(ids)
                    line([ids(j) ids(j)] , [0 1] , 'LineWidth' , 3 , 'color' , 'c' , 'LineStyle' , ':')
                end
            end
            legend({'Day1' 'Day2' , 'Day3' , 'Day4' , 'Day5'})
            title(['Average estimated Chunk boundries - structure ' , titl{s}])
            cnt = cnt+1;
        end
        
        
        h1 = figure;
        [CO_off , PL_off , ER_off] = lineplot(offdiag_CnkBndry_es_es(:,2) , offdiag_CnkBndry_es_es(:,1) , 'plotfcn' , 'nanmean');
        hold on
        [CO_on , PL_on , ER_on] = lineplot(diag_CnkBndry_es_es(:,2) , diag_CnkBndry_es_es(:,1) , 'plotfcn' , 'nanmean');
        close(h1)
        
        h1 = figure;
        [CO1_off , PL1_off , ER1_off] = lineplot(offdiag_CnkBndry_es_ac(:,2) , offdiag_CnkBndry_es_ac(:,1) , 'plotfcn' , 'nanmean');
        hold on
        [CO1_on , PL1_on , ER1_on] = lineplot(diag_CnkBndry_es_ac(:,2) , diag_CnkBndry_es_ac(:,1) , 'plotfcn' , 'nanmean');
        close(h1)
        
        figure('color' , 'white')
        
        subplot(2,1,1)
        h1 = plotshade(CO_off' , PL_off , ER_off,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
        hold on
        h2 = plotshade(CO_on' , PL_on , ER_on,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':');
        legend([h1 h2] ,{'Between structures' , 'Within structures'})
        grid on
        title('Average Dissimilarity Between and Within estimated Chunking Structures in City Block Chunk boundary distance')
        ax = gca;
        ax.XLim = [0 6];
        ax.XTick = [1:5];
        ax.XTickLabel = {'Day1' 'Day2' , 'Day3' , 'Day4' , 'Day5'};
        ax.FontSize = 20;
        ylabel('Average distance');
        
        subplot(2,1,2)
        h1 = plotshade(CO1_off' , PL1_off , ER1_off,'transp' , .2 , 'patchcolor' , 'b' , 'linecolor' , 'b' , 'linewidth' , 3 , 'linestyle' , ':');
        hold on
        h2 = plotshade(CO1_on' , PL1_on , ER1_on,'transp' , .2 , 'patchcolor' , 'r' , 'linecolor' , 'r' , 'linewidth' , 3 , 'linestyle' , ':');
        legend([h1 h2] ,{'Between structures' , 'Within structures'})
        grid on
        title('Average Dissimilarity Between and Within estimated/actual Chunking Structures in City Block Chunk boundary distance')
        ax = gca;
        ax.XLim = [0 6];
        ax.XTick = [1:5];
        ax.XTickLabel = {'Day1' 'Day2' , 'Day3' , 'Day4' , 'Day5'};
        ax.FontSize = 20;
        ylabel('Average distance');
        
        
        
        
        
        out = [];
end


