function sig = se2_vizModelz(what)

% what  is 's' for shadeplot and 'b' for barplot



% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
baseDir = '/Volumes/MotorControl/data/SeqEye2/analyze';  % server
%baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';          %iMac


iN1 = input('Visualize Chunked or Random? (c/r)'  , 's');
iN3 = input('Fit or Crossvalidated? (f/c)' , 's');
switch iN3
    case 'c'
        iN2 = input('Ridge or OLS? (r/o)' , 's');
end

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

%% load proper model

switch iN1
    case 'c'
        titleSuffix = 'Chunked';
    case 'r'
        titleSuffix = 'Random';
end
switch iN3
    case 'c'
        switch iN2
            case 'o'
                load([baseDir , '/se2_CrossvalIPI_',titleSuffix,'-norm_OLS.mat'])
            case 'r'
                load([baseDir , '/se2_CrossvalIPI_',titleSuffix,'-norm_RR.mat'])
        end
    case 'f'
        load([baseDir , '/se2_fitIPI_',titleSuffix,'-norm_OLS.mat'])
end


%% set up legends and labels
clear xp_dev pp_dev ep_dev xp_r2 pp_r2 ep_r2 xp_r2a pp_r2a ep_r2a dev_image R2_image xp_dn  pp_dn ep_dn dvn_image xp_cor  pp_cor  ep_cor sig
horzSize = {1,2,3,[4:13]}; 
legenddHor = {'H = 1', 'H = 2','H = 3','H = 4','H = 5','H = 6','H = 7','H = 8','H = 13' , 'H = 5-13'};
% x labels for shadeplot
label = {'R'  'C+R',     '1st+R' ,'1st+2nd+R'    '1st+2nd+3rd+R' , 'C+1st+R'  ,  'C+1st+2nd+R'  ,  'Full'};
dayz = {[2 3] [4] [5]};
plotIND = [2 5 , 8]; % the model indices to include in bar plot
% do significance tests:

% 'right' tail testing that A > B
% 'left' tail testing that A < B
subj = unique(Mdl.SN);
count = 1;
for h = 1:length(horzSize)
    for dd = 1:length(dayz)
        for sb = 1:length(subj)
            ido = (Mdl.Horizon == h & ismember(Mdl.Day , dayz{dd}) & ismember(Mdl.SN , subj(sb)));
            A = getrow(Mdl , Mdl.modelNum == 2 & ido);
            B = getrow(Mdl , Mdl.modelNum == 5 & ido);
            C = getrow(Mdl , Mdl.modelNum == 8 & ido);
            N = getrow(Mdl , Mdl.modelNum == 1 & ido);
            
            [~ , sig.R2_chunkVprob(count , 1)] = ttest2(A.R2 , B.R2);
            [~ , sig.R2_chunkVfull(count , 1)] = ttest2(A.R2 , C.R2);
            [~ , sig.R2_probVfull(count , 1)] = ttest2(B.R2 , C.R2);
            [~ , sig.R2_NVprob(count , 1)] = ttest2(N.R2 , B.R2);
            [~ , sig.R2_NVchunk(count , 1)] = ttest2(A.R2 , N.R2);
            [~ , sig.R2_NVfull(count , 1)] = ttest2(N.R2 , C.R2);
            
            [~ , sig.rel_R_chunkVprob(count , 1)] = ttest2(A.rel_R , B.rel_R);
            [~ , sig.rel_R_chunkVfull(count , 1)] = ttest2(A.rel_R , C.rel_R);
            [~ , sig.rel_R_probVfull(count , 1)] = ttest2(B.rel_R , C.rel_R);
            [~ , sig.rel_R_NVprob(count , 1)] = ttest2(N.rel_R , B.rel_R);
            [~ , sig.rel_R_NVchunk(count , 1)] = ttest2(A.rel_R , N.rel_R);
            [~ , sig.rel_R_NVfull(count , 1)] = ttest2(N.rel_R , C.rel_R);
            
            
            [~ , sig.R_chunkVprob(count , 1)] = ttest2(A.R , B.R);
            [~ , sig.R_chunkVfull(count , 1)] = ttest2(A.R , C.R);
            [~ , sig.R_probVfull(count , 1)] = ttest2(B.R , C.R);
            [~ , sig.R_NVprob(count , 1)] = ttest2(N.R , B.R);
            [~ , sig.R_NVchunk(count , 1)] = ttest2(A.R , N.R);
            [~ , sig.R_NVfull(count , 1)] = ttest2(N.R , C.R);
            
            [~ , sig.aic_chunkVprob(count , 1)] = ttest2(A.rel_AIC , B.rel_AIC);
            [~ , sig.aic_chunkVfull(count , 1)] = ttest2(A.rel_AIC , C.rel_AIC);
            [~ , sig.aic_probVfull(count , 1)] = ttest2(B.rel_AIC , C.rel_AIC);
            [~ , sig.aic_NVprob(count , 1)] = ttest2(N.rel_AIC , B.rel_AIC);
            [~ , sig.aic_NVchunk(count , 1)] = ttest2(A.rel_AIC , N.rel_AIC);
            [~ , sig.aic_NVfull(count , 1)] = ttest2(N.rel_AIC , C.rel_AIC);
            
            
            sig.day(count , 1) = dd;
            sig.h(count , 1) = h;
            sig.sn(count , 1) = sb;
            count  = count +1;
        end
    end
end

S = tapply(sig , {'h' , 'day'} , {'rel_R_chunkVprob' , 'nanmean'} , {'rel_R_chunkVfull' , 'nanmean'} ); % average over cv loops
S = getrow(sig , sig.h == 1 & sig.day ==3);


% Model descriptions
% cleanLabel = {'within/between Chunk', '1st order probability' ,'1st + 2nd order probability' ,'1st + 2nd + 3rd order probability' ,'Full Model'};
cleanLabel = {'within/between Chunk' ,'1st + 2nd + 3rd order probability' ,'Full Model'};

cat_cleanLabel = categorical(repmat(cleanLabel , length(horzSize)-1 , 1));

% =================== % =================== % =================== % =================== summarize

K = tapply(Mdl , {'Horizon' , 'SN' , 'Day' , 'modelNum'} , {'R2' , 'nanmean'} , {'rel_R' , 'nanmean'}, {'rel_AIC' , 'nanmean'},{'R' , 'nanmean'} ); % average over cv loops


h1 = figure;
% I find color setting and legending  in lineplot a pain, so I generate the data using lineplot, close that figure and pop them into shadeplot and barplot
for h = 1:length(horzSize)
    for dd = 1:length(dayz)
        ido = (K.Horizon == h & ismember(K.Day , dayz{dd}));
        [xp_cor_rel{dd}(h, :) , pp_cor_rel{dd}(h,:) , ep_cor_rel{dd}(h,:)]  = lineplot(K.modelNum, K.rel_R , 'plotfcn','nanmean' ,  'subset' , ido);
        hold on
        [xp_r2{dd}(h, :) , pp_r2{dd}(h,:) , ep_r2{dd}(h,:)]  = lineplot(K.modelNum, K.R2 , 'plotfcn','nanmean' , 'subset', ido);
        [xp_aic{dd}(h, :) , pp_aic{dd}(h,:) , ep_aic{dd}(h,:)]  = lineplot(K.modelNum, K.rel_AIC , 'plotfcn','nanmean' , 'subset', ido);
        [xp_cor{dd}(h, :) , pp_cor{dd}(h,:) , ep_cor{dd}(h,:)]  = lineplot(K.modelNum, K.R , 'plotfcn','nanmean' ,  'subset' , ido);
        hold on
    end
end
close(h1)


%% do the plotting
switch what
    
    case 's'
        lo = min(min(min(pp_r2 + ep_r2)));
        hi = max(max(max(pp_r2 - ep_r2)));
        ylim = [1.1*lo   1.1*hi];
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons
        figure('color' , 'white')
        
        for dd = 1:length(dayz)
            subplot(2,3,dd)
            for h = 1:length(horzSize)-1
                hold on
                eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_r2{dd}(h,plotIND) , ep_r2{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                hold on
            end
            ylabel('R^2')
            set(gca ,'YLim' , ylim, 'XLim' , [0 length(plotIND)+1] ,'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' R^2, Days ' , num2str(dayz{dd})])
             
            grid on
        end
        legend([h1,h2,h3,h4,h5,h6,h7,h8,h9] ,legenddHor(1:end-1), 'Box' , 'off')
        
        
        for dd = 1:length(dayz)
            subplot(2,3,dd+length(dayz))
            for h = length(horzSize)
                hold on
                eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_r2{dd}(h,plotIND) , ep_r2{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                hold on
            end
            ylabel('R^2')
            set(gca ,'YLim' , ylim, 'XLim' , [0 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , cleanLabel , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' R^2, Days ' , num2str(dayz{dd})])
            grid on
        end
        legend([h10] ,legenddHor(end), 'Box' , 'off')
        
        % =================== % =================== % =================== % =================== Visualize model correlation comparisons
        lo = min(min(min(pp_cor + ep_cor)));
        hi = max(max(max(pp_cor - ep_cor)));
        ylim = [1.1*lo   1.1*hi];
        
        figure('color' , 'white')
        
        for dd = 1:length(dayz)
            subplot(2,3,dd)
            for h = 1:length(horzSize)-1
                hold on
                eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_cor{dd}(h,plotIND) , ep_cor{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                hold on
            end
            ylabel('correlation')
            set(gca , 'XLim' , [1 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , label(1:end) , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' Model Correlations, Days ' , num2str(dayz{dd})])
            
            grid on
        end
        legend([h1,h2,h3,h4,h5,h6,h7,h8,h9] ,legenddHor(1:end-1), 'Box' , 'off')
        
        
        for dd = 1:length(dayz)
            subplot(2,3,dd+length(dayz))
            for h = length(horzSize)
                hold on
                eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_cor{dd}(h,plotIND) , ep_cor{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                hold on
            end
            ylabel('correlation')
            set(gca , 'XLim' , [1 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , label(1:end) , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' Model Correlations, Days ' , num2str(dayz{dd})])
            grid on
        end
        legend([h10] ,legenddHor(end), 'Box' , 'off')
        
        % =================== % =================== % =================== % =================== Visualize model AIC comparisons
        lo = min(min(min(pp_aic + ep_aic)));
        hi = max(max(max(pp_aic - ep_aic)));
        ylim = [1.1*lo   1.1*hi];
        
        figure('color' , 'white')
        
        for dd = 1:length(dayz)
            subplot(2,3,dd)
            for h = 1:length(horzSize) - 1
                hold on
                eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_aic{dd}(h,plotIND) , ep_aic{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                hold on
            end
            set(gca , 'XLim' , [1 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , label(1:end) , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' Model Relative AIC, Days ' , num2str(dayz{dd})])
            
            grid on
        end
        legend([h1,h2,h3] ,legenddHor(1:end-1), 'Box' , 'off')
        
        
        for dd = 1:length(dayz)
            subplot(2,3,dd+length(dayz))
            for h = length(horzSize)
                hold on
                eval(['h' , num2str(h) , ' = plotshade([1:length(plotIND)] , pp_cor{dd}(h,plotIND) , ep_cor{dd}(h,plotIND),''transp'' , .2 , ''patchcolor'' , colors(h,:) , ''linecolor'' , colors(h,:) , ''linewidth'' , 3 , ''linestyle'' , '':'')']);
                hold on
            end
            set(gca , 'XLim' , [1 length(plotIND)+1],'XTick' , [1: length(plotIND)] , 'XTickLabels' , label(1:end) , 'FontSize' , 20 ,...
                'XTickLabelRotation',45,'Box' , 'off' , 'GridAlpha' , 1)
            title([titleSuffix , ' Model Relative AIC, Days ' , num2str(dayz{dd})])
            grid on
        end
        legend([h4] ,legenddHor(end), 'Box' , 'off')
        
    case 'b'
        % *******************barplot
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons
        xp_ = reshape(cell2mat(xp_r2) , size(xp_r2{1} , 1) , size(xp_r2{1} , 2) , length(xp_r2));
        pp_ = reshape(cell2mat(pp_r2) , size(pp_r2{1} , 1) , size(pp_r2{1} , 2) , length(pp_r2));
        ep_ = reshape(cell2mat(ep_r2) , size(ep_r2{1} , 1) , size(ep_r2{1} , 2) , length(ep_r2));
        
        
        
        lo = min(min(min(pp_)));
        hi = max(max(max(pp_)));
        ylim = [1.1*lo   1.1*hi];
        figure('color' , 'white')
        hlabel = {'1' ,'2' , '3' ,'4 - 13 (Full)'};
        fc = 1;
        for i = [1 2 3 4]
            subplot(4 ,1, i)
%             bar(squeeze(pp_(1:9,plotIND , i)));
%             set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(legenddHor)-1] ,...
%                 'YLim' , ylim,'XTickLabel' , legenddHor(1:end-1))
            bar(squeeze(pp_(i,plotIND , :))');
            set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] ,...
                'YLim' , ylim,'XTickLabel' , {'Day 1' , 'Days 2, 3' , 'Days 4, 5'})
            grid on
            title([titleSuffix , ' - Model Relative R^2 on Horizon Size(s) ' , hlabel{fc}]);
            fc = fc+1;
        end
        legend(cleanLabel)
        % =================== % =================== % =================== % =================== Visualize model correlation comparisons
        xp_ = reshape(cell2mat(xp_cor) , size(xp_cor{1} , 1) , size(xp_cor{1} , 2) , length(xp_cor));
        pp_ = reshape(cell2mat(pp_cor) , size(pp_cor{1} , 1) , size(pp_cor{1} , 2) , length(pp_cor));
        ep_ = reshape(cell2mat(ep_cor) , size(ep_cor{1} , 1) , size(ep_cor{1} , 2) , length(ep_cor));
        
       
        
        lo = min(min(min(pp_ - ep_)));
        hi = max(max(max(pp_ + ep_)));
        ylim = [1.1*lo   1.1*hi];
        figure('color' , 'white')
        hlabel = {'1' ,'2' , '3' ,'4 - 13 (Full)'};
        fc = 1;
        for i = [1 2 3 4]
            subplot(4 ,1, i)
%             bar(squeeze(pp_(1:9,plotIND , i)));
%             set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(legenddHor)-1] ,...
%                 'YLim' , ylim,'XTickLabel' , legenddHor(1:end-1))
%             bar(squeeze(-pp_(i,plotIND , :))');            
            barwitherr(squeeze(ep_(i,plotIND ,:)) , squeeze(pp_(i,plotIND ,:)))
            set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] ,...
                'YLim' , ylim,'XTickLabel' , {'Day 2' , 'Day 3' , 'Day 4'})
            grid on
            title([titleSuffix , ' - Model Prediction - Output Correlation on Horizon Size(s) ' , hlabel{fc}]);
            fc = fc+1;
        end
        legend(cleanLabel)
        
        % =================== % =================== % =================== %
        % =================== Visualize model AIC
        xp_ = reshape(cell2mat(xp_aic) , size(xp_aic{1} , 1) , size(xp_aic{1} , 2) , length(xp_aic));
        pp_ = reshape(cell2mat(pp_aic) , size(pp_aic{1} , 1) , size(pp_aic{1} , 2) , length(pp_aic));
        ep_ = reshape(cell2mat(ep_aic) , size(ep_aic{1} , 1) , size(ep_aic{1} , 2) , length(ep_aic));
        
        lo = min(min(min(-pp_ - ep_)));
        hi = max(max(max(-pp_ + ep_)));
        ylim = [1.1*lo   1.1*hi];
        figure('color' , 'white')
        hlabel = {'1' ,'2' , '3' ,'4 - 13 (Full)'};
        fc = 1;
        for i = [1 2 3 4]
            subplot(4 ,1, i)
%             bar(squeeze(pp_(1:9,plotIND , i)));
%             set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(legenddHor)-1] ,...
%                 'YLim' , ylim,'XTickLabel' , legenddHor(1:end-1))
%             bar(squeeze(pp_(i,plotIND , :))');
            barwitherr(squeeze(ep_(i,plotIND ,:)) , -squeeze(pp_(i,plotIND ,:)))
            set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(dayz)] ,...
                'YLim' , ylim,'XTickLabel' , {'Day 1' , 'Days 2, 3' , 'Days 4, 5'})
            grid on
            title([titleSuffix , ' - Model Prediction - Relative AIC on Horizon Size(s) ' , hlabel{fc}]);
            fc = fc+1;
        end
        legend(cleanLabel)
        
end


