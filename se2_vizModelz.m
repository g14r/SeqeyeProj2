function sig = se2_vizModelz(what)

% what  is 's' for shadeplot and 'b' for barplot



% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
baseDir = '/Volumes/MotorControl/data/SeqEye2/analyze';  % server
baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye2/analyze';          %iMac


iN1 = input('Visualize Chunked or Random? (c/r)'  , 's');
iN2 = input('Ridge or OLS? (r/o)' , 's');
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
        switch iN2
            case 'o'
                load([baseDir , '/se2_CrossvalIPI_',titleSuffix,'-norm_OLS.mat'])
            case 'r'
                load([baseDir , '/se2_CrossvalIPI_',titleSuffix,'-norm_RR.mat']) 
        end
        
        
    case 'r'
        titleSuffix = 'Random';
        switch iN2
            case 'o'
                load([baseDir , '/se2_CrossvalIPI_',titleSuffix,'-norm_OLS.mat'])
            case 'r'
                load([baseDir , '/se2_CrossvalIPI_',titleSuffix,'-norm_RR.mat']) 
        end
        
end
%% set up legends and labels
clear xp_dev pp_dev ep_dev xp_r2 pp_r2 ep_r2 xp_r2a pp_r2a ep_r2a dev_image R2_image xp_dn  pp_dn ep_dn dvn_image xp_cor  pp_cor  ep_cor sig
horzSize = {1,2,3,4,5,6,7,8,13,[5:13]}; 
legenddHor = {'H = 1', 'H = 2','H = 3','H = 4','H = 5','H = 6','H = 7','H = 8','H = 13' , 'H = 5-13'};
% x labels for shadeplot
label = {'R'  'C+R',     '1st+R' ,'1st+2nd+R'    '1st+2nd+3rd+R' , 'C+1st+R'  ,  'C+1st+2nd+R'  ,  'Full'};
dayz = {[1] [2 3] [4 5]};
plotIND = [2:5 , 8]; % the model indices to include in bar plot
% do significance tests:

% 'right' tail testing that A > B
% 'left' tail testing that A < B

count = 1;
for h = 1:length(horzSize)
    for dd = 1:length(dayz)
        ido = (Mdl.hor == h & ismember(Mdl.day , dayz{dd}));
        A = Mdl.R2(Mdl.xx == 2 & ido);
        B = Mdl.R2(Mdl.xx == 5 & ido);
        C = Mdl.R2(Mdl.xx == 8 & ido);
        [~ , sig.chunkVprob(count , 1)] = ttest2(A , B);
        [~ , sig.chunkVfull(count , 1)] = ttest2(A , C);
        [~ , sig.probVfull(count , 1)] = ttest2(B , C);
        sig.day(count , 1) = dd;
        sig.h(count , 1) = h;
        count  = count +1;
    end
end


% Model descriptions
cleanLabel = {'within/between Chunk', '1st order probability' ,'1st + 2nd order probability' ,'1st + 2nd + 3rd order probability' ,'Full Model'};
cat_cleanLabel = categorical(repmat(cleanLabel , length(horzSize)-1 , 1));

% =================== % =================== % =================== % =================== summarize

K = tapply(Mdl , {'hor' , 'subj' , 'day' , 'xx'} , {'R2' , 'nanmean'} , {'corYY' , 'nanmean'} ); % average over cv loops


h1 = figure;
% I find color setting and legending  in lineplot a pain, so I generate the data using lineplot, close that figure and pop them into shadeplot and barplot
count = 1;
for h = 1:length(horzSize)
    for dd = 1:length(dayz)
        ido = (K.hor == h & ismember(K.day , dayz{dd}));
        [xp_cor{dd}(h, :) , pp_cor{dd}(h,:) , ep_cor{dd}(h,:)]  = lineplot(K.xx, K.corYY , 'plotfcn','nanmean' ,  'subset' , ido);
        hold on
        [xp_r2{dd}(h, :) , pp_r2{dd}(h,:) , ep_r2{dd}(h,:)]  = lineplot(K.xx, K.R2 , 'plotfcn','nanmean' , 'subset', ido);
        hold on
        A = K.R2(K.xx == 2 & ido);
        B = K.R2(K.xx == 5 & ido);
        C = K.R2(K.xx == 8 & ido);
        [~ , sig.chunkVprob(count , 1)] = ttest2(A , B);
        [~ , sig.chunkVfull(count , 1)] = ttest2(A , C);
        [~ , sig.probVfull(count , 1)] = ttest2(B , C);
        sig.day(count , 1) = dd;
        sig.h(count , 1) = h;
        count  = count +1;
    end
end
close(h1)


%% do the plotting
switch what
    case 's'
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
        
        
        
    case 'b'
        % *******************barplot
        % =================== % =================== % =================== % =================== Visualize model R2 comparisons
se2_linearIPIModels.m        xp_ = reshape(cell2mat(xp_r2) , size(xp_r2{1} , 1) , size(xp_r2{1} , 2) , length(xp_r2));
        pp_ = reshape(cell2mat(pp_r2) , size(pp_r2{1} , 1) , size(pp_r2{1} , 2) , length(pp_r2));
        ep_ = reshape(cell2mat(ep_r2) , size(ep_r2{1} , 1) , size(ep_r2{1} , 2) , length(ep_r2));
        

        figure('color' , 'white')
        for i = 1:length(dayz)
            subplot(length(dayz) ,1, i)
            bar(squeeze(pp_(1:9,plotIND , i)));
            grid on
            set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(legenddHor)-1] , 'XTickLabel' , legenddHor(1:end-1),'YLim' , ylim)
            title([titleSuffix ,' - Model Crossvalidated R^2 on Day(s) ' , num2str(dayz{i})])
        end
        legend(cleanLabel)
        % =================== % =================== % =================== % =================== Visualize model correlation comparisons
        xp_ = reshape(cell2mat(xp_cor) , size(xp_cor{1} , 1) , size(xp_cor{1} , 2) , length(xp_cor));
        pp_ = reshape(cell2mat(pp_cor) , size(pp_cor{1} , 1) , size(pp_cor{1} , 2) , length(pp_cor));
        ep_ = reshape(cell2mat(ep_cor) , size(ep_cor{1} , 1) , size(ep_cor{1} , 2) , length(ep_cor));
        
        
        figure('color' , 'white')
        for i = 1:length(dayz)
            subplot(length(dayz) ,1, i)
            bar(squeeze(pp_(1:9,plotIND , i)));
            grid on
            set(gca , 'FontSize' , 20 ,'Box' , 'off' , 'GridAlpha' , 1 , 'XTick' , [1:length(legenddHor)-1] , 'XTickLabel' , legenddHor(1:end-1),...
                'YLim' , [-.1 .4])
            title([titleSuffix , ' - Model Prediction - Output Correlation on Day(s) ' , num2str(dayz{i})])
        end
        legend(cleanLabel)
end


