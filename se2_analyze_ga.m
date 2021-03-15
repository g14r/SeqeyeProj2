function [varargout] = se2_analyze_ga(what, varargin)
%% [varargout] = se2_analyze_ga(what, varargin)
%
%
%       [D] = se2_analyze_ga('horizon_MT');
%       [D] = se2_analyze_ga('horizon_IPI');
%       [D] = se2_analyze_ga('learning_RT_MT');
%       [D] = se2_analyze_ga('learning_IPI');
%       [D] = se2_analyze_ga('RTIPI_scatter');
%
% --
% gariani@uwo.ca - 2020.02.13

%% paths
path_to_data = '/Volumes/MotorControl/data/SeqEye2/analyze';

%% global vars
% subjects
subj_name = {'s01_2', 's02_2', 's03_2', 's04_2', 's05_2', 's06_2', 's07_2', 's08_2', 's09_2', 's10_2', 's11_2', 's12_2', 's13_2'};
subjnum = 1:length(subj_name); % all subjects

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray2 = [100,100,100]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
white = [255,255,255]/255;

% plot defaults
fs = 24; %default font size for all figures
lw = 3; %default line width for all figures
ms = 8; %default marker size for all figures

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,gray2,lightgray,green,lightgreen,black,silver,white,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
daysty = style.custom({lightgray, gray, gray2, darkgray, black}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
horsty = style.custom({black, darkgray, gray2, gray, lightgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
lrnsty = style.custom({lightgray, darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
lightsty = style.custom({lightgray}, 'markertype','none', 'linewidth',1, 'errorbars','none');
dgsty = style.custom({darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
bksty = style.custom({black}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
grsty = style.custom({gray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');

% legends
dayleg = {'Day 1', 'Day 2', 'Day 3', 'Day 4', 'Day 5'};
horleg = {'W = 1', 'W = 2', 'W = 3', 'W = 4', 'W = 5+'};
lrnleg = {'Day 1 | 2', 'Day 4 | 5'};

%% types of analysis
switch (what)
    case 'horizon_MT' % MT vs horizon
        % load data
        Dall = load( fullfile(path_to_data, 'se2_ds_N=13.mat') );
        Dall = Dall.Dall3;
        Dall.RT = Dall.AllPressTimes(:,1)-1500;
        ANA = getrow(Dall, ismember(Dall.SN, subjnum) & Dall.isgood==1 & Dall.isError==0 & ismember(Dall.seqNumb, 0) & Dall.RT <= 9000);
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500, ANA.IPI];
        
        % pool days
        ANA.Day(ANA.Day==2) = 1;
        ANA.Day(ANA.Day==4) = 5;
        
        % open multi-panel figure
        figure('Name', 'Horizon MT'); set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        % ------------------------------------------------------------------------------------------------------
        % summarize data
        T = tapply(ANA, {'SN', 'Horizon'}, ...
            {ANA.MT, 'nanmedian', 'name', 'MT'});
        T = normData(T, {'MT'});
        % plot data
        subplot(2,2,[1,3]);
        plt.scatter(T.Horizon, T.normMT/1000, 'regression','none', 'style',dgsty, 'markersize',1);
        hold on;
        plt.line(T.Horizon, T.normMT/1000, 'plotfcn','nanmean', 'style',bksty);
        xlabel('Viewing window (W)'); ylabel('Movement time (s)'); set(gca,'fontsize',fs);
        %stats
        T.ANOVA = anovaMixed(T.MT, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 1:13));
        
        % ------------------------------------------------------------------------------------------------------
        % summarize data
        T = tapply(ANA, {'SN', 'Horizon', 'Day'}, ...
            {ANA.MT, 'nanmedian', 'name', 'MT'});
        T = normData(T, {'MT'});
        % plot data
        subplot(222);
        plt.scatter(T.Horizon, T.normMT/1000, 'split',T.Day, 'style',lrnsty, 'leg',lrnleg, 'regression','none', 'subset',ismember(T.Day,[1,2, 4,5]), 'markersize',1);
        hold on;
        plt.line(T.Horizon, T.normMT/1000, 'split',T.Day, 'style',lrnsty, 'leg',lrnleg, 'plotfcn','nanmean', 'subset',ismember(T.Day,[1,2, 4,5]));
        ylabel('MT (s)'); set(gca,'fontsize',fs); ylim([2.8 8.2]); axis square;
        
        % ------------------------------------------------------------------------------------------------------
        % summarize data
        T = tapply(ANA, {'SN', 'Horizon'}, ...
            {ANA.MT, 'nanmedian', 'name', 'MT'}, ...
            {ANA.MT, 'nanmedian', 'name', 'MTd1', 'subset',ismember(ANA.Day,1)}, ...
            {ANA.MT, 'nanmedian', 'name', 'MTd5', 'subset',ismember(ANA.Day,5)});
        T = normData(T, {'MT', 'MTd1', 'MTd5'});
        % plot data
        subplot(224);
        plt.scatter(T.Horizon, (nanplus(T.normMTd1,-T.normMTd5)./T.normMT)*100, 'regression','none', 'style',dgsty, 'markersize',1);
        hold on;
        plt.line(T.Horizon, (nanplus(T.normMTd1,-T.normMTd5)./T.normMT)*100, 'plotfcn','nanmean', 'style',bksty);
        xlabel('Viewing window (W)'); ylabel('MT difference (% of avg MT)'); set(gca,'fontsize',fs); axis square;
        drawline(0, 'dir','horz', 'linestyle','--');
        %stats
        % day 1-2
        T.ANOVA = anovaMixed(T.MTd1, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 1:13));
        T.ANOVA = anovaMixed(T.MTd1, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 2:13));
        T.ANOVA = anovaMixed(T.MTd1, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 3:13));
        % day 4-5
        T.ANOVA = anovaMixed(T.MTd5, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 1:13));
        T.ANOVA = anovaMixed(T.MTd5, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 2:13));
        T.ANOVA = anovaMixed(T.MTd5, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 3:13));
        T.ANOVA = anovaMixed(T.MTd5, T.SN,'within', [T.Horizon], {'Horizon'}, 'subset',ismember(T.Horizon, 4:13));
        
        % out
        varargout = {T}; %return main structure
        
    case 'horizon_IPI' % IPI vs horizon
        % load data
        Dall = load( fullfile(path_to_data, 'se12_ds_N=25.mat') );
        Dall = Dall.D;
        ANA = getrow(Dall, ismember(Dall.SN, subjnum) & Dall.isgood==1 & Dall.isError==0 & ismember(Dall.seqNumb, 0));
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500, ANA.IPI];
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        for tn  = 1:length(ANA.TN)
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,14);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,14);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,14);
            ANA.IPI_prsnumb(tn , :) = 0:13;
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,14);
            ANA.IPI_BN(tn , :) = ANA.BN(tn)*ones(1,14);
        end
        for tn  = 1:length(ANA.TN)
            ANA.badpress{tn,1} = 0;
            for prs=1:13
                ANA.badpress{tn,prs+1} = double(sum(ANA.badPress(tn,prs:prs+1))>0);
            end
        end
        
        IPItable.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPItable.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        IPItable.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        IPItable.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        IPItable.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        IPItable.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        IPItable.BN = reshape(ANA.IPI_BN , numel(ANA.IPI) , 1);
        IPItable.badpress = cell2mat(reshape(ANA.badpress , numel(ANA.IPI) , 1));
        
        %%%%
        IPItable.Horizon(IPItable.Horizon >= 5) = 5;
        %%%%
        
        % summarize data
        IPI = getrow(IPItable , ~IPItable.badpress);
        IPI  = tapply(IPI , {'Horizon', 'SN', 'prsnumb', 'Day'} , {'IPI', 'nanmedian(x)'});
        IPI = normData(IPI , {'IPI'});
        
        % plot data
        figure('Name', 'Horizon IPI'); set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        subplot(2,2,3:4);
        plt.line([IPI.Horizon IPI.prsnumb], IPI.normIPI, 'style',dgsty, 'subset',IPI.prsnumb>0 & ismember(IPI.Day,[1,2,3,4,5]));
        xlabel('Transition number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); %axis square;
        xtl = repmat([{'1'}, repmat({''},1,numel(unique(IPI.prsnumb))-2), {'13'}],1,numel(unique(IPI.Horizon))); xticklabels(xtl);
        %         figure('Name', 'IPI'); set(gcf, 'Renderer','painters');
        %         plt.line([IPI.prsnumb], IPI.IPI, 'split',[IPI.Horizon], 'style',horsty, 'leg',horleg, 'leglocation','east', 'subset',IPI.prsnumb>0 & IPI.prsnumb<11 & ismember(IPI.Day, 1:5));
        %         xlabel('Transition number'); ylabel('IPI (ms)'); set(gca,'fontsize',fs); ylim([225 575]); axis square;
        
        % out
        varargout = {IPI}; %return main structure
        
    case 'learning_RT_MT' % learning effects on RT and MT
        % load data
        Dall = load( fullfile(path_to_data, 'se12_ds_N=25.mat') );
        Dall = Dall.D;
        Dall.RT = Dall.AllPressTimes(:,1)-1500;
        ANA = getrow(Dall, ismember(Dall.SN, subjnum) & Dall.isgood==1 & Dall.isError==0 & ismember(Dall.seqNumb, 0) & Dall.RT <= 9000);
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500, ANA.IPI];
        
        % summarize data
        T = tapply(ANA, {'SN','Horizon','Day'},...
            {'RT', 'nanmedian(x)'},...
            {'MT', 'nanmedian(x)'});
        T = normData(T, {'MT', 'RT'});
        
        % plot data
        figure('Name', 'Learning RT MT'); set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        subplot(221);
        plt.line([T.Horizon], T.normMT, 'split',T.Day, 'style',trsty, 'leg',trleg, 'subset',ismember(T.Day,[1,5]));
        xlabel('Horizon (W)'); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); axis square;
        %         plt.line([T.Day T.Horizon], T.normRT, 'split',T.Day, 'style',daysty, 'leg',dayleg, 'leglocation','NorthWest');
        %         xlabel('Horizon (W)'); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); %axis square;
        %         ylim([500 900]);
        %         xtl = repmat([{'1'}, repmat({''},1,numel(unique(T.Horizon))-2), {'13'}],1,numel(unique(T.Day))); xticklabels(xtl);
        %         subplot(223);
        %         plt.line([T.Day T.Horizon], T.normMT, 'split',T.Day, 'style',daysty, 'leg',dayleg);
        %         xlabel('Horizon (W)'); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); %axis square;
        %         ylim([3000 8000]);
        %         xtl = repmat([{'1'}, repmat({''},1,numel(unique(T.Horizon))-2), {'13'}],1,numel(unique(T.Day))); xticklabels(xtl);
        %%%%
        % T.Horizon(T.Horizon >= 5) = 5;
        %%%%
        subplot(222);
        plt.line([T.Horizon], T.normRT, 'split',T.Day, 'style',trsty, 'leg',trleg, 'subset',ismember(T.Day,[1,5]));
        xlabel('Horizon (W)'); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
        ylim([500 900]);
        %         plt.line([T.Horizon T.Day], T.normRT, 'split',T.Horizon, 'style',horsty, 'leg',horleg, 'leglocation','NorthWest');
        %         xlabel('Day'); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); %axis square;
        %         ylim([500 900]);
        %         subplot(224);
        %         plt.line([T.Horizon T.Day], T.normMT, 'split',T.Horizon, 'style',horsty, 'leg',horleg);
        %         xlabel('Day'); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); %axis square;
        %         ylim([3000 8000]);
        
        % out
        varargout = {T}; %return main structure
        
    case 'learning_IPI' % IPI vs horizon
        % load data
        Dall = load( fullfile(path_to_data, 'se12_ds_N=25.mat') );
        Dall = Dall.D;
        ANA = getrow(Dall, ismember(Dall.SN, subjnum) & Dall.isgood==1 & Dall.isError==0 & ismember(Dall.seqNumb, 0));
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500, ANA.IPI];
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        for tn  = 1:length(ANA.TN)
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,14);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,14);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,14);
            ANA.IPI_prsnumb(tn , :) = 0:13;
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,14);
            ANA.IPI_BN(tn , :) = ANA.BN(tn)*ones(1,14);
        end
        for tn  = 1:length(ANA.TN)
            ANA.badpress{tn,1} = 0;
            for prs=1:13
                ANA.badpress{tn,prs+1} = double(sum(ANA.badPress(tn,prs:prs+1))>0);
            end
        end
        
        IPItable.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPItable.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        IPItable.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        IPItable.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        IPItable.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        IPItable.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        IPItable.BN = reshape(ANA.IPI_BN , numel(ANA.IPI) , 1);
        IPItable.badpress = cell2mat(reshape(ANA.badpress , numel(ANA.IPI) , 1));
        
        %%%%
        IPItable.Horizon(IPItable.Horizon >= 5) = 5;
        %%%%
        
        % summarize data
        IPI = getrow(IPItable , ~IPItable.badpress);
        IPI  = tapply(IPI , {'Horizon', 'SN', 'prsnumb', 'Day'} , {'IPI', 'nanmedian(x)'});
        IPI = normData(IPI , {'IPI'});
        
        % plot data
        figure('Name', 'Learning IPI'); set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        subplot(2,2,3:4);
        plt.line([IPI.Horizon IPI.prsnumb], IPI.normIPI, 'split',IPI.Day, 'style',lrnsty, 'leg',lrnleg, 'subset',IPI.prsnumb>0 & ismember(IPI.Day,[1,5]));
        xlabel('Transition number'); ylabel('Inter-press interval (ms)'); set(gca,'fontsize',fs); %axis square;
        xtl = repmat([{'1'}, repmat({''},1,numel(unique(IPI.prsnumb))-2), {'13'}],1,numel(unique(IPI.Horizon))); xticklabels(xtl);
        
        % out
        varargout = {IPI}; %return main structure
        
    case 'RTIPI_scatter' % scatterplot of RT vs IPI
        % load data
        Dall = load( fullfile(path_to_data, 'se12_ds_N=25.mat') );
        Dall = Dall.D;
        ANA = getrow(Dall, ismember(Dall.SN, subjnum) & Dall.isgood==1 & Dall.isError==0 & ismember(Dall.seqNumb, 0));
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500, ANA.IPI];
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        for tn  = 1:length(ANA.TN)
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,14);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,14);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,14);
            ANA.IPI_prsnumb(tn , :) = 0:13;
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,14);
            ANA.IPI_BN(tn , :) = ANA.BN(tn)*ones(1,14);
        end
        for tn  = 1:length(ANA.TN)
            ANA.badpress{tn,1} = 0;
            for prs=1:13
                ANA.badpress{tn,prs+1} = double(sum(ANA.badPress(tn,prs:prs+1))>0);
            end
        end
        
        IPItable.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPItable.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        IPItable.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        IPItable.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        IPItable.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        IPItable.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        IPItable.BN = reshape(ANA.IPI_BN , numel(ANA.IPI) , 1);
        IPItable.badpress = cell2mat(reshape(ANA.badpress , numel(ANA.IPI) , 1));
        
        %%%%
        % load data
        Dall = load( fullfile(path_to_data, 'se12_ds_N=25.mat') );
        Dall = Dall.D;
        Dall.RT = Dall.AllPressTimes(:,1)-1500;
        ANA = getrow(Dall, ismember(Dall.SN, subjnum) & Dall.isgood==1 & Dall.isError==0 & ismember(Dall.seqNumb, 0));% & Dall.RT <= 9000);
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500, ANA.IPI];
        
        % summarize data
        T = tapply(ANA, {'SN','Horizon','Day'},...
            {'RT', 'nanmedian(x)'});
        
        % summarize data
        IPI = getrow(IPItable , ~IPItable.badpress);
        IPI  = tapply(IPI , {'Horizon', 'SN', 'Day'} , ...
            {'IPI', 'nanmedian(x)', 'name','IPI1', 'subset',IPI.prsnumb==1}, ...
            {'IPI', 'nanmedian(x)', 'name','IPI2', 'subset',IPI.prsnumb==2}, ...
            {'IPI', 'nanmedian(x)', 'name','IPI3', 'subset',IPI.prsnumb==3}, ...
            {'IPI', 'nanmedian(x)', 'name','IPI4', 'subset',IPI.prsnumb==4});
        
        % plot data
        figure('Name', 'RT IPI'); set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        subplot(2,2,1); plt.scatter(T.RT, IPI.IPI1, 'subset',ismember(T.Horizon,2:4) & T.Day==1); xlabel('RT (ms)'); ylabel('IPI 1 (ms)'); axis image; title('Day 1');
        subplot(2,2,2); plt.scatter(T.RT, IPI.IPI1, 'subset',ismember(T.Horizon,2:4) & T.Day==2); xlabel('RT (ms)'); ylabel('IPI 1 (ms)'); axis image; title('Day 5');
        
        subplot(2,2,3); plt.scatter(T.RT, IPI.IPI1+IPI.IPI2, 'subset',ismember(T.Horizon,2:4) & T.Day==1); xlabel('RT (ms)'); ylabel('IPI 1-2 (ms)'); axis image; title('Day 1');
        subplot(2,2,4); plt.scatter(T.RT, IPI.IPI1+IPI.IPI2, 'subset',ismember(T.Horizon,2:4) & T.Day==2); xlabel('RT (ms)'); ylabel('IPI 1-2 (ms)'); axis image; title('Day 5');
        
        plt.match('both');
        
        % out
        varargout = {T}; %return main structure
        
    otherwise
        error('no such case!')
end