function se2_compareExp(Dall1 , Dall2 , what)
% compare experiments
%% adding eperiments 1 and 2 together
% load('se1_all.mat') --> Dall1
% load('se2_alldata.mat') --> Dall2

switch what
    case 'MT'
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANDOM
        D1 = getrow(Dall1 , Dall1.Horizon == 13 & ~Dall1 .isError & Dall1.seqNumb == 0);
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D1 = normData(D1 , {'MT'});
        D1.exp = ones(size(D1.BN ));
        
        Dall2.exp = 2*ones(size(Dall2.TN));
        D = getrow(Dall2 , Dall2.Horizon == 13 & ~Dall2 .isError & Dall2.seqNumb == 0);
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D = normData(D , {'MT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        figure('color' , 'white')
        subplot(221)
        title('Random Sequences')
        colorz = {[255, 112, 77]/255' , [102, 153, 153]/255'};
        lineplot([T.BN] , T.normMT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [3000 6000],'YTick' , [3000:500:5500] ,...
            'YTickLabels' , [3:.5:5.5] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Execution time [s]')
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STRUCTURED
        D1 = getrow(Dall1 , Dall1.Horizon == 13 & ~Dall1 .isError & ismember(Dall1.seqNumb  , [1:6]));
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D1 = normData(D1 , {'MT'});
        D1.exp = ones(size(D1.BN ));
        % load('se2_alldata.mat')
        D = getrow(Dall2 , Dall2.Horizon == 13 & ~Dall2 .isError & ismember(Dall2.seqNumb  , [1:6]));
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D = normData(D , {'MT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        subplot(222)
        title('Structured Sequences')
        lineplot([T.BN] , T.normMT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [3000 6000],'YTick' , [3000:500:5500] ,...
            'YTickLabels' , [3:.5:5.5] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Execution time [s]')
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Triplets
        D1 = getrow(Dall1 ,  ~Dall1 .isError & ismember(Dall1.seqNumb  , [103 203 303]));
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D1 = normData(D1 , {'MT'});
        D1.exp = ones(size(D1.BN ));
        % load('se2_alldata.mat')
        D = getrow(Dall2 ,  ~Dall2 .isError & ismember(Dall2.seqNumb  ,  [103 203 303]));
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D = normData(D , {'MT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        subplot(223)
        title('Triplets')
        lineplot([T.BN] , T.normMT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [300 1300],'YTick' , [400:200:1200] ,...
            'YTickLabels' , [0.4:.2:1.2] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Execution time [s]')
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quadruples
        D1 = getrow(Dall1 ,  ~Dall1 .isError & ismember(Dall1.seqNumb  , [104 204 304]));
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D1 = normData(D1 , {'MT'});
        D1.exp = ones(size(D1.BN ));
        % load('se2_alldata.mat')
        D = getrow(Dall2 ,  ~Dall2 .isError & ismember(Dall2.seqNumb  ,  [104 204 304]));
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'MT' , 'nanmean'});
        D = normData(D , {'MT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        subplot(224)
        title('Quadruples')
        lineplot([T.BN] , T.normMT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [300 1300],'YTick' , [400:200:1200] ,...
            'YTickLabels' , [0.4:.2:1.2] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Execution time [s]')  
    case 'RT'
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RANDOM
        Dall1.RT = Dall1.AllPressTimes(:,1);
        Dall2.RT = Dall1.AllPressTimes(:,1);
        
        
        D1 = getrow(Dall1 , Dall1.Horizon == 13 & ~Dall1 .isError & Dall1.seqNumb == 0);
        
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D1 = normData(D1 , {'RT'});
        D1.exp = ones(size(D1.BN ));
        
        Dall2.exp = 2*ones(size(Dall2.TN));
        
        D = getrow(Dall2 , Dall2.Horizon == 13 & ~Dall2 .isError & Dall2.seqNumb == 0);
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D = normData(D , {'RT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        figure('color' , 'white')
        subplot(221)
        title('Random Sequences')
        colorz = {[255, 112, 77]/255' , [102, 153, 153]/255'};
        lineplot([T.BN] , T.normRT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [2000 3200],'YTick' , [2200:200:3000] ,...
            'YTickLabels' , [2.2:.2:3] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Reaction time [s]')
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STRUCTURED
        D1 = getrow(Dall1 , Dall1.Horizon == 13 & ~Dall1 .isError & ismember(Dall1.seqNumb  , [1:6]));
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D1 = normData(D1 , {'RT'});
        D1.exp = ones(size(D1.BN ));
        % load('se2_alldata.mat')
        D = getrow(Dall2 , Dall2.Horizon == 13 & ~Dall2 .isError & ismember(Dall2.seqNumb  , [1:6]));
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D = normData(D , {'RT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        subplot(222)
        title('Structured Sequences')
        lineplot([T.BN] , T.normRT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [2000 3200],'YTick' , [2200:200:3000] ,...
            'YTickLabels' , [2.2:.2:3] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Reaction time [s]')
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Triplets
        D1 = getrow(Dall1 ,  ~Dall1 .isError & ismember(Dall1.seqNumb  , [103 203 303]));
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D1 = normData(D1 , {'RT'});
        D1.exp = ones(size(D1.BN ));
        % load('se2_alldata.mat')
        D = getrow(Dall2 ,  ~Dall2 .isError & ismember(Dall2.seqNumb  ,  [103 203 303]));
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D = normData(D , {'RT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        subplot(223)
        title('Triplets')
        lineplot([T.BN] , T.normRT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [2000 2500],'YTick' , [2000:100:2500] ,...
            'YTickLabels' , [2:.1:2.5] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Reaction time [s]')
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Quadruples
        D1 = getrow(Dall1 ,  ~Dall1 .isError & ismember(Dall1.seqNumb  , [104 204 304]));
        D1 = tapply(D1 , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D1 = normData(D1 , {'RT'});
        D1.exp = ones(size(D1.BN ));
        % load('se2_alldata.mat')
        D = getrow(Dall2 ,  ~Dall2 .isError & ismember(Dall2.seqNumb  ,  [104 204 304]));
        D = tapply(D , {'SN' , 'BN' , 'Day'} , {'RT' , 'nanmean'});
        D = normData(D , {'RT'});
        D.exp = 2*ones(size(D.BN ));
        
        T = addstruct(D1 , D);
        
        subplot(224)
        title('Quadruples')
        lineplot([T.BN] , T.normRT , 'plotfcn' , 'nanmean',...
            'split', T.exp , 'linecolor' , colorz,...
            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'exp1'  , 'exp2'} );
        set(gca,'FontSize' , 18 , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [2000 2500],'YTick' , [2000:100:2500] ,...
            'YTickLabels' , [2:.1:2.5] , 'YGrid' , 'on');
        xlabel('Training Block')
        ylabel('Reaction time [s]')
end

