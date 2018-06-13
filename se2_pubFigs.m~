
function out  = se2_pubFigs(Dall , what, nowWhat , varargin)

subj_name = {'AT1' , 'CG1' , 'HB1' , 'JT1' , 'CB1' , 'YM1' , 'NL1' , 'SR1' , 'IB1' , 'MZ1' , 'DW1','RA1' ,'CC1' 'DK1' , 'JM1' , 'All'};
%% Define defaults
subjnum = length(subj_name); % all subjects
Repetition = [1 2];
poolDays = 0;
MaxIter = 100;
%% Deal with inputs
c = 1;
while(c<=length(varargin))
    switch(varargin{c})
        case {'subjnum'}
            % define the subject
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Repetition'}
            % Repetitions to include
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'poolDays'}
            % pool together days 2,3 and days 4 5
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'MaxIter'}
            % maximum number of iterations for exponential fitting
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'isSymmetric'}
            % gaze filed around a digit symmetric or not(0 no/1 yes);
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'dayz'}
            % days to consider -> cell;
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

%%
% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye2/analyze';          %iMac



if subjnum == length(subj_name)
    subjnum = 1:length(subj_name)-1;
else
    Dall  = getrow(Dall , ismember(Dall.SN , subjnum));
end


days  = {1 ,2 ,3 ,4 ,5,[1:5] ,[2:5] [2:3] [4:5],[3:5]};
if ~exist('dayz')
    if poolDays
        dayz = {1 [3] [5]};
        daylab = {'Day 1' , 'Days 2,3' , 'Days 4,5'}
    else
        dayz = {[1] [2] [3] [4] [5]};
        daylab = {'Day 1' , 'Day 2' , 'Day 3' , 'Day 4' , 'Day 5'};
    end
else
    daylab = {};
    for dd = 1:length(dayz)
       daylab = [daylab ['Day(s) ' , num2str(dayz{dd})]]; 
    end
end

clear tempcol
c1 = [255, 153, 179]/255; % Random Red Tones
ce = [153, 0, 51]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 7);
end
for i = 1:length(tempcol)
    colz{i,1} = tempcol(: , i)';
end

clear tempcol
c1 = [153, 194, 255]/255; % structured blue tones
ce = [0, 0, 102]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 7);
end
for i = 1:length(tempcol)
    colz{i,2} = tempcol(: , i)';
    avgCol{i} = mean([colz{i,2} ; colz{i,1}],1);
end


colIPI(:,1) = colz(:,1);    % Random IPIs, same as Random sequences
clear tempcol
c1 = [214, 153, 255]/255;   % Between chunk Fuscia tones
ce = [61, 0, 102]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 7);
end
for i = 1:length(tempcol)
    colIPI{i,2} = tempcol(: , i)';
end
clear tempcol
c1 = [179, 255, 153]/255; % middle chunk green tones
ce = [19, 77, 0]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 7);
end
for i = 1:length(tempcol)
    colIPI{i,3} = tempcol(: , i)';
end

clear tempcol
c1 = [255, 179, 255]/255; % last chunk fuscia tones
ce = [102, 0, 102]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 7);
end
for i = 1:length(tempcol)
    colIPI{i,4} = tempcol(: , i)';
end


scol = [200 200 200]/255;
ecol = [0, 0 0]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(scol(rgb),ecol(rgb) , 7);
end

for i = 1:length(tempcol)
    horzcolor{i,1} = tempcol(: , i)';
end

switch what
    
    case 'MT'
        ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1:6]) & ~Dall.isError);
        ANA.seqNumb(ANA.seqNumb >=1) = 1;
        
        ANA = getrow(ANA , ANA.MT <= 9000 );
        
        for d = 1:length(dayz)
            ANA.Day(ismember(ANA.Day , dayz{d})) = d;
        end
        MT = tapply(ANA , {'Horizon' , 'BN' , 'seqNumb' , 'SN' , 'Day'} , {'MT'});
        MT = normData(MT , {'MT'});
        % segments
        ANA = getrow(Dall , Dall.isgood & ~ismember(Dall.seqNumb , [0 1:6]) & ~Dall.isError);
        ANA.seqNumb(ismember(ANA.seqNumb ,[103 203 303])) = 3;
        ANA.seqNumb(ismember(ANA.seqNumb ,[104 204 304])) = 4;
        
        
        ANA = getrow(ANA , ANA.MT <= 9000 );
        
        for d = 1:length(dayz)
            ANA.Day(ismember(ANA.Day , dayz{d})) = d;
        end
        MTseg = tapply(ANA , {'BN' , 'seqNumb' , 'SN'} , {'MT'});
        MTseg = normData(MTseg , {'MT'});
        
       
        switch nowWhat
            case 'RandvsStructCommpare'
                figure('color' , 'white');
                H = unique(MT.Horizon);
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d);
                    colorz = colz(d,1:2);
                    lineplot([MT.Horizon] , MT.normMT , 'plotfcn' , 'nanmean',...
                        'split', MT.seqNumb , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 2 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                        'markersize' , 6, 'markercolor' , colorz , 'leg' , {'Random'  , 'Structured'}  , 'subset' , ismember(MT.Day , dayz{d}) &  ismember(MT.seqNumb , 0));
                    set(gca,'FontSize' , 18 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [3000 7500],'YTick' , [4000 5000 6000 7000] ,...
                        'YTickLabels' , [4 5 6 7] , 'YGrid' , 'on','XLim' , [1 13]);
                    xlabel('Viewing window size' )
                    ylabel('Execution time [s]')
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end
            case 'RandStructAcrossDays'
                figure('color' , 'white');
                H = unique(MT.Horizon);
                for sn = 0:1
                    subplot(2,1 , sn+1);
                    colorz = colz(:,sn+1);
                    lineplot([MT.Horizon] , MT.normMT , 'plotfcn' , 'nanmean',...
                        'split', MT.Day , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 1 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                        'markersize' , 5, 'markercolor' , colorz , 'leg' , daylab  , 'subset' , ismember(MT.seqNumb , sn));
                    set(gca,'FontSize' , 7 , ...
                        'GridAlpha' , .1 , 'Box' , 'off' , 'YLim' , [3500 7500],'YTick' , [4 :1000:7000 ] ,...
                        'YTickLabels' , [4 :7]);% , 'YGrid' , 'on');
                    xlabel('Viewing window size','FontSize' , 9 )
                    ylabel('Execution time [s]','FontSize' , 9)
                end
            case 'BoxFirstLastDays'
                %% THE BOX PLOTS
                dayz = {[1] [5]};
                horz = {[1] [2] [3] [4] [5] [6;13]};
                MT.Horizon(MT.Horizon>6) = 6;
                
                figure('color' , 'white');
                M = getrow(MT , ismember(MT.Day , [1 5]));
                M.Day(ismember(M.Day , [5])) = 2;
                for d = 1:2
                    subplot(2,1,d)
                    myboxplot(M.Horizon ,M.MT, 'notch',1 ,'plotall',0, 'fillcolor',{colz{dayz{d},1} colz{dayz{d},2}},...
                        'linecolor',{colz{dayz{d},1} colz{dayz{d},2}},...
                        'whiskerwidth',3,'split' , [M.seqNumb] , 'subset' , M.Day == d);
                    
                    set(gca,'FontSize' , 18 ,  ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [1500 9000],...
                        'YTick' , [3000 4000 5000 6000 7000 8000] , 'YTickLabels' , [3 4 5 6 7 8]);%,'XTick' , xs , 'XTickLabel' , xs_labs ,);
                    ylabel('Sec','FontSize' , 21)
                end

            case 'LearningEffectShade'
                H = unique(MT.Horizon);
                for sn = 0:1
                    figure('color' , 'white');
                    colorz = colz(4,sn+1);
                    hcolorz = horzcolor(3,1);
                    lineplot([MT.Day] , MT.normMT , 'plotfcn' , 'nanmean',...
                        'split', MT.Horizon , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,hcolorz,...
                        'linewidth' , 3 , 'markertype' , {'o' , 's' , '<' , '*' , 'x' , '+' , 'd' , '>' , 'p'} , 'markerfill' , colorz,...
                        'markersize' , 15, 'markercolor' , colorz , 'leg' , 'auto'  , 'subset' , ismember(MT.seqNumb , sn));
                    set(gca,'FontSize' , 18 , 'XTick' , [1:length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YLim' , [3500 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                    xlabel('Viewing window size' )
                    ylabel('Execution time [s]')
                end
            case 'compareLearning'
                figure('color' , 'white');
                H = unique(MT.Horizon);
%                 MT.Horizon(MT.Horizon>6) = 6;
                colorz = colz(3,:);
                lineplot([ MT.Horizon MT.Day] , MT.normMT , 'plotfcn' , 'nanmean',...
                    'linecolor' , colorz,...
                    'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                    'linewidth' , 1 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                    'markersize' , 3, 'markercolor' , colorz , 'leg' , 'auto' ,'subset', MT.seqNumb==0 );% ,'split', MT.seqNumb, );
                set(gca,'FontSize' , 7  , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [3500 7500],'YTick' ,...
                    [4000 5000 6000 7000] , 'YTickLabels' , [4 5 6 7]);% , 'YGrid' , 'on');
                xlabel('Viewing window size' ,'FontSize' , 9)
                ylabel('Execution time [s]','FontSize' , 9)
            case 'subjEffectiveHorizon'
                seqN = {[0] , [1 2]};
                allcount = 1;
                dayz = {[1] [2 3] [4 5]};
                daylab = {'Early training (day 1)' 'Mid-training (day 3)' 'Trained (days 4, 5)'};
                subjs  = unique(Dall.SN);
                poolpress = {[1:14]};
                for pp = 1:length(poolpress)
                    Dall.MT = Dall.AllPressTimes(:,poolpress{pp}(end)) -Dall.AllPressTimes(:,poolpress{pp}(1));
                    for sn = 1:length(subjs)
                        dcount = 1;
                        for d  = 1:length(dayz)
                            for sq = 1:length(seqN)
                                EH.Day(allcount,1) = d;
                                EH.SN(allcount,1) = sn;
                                EH.sq(allcount,1) = sq;
                                EH.poolpress(allcount,1) = pp;
                                for h = 1:6
                                    stats = se2_SigTest(Dall , 'MT' , 'seqNumb' , seqN{sq} , 'Day' , dayz{d} , 'Horizon' , [h:13],...
                                        'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
                                        'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , subjs(sn));
                                    pval{sq}(sn,h,dcount) = stats(1);
                                end
                                temp = squeeze(pval{sq}(sn,:,dcount));
                                if ~isempty(find(temp>0.05 ,1 , 'first'))
                                    EH.effH(allcount,1) = find(temp>0.05 ,1 , 'first');
                                else
                                    EH.effH(allcount,1) = 7;
                                end
                                allcount = allcount+1;
                            end
                            dcount = dcount+1;
                        end
                    end
                end
                
                for sn = 1:2
                    figure('color' , 'white')
                    colorz = colz(1:length(dayz) , sn);
                    barplot([EH.sq] ,EH.effH , 'split' ,  EH.Day , 'plotfcn' , 'mean',...
                        'facecolor' , colorz,'edgecolor' , 'none',...
                        'errorwidth' , 1 ,'leg' , daylab , 'subset' ,EH.sq == sn);
                    ylabel('Plannig horizon')
                    set(gca , 'FontSize' , 7 , 'YLim' , [1 4] , 'YTick' , [1 2 3 ] , 'XtickLabels' , {})
                end
                
                for sn = 1:2
                    figure('color' , 'white')
                    colorz = colz(1:length(poolpress) , sn);
                    barplot([EH.Day] ,EH.effH , 'split' ,  EH.poolpress , 'plotfcn' , 'mean',...
                        'facecolor' , colorz,'edgecolor' , 'none',...
                        'errorwidth' , 1 ,'leg' , daylab , 'subset' ,EH.sq == sn);
                    hold on
                    ylabel('Plannig horizon')
                    set(gca , 'FontSize' , 18 , 'YLim' , [1 4] , 'YTick' , [1 2 3 ] , 'XtickLabels' , {})
                end
            case 'subjEffectiveHorizonThresh'
                ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1:6]) & ~Dall.isError);
                ANA.seqNumb(ANA.seqNumb >=1) = 1;
                seqN = {[0] , [1 2]};
                ANA = getrow(ANA , ANA.MT <= 9000 );
                for d = 1:length(dayz)
                    ANA.Day(ismember(ANA.Day , dayz{d})) = d;
                end
%                 dayz = {[1] [3] [5]};
                daylab = {'Early training (day 1)' 'Mid-training (day 3)' 'Trained (day 5)'};
                subjs  = unique(Dall.SN);
                
                Temp = ANA;
                Temp.Horizon(Temp.Horizon>7) = 7;
                Temp = tapply(Temp , {'Horizon' , 'SN' , 'seqNumb' , 'Day'} , {'MT' , 'nanmedian(x)'});
                thresh = [.01 .03 .05 .07  .09 .11];
                figure('color' , 'white')
                fcount = 1;
                for t = 1:length(thresh)
                    EH = [];
                    allcount = 1;
                    for sn = 1:length(subjs)
                        for d  = 1:length(dayz)
                            for sq = 1:length(seqN)
                                EH.Day(allcount,1) = d;
                                EH.SN(allcount,1) = sn;
                                EH.sq(allcount,1) = sq;
                                Th = getrow(Temp , Temp.Horizon~=7 & Temp.SN == sn & ismember(Temp.seqNumb , seqN{sq}) & ismember(Temp.Day ,d));
                                Tfull = getrow(Temp , Temp.Horizon==7 & Temp.SN == sn & ismember(Temp.seqNumb , seqN{sq}) & ismember(Temp.Day ,d));
                                h = find(Th.MT<(1+thresh(t))*Tfull.MT , 1 , 'first');
                                if isempty(h)
                                    EH.effH(allcount,1) = 6;
                                else
                                    EH.effH(allcount,1) = h-1;
                                end
                                allcount = allcount+1;
                            end
                        end
                    end
                    for sn = 1:2
                        subplot(length(thresh),2,fcount)
                        colorz = colz(1:length(dayz) , sn);
                        barplot([EH.sq] ,EH.effH , 'split' ,  EH.Day , 'plotfcn' , 'mean',...
                            'facecolor' , colorz,'edgecolor' , 'none',...
                            'errorwidth' , 1  , 'subset' ,EH.sq == sn);%,'leg' , daylab);% & ~ismember(EH.SN , [3,9]));
                        hold on
                        ylabel('Plannig horizon')
                        set(gca , 'FontSize' , 18 , 'YLim' , [1 4.5] , 'YTick' , [1 2 3 ] , 'XtickLabels' , {})
                        fcount = fcount+1;
                    end
                end
        end
        out = [];
    case 'IPI'
        structNumb = [1 2];
        out = [];
        %         plotfcn = input('nanmean or nanmean?' , 's');
        %% IPIs vs horizon
        % this is the output of the case: 'transitions_All' that is saved to disc
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 , structNumb])  & ~Dall.isError);
        ANA.IPI = [ANA.AllPressTimes(:,1) - 1500 , ANA.IPI];
        ANA.seqNumb(ANA.seqNumb>=2) = 1;
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
            if ismember(ANA.seqNumb(tn) , [1:6])
                ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                a = find(ANA.ChunkBndry(tn , :));
                ANA.ChunkBndry(tn , a(2:end)-1) = 3;
                ANA.ChunkBndry(tn ,end) = 3;
                ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
            else
                ANA.ChunkBndry(tn , :) = zeros(size(ANA.ChnkArrang(tn , :)));
            end
            ANA.ChunkBndry(tn , 1) = -1; % -1 for IPI1 which is RT
        end
        for tn  = 1:length(ANA.TN)
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,14);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,14);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,14);
            ANA.IPI_prsnumb(tn , :) = [0 :13];
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,14);
            ANA.IPI_BN(tn , :) = ANA.BN(tn)*ones(1,14);
        end
        IPItable.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPItable.ChunkBndry = reshape(ANA.ChunkBndry , numel(ANA.IPI) , 1);
        IPItable.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        IPItable.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        IPItable.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        IPItable.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        IPItable.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        IPItable.BN = reshape(ANA.IPI_BN , numel(ANA.IPI) , 1);
        
        A = [];
        for d = 1:length(dayz)
            T = getrow(IPItable , ismember(IPItable.Day , dayz{d}));
            T.Day = ones(size(T.Day))*d;
            A = addstruct(A , T);
        end
        IPItable = A;
        
        switch nowWhat
            case 'sigIPIvssteadystate'
                H = {[1] [2] [3] [4] [5] [6] [7] [8] [13]};
                for d = 1:length(dayz)
                    pval{d} = nan(length(H),13);
                    for h = 1:length(H)
                        for pn = [1:4 , 10:13]
                            stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , [0] , 'Day' , dayz{d} , 'Horizon' , [H{h}],...
                                'PoolDays' , 1,'whatIPI','ipistoEachother','PoolSequences' , 0 ,...
                                'PoolHorizons' , [],'ipiOfInterest' , {[pn] [5] [6] [7] [8] [9]} , 'poolIPIs' , 0 , 'subjnum' , [1:13]);
                            if stats.eff(2).p<=0.05
                                pval{d}(h , pn) = stats.eff(2).p;
                            end
                            
                        end
                    end
                end
                figure('color' , 'white')
                colormap cool
                for d = 1:length(dayz)
                    temp = 100 * ~isnan(pval{d});
                    temp(:,5:9) = 200;
                    subplot(1,3,d)
                    imagesc(temp);
                    y = [1.5:1:8.5];
                    x = [1.5:1:12.5];
                    hold on
                    set(gca,'xtick', x, 'XTickLabel' , {} , 'ytick', y , 'YTickLabel', {});
                    set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'k', 'ycolor', 'k');
                    grid on
                    set(gca , 'GridAlpha' , 1)
                    xlabel('IPI number')
                    ylabel('Window size')
                end
            case 'IPIFullDispsplitDay'
                
                horz = {[1] [2] [3] [4] [5] [6:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5' , '6 - 13'} , 1 , length(dayz));
                
                IPIs  = IPItable;
                % pool last and within
                IPIs.ChunkBndry(IPIs.ChunkBndry == 3) = 2;
                for d = 1:length(dayz)
                    IPIs.Day(ismember(IPIs.Day , dayz{d})) = d;
                end
%                 IPIs.prsnumb(ismember(IPIs.prsnumb , [4:10])) = 4;
                IPIs  = tapply(IPIs , {'Horizon' , 'Day' ,'SN' , 'prsnumb' , 'seqNumb'} , {'IPI' , 'nanmean(x)'});
                IPIs = normData(IPIs , {'IPI'});
                IPIs = getrow(IPIs , ismember(IPIs.Day , [1,length(dayz)]));
                for sqn = 0:1
                    colorz = colz([1,end] , sqn+1);
                    figure('color' , 'white');
                    lineplot([ IPIs.Horizon IPIs.prsnumb] , IPIs.normIPI , 'plotfcn' , 'nanmean',...
                        'split', [ IPIs.Day ] , 'subset', ismember(IPIs.seqNumb  , sqn) & ismember(IPIs.prsnumb ,  [1:13]), 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
                        'linewidth' , .5 , 'markertype' , {'.' , '.'}  , 'markerfill' , colorz,...
                        'markersize' , 10, 'markercolor' , colorz , 'leg' , daylab([1,end]));
                    set(gca,'FontSize' , 7,'GridAlpha' , .2 , 'Box' , 'off',...
                        'YLim' , [150 650] , 'YTick' , [200 :100:600],...
                        'YTickLabel' , [0.2 :0.1: 0.6] ,'XTickLabel' , []);
                    ylabel('Inter-press interval time [s]','FontSize' , 7)
                end
            case 'IPILearningPlacement'
                horz = {[1] [2] [3] [4] [5] [6:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5' , '6 - 13'} , 1 , length(dayz));
                
                IPIs  = IPItable;
                IPIs.Horizon(IPIs.Horizon==13) = 9;
                % pool last and within
                IPIs.ChunkBndry(IPIs.ChunkBndry == 3) = 2;
                for d = 1:length(dayz)
                    IPIs.Day(ismember(IPIs.Day , dayz{d})) = d;
                end
                ipiOfInterest = {[1:2] , [5:9] [12:13]};
                A = [];
                for n = 1:length(ipiOfInterest)
                    temp = getrow(IPIs , ismember(IPIs.prsnumb , ipiOfInterest{n}));
                    temp.IPIPlace = n*ones(size(temp.prsnumb));
                    A = addstruct(A , temp);
                end
                IPIs  = tapply(A , {'Horizon' , 'Day' ,'SN' , 'IPIPlace' , 'seqNumb'} , {'IPI' , 'nanmean(x)'});
                IPIs = normData(IPIs , {'IPI'});
%                 IPIs = getrow(IPIs , ismember(IPIs.Day , [1,length(dayz)]));
                for sqn = 0
                    colo = colIPI([1 5] , 3:-1:1);
                    figure('color' , 'white');
                    hold on
                    for d = 1:length(dayz)
                        if d >1
                            IPIs.Horizon = IPIs.Horizon + 9.5;
                        end
                        colorz = colo(d,:);
                        lineplot([ IPIs.Horizon ] , IPIs.normIPI , 'plotfcn' , 'nanmean',...
                            'split', [  IPIs.IPIPlace] , 'subset', ismember(IPIs.seqNumb  , sqn) & ismember(IPIs.Day  , d), 'linecolor' , colorz,...
                            'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
                            'linewidth' , 1 , 'markertype' , {'o' , '>' , 'd'}  , 'markerfill' , colorz,...
                            'markersize' , 4, 'markercolor' , colorz , 'leg' , 'auto');
                    end
                    set(gca,'FontSize' , 7,'GridAlpha' , .2 , 'Box' , 'off',...
                        'YLim' , [200 600] , 'YTick' , [300 :100:500],...
                        'YTickLabel' , [0.3 :0.1: 0.5] ,'XTick' , [1:18],...
                        'XTickLabel' , repmat({'W=1','W=2','W=3','W=4','W=5','W=6','W=7','W=8','W=13'} ,1 , length(dayz)),...
                        'XTickLabelRotation' , 30);
                    ylabel('Inter-press interval time [s]','FontSize' , 9)
                end
            case 'IPIFullDispsplitseqNumb'
                
                horz = {[1] [2] [3] [4] [5] [6:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5' , '6 - 13'} , 1 , length(dayz));
                IPIs  = IPItable;
                % pool last and within
                IPIs.ChunkBndry(IPIs.ChunkBndry == 3) = 2;
                for d = 1:length(dayz)
                    IPIs.Day(ismember(IPIs.Day , dayz{d})) = d;
                end
                IPIs  = tapply(IPIs , {'Horizon' , 'Day' ,'SN' , 'prsnumb' , 'seqNumb'} , {'IPI' , 'nanmean(x)'});
                IPIs = normData(IPIs , {'IPI'});
                figure('color' , 'white');
                for  d = 1:length(dayz)
                    colorz = colz(d , 1:2);
                    for h = 1:length(horz)
                        subplot(length(dayz),length(horz),(d-1)*length(horz) + h)
                        lineplot([IPIs.prsnumb] , IPIs.normIPI , 'plotfcn' , 'nanmean',...
                            'split', [ IPIs.seqNumb ] , 'subset',ismember(IPIs.Horizon  , horz{h}) & ismember(IPIs.Day  , d), 'linecolor' , colorz,...
                            'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
                            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                            'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Random' , 'Strutured'});
                        set(gca,'FontSize' , 18,'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTick'  , [1:13],...
                            'YLim' , [150 590] , 'YTick' , [200 :100:500],...
                            'YTickLabel' , [0.2 :0.1: 0.5]);
                        if h >1
                            set(gca,'YColor' , 'none');
                        end
                        ylabel('Inter-press interval time [s]','FontSize' , 20)
                        title(['Window size = ' , hlab{h}] ,'FontSize' , 24)
                    end
                end
            case 'IPIFullDispsplitHorizon'
                
                horz = {[1] [2] [3] [4] [5] [6] [7:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5' , '6 - 13'} , 1 , length(dayz));
                IPIs  = IPItable;
                % pool last and within
                IPIs.ChunkBndry(IPIs.ChunkBndry == 3) = 2;
                for d = 1:length(dayz)
                    IPIs.Day(ismember(IPIs.Day , dayz{d})) = d;
                end
                IPIs  = tapply(IPIs , {'Horizon' , 'Day' ,'SN' , 'prsnumb' , 'seqNumb'} , {'IPI' , 'nanmean(x)'});
                IPIs = normData(IPIs , {'IPI'});
                figure('color' , 'white');
                for  sqn = 0:1
                    colorz = colz(: , sqn+1);
                    for d = 1:length(dayz)
                        subplot(2,length(dayz),sqn*length(dayz) + d)
                        lineplot([IPIs.prsnumb] , IPIs.normIPI , 'plotfcn' , 'nanmean',...
                            'split', [ IPIs.Horizon ] , 'subset',ismember(IPIs.seqNumb  , sqn) & ismember(IPIs.Day  , d), 'linecolor' , colorz,...
                            'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
                            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                            'markersize' , 10, 'markercolor' , colorz);% , 'leg' , 'auto');
                        set(gca,'FontSize' , 18,'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTick'  , [1:13],...
                            'YLim' , [150 590] , 'YTick' , [200 :100:500],...
                            'YTickLabel' , [0.2 :0.1: 0.5]);
                        if d >1
                            set(gca,'YColor' , 'none');
                        end
                        ylabel('Inter-press interval time [s]','FontSize' , 20)
                        title(daylab{d} ,'FontSize' , 24)
                    end
                end
            case 'compareLearning'
                horz = {[1] [2] [3] [4] [5] [7:13]};
                
                IPIs=  getrow(IPItable,ismember(IPItable.prsnumb , [4:10]));
                % pool last and within
                IPIs.ChunkBndry(IPIs.ChunkBndry == 3) = 2;
                IPIs  = tapply(IPIs , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI' , 'nanmean(x)'});
                IPIs = normData(IPIs , {'IPI'});
                IPIs.Horizon(IPIs.Horizon>7) = 7;
                colorz = colIPI(3,1:3);
                figure('color' , 'white');
                hlab = {'1' , '2' , '3' , '4'  , '5' , '6' , '7-13'};
                lineplot([IPIs.Day IPIs.Horizon] , IPIs.normIPI , 'plotfcn' , 'nanmean',...
                    'split', [ IPIs.ChunkBndry ] , 'linecolor' , colorz,...
                    'errorcolor' , colorz , 'errorbars' , {'shade'}  , 'shadecolor' ,colorz,...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                    'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Random'  , 'First' , 'Middle'});
                
                set(gca,'FontSize' , 18,'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTickLabel' , ...
                    hlab , 'YLim' , [250 600] , 'YTick' , [300 400 500 600],...
                    'YTickLabel' , [0.3 0.4 0.5 0.6]);
                ylabel('Inter-press interval time [s]','FontSize' , 20)
                xlabel('Window size' ,'FontSize' , 20)
                title('Days  1      2       3      4       5')
            case 'compareLearning_histogram'
                horz = {[1] [2] [3] [4] [5] [6:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5' , '6 - 13'} , 1 , length(dayz));
                IPIs = IPItable;
                IPIs.Horizon(IPIs.Horizon>6) = 6;
                
                figure('color' , 'white')
                fcount = 1;
                for d = 1:length(dayz)
                    for h = 1:length(horz)
                        subplot(length(dayz) , length(horz) , fcount)
                        hold on
                        for ch = 0:2
                            temp = getrow(IPIs , ismember(IPIs.ChunkBndry ,ch) & ismember(IPIs.Horizon , horz{h}) & ismember(IPIs.Day , dayz{d}));
                            histogram(temp.IPI);
                            M(ch+1) = nanmean(temp.IPI);
                        end
                        legend({'r' , 'b' , 'w'})
                        fcount = fcount + 1;
                        title(['H = ' , num2str(h) , 'd = ' , num2str(dayz{d}) , ' ' , num2str(M)])
                    end
                end
            case 'percentTotalLearning_IPIplacement'
                Hz = {[1] [2] [3] , [4] [5] [6] [7:9]};
                ANA = IPItable;
%                 ANA.Horizon(ANA.Horizon>Hz{end}(1)) = Hz{end}(1);
                ANA.IPIPlace = ones(size(ANA.prsnumb));
                ipiOfInterest = {[1:3] , [5:9]};
                for n = 1:length(ipiOfInterest)
                    ANA.IPIPlace(ismember(ANA.prsnumb , ipiOfInterest{n})) = n;
                end
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'IPIPlace' , 'seqNumb'} , {'IPI' , 'nanmedian(x)'});
                ANA.percChangeIPI = zeros(size(ANA.IPI));
                Daybenefit = [];
                for d = length(dayz)
                    for sn = 0:1
                        for ip = 1:length(ipiOfInterest)
                            Db1= getrow(ANA , ANA.Day == 1 & ANA.IPIPlace==ip & ANA.seqNumb == sn);
                            Db_d = getrow(ANA , ANA.Day == d & ANA.IPIPlace==ip & ANA.seqNumb == sn);
                            Db_d.percChangeIPI = 100*abs((Db_d.IPI - Db1.IPI)./Db1.IPI);
                            Daybenefit = addstruct(Daybenefit , Db_d);
                        end
                    end
                end
                Daybenefit  = tapply(Daybenefit , {'Horizon' , 'Day' ,'SN' , 'IPIPlace' , 'seqNumb'} , {'percChangeIPI' , 'nanmedian(x)'} );
                Daybenefit = normData(Daybenefit , {'percChangeIPI'});
                for sn = 0:1
                    figure('color' , 'white')
                    colorz = colz(1:2:end,sn+1);
                    lineplot([Daybenefit.Day Daybenefit.Horizon ] , Daybenefit.percChangeIPI , 'plotfcn' , 'nanmean',...
                        'split', Daybenefit.IPIPlace , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 2 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                        'markersize' , 5, 'markercolor' , colorz , 'leg' , {'first 3'  , 'steady state' , 'last 3'} , 'subset' ,  Daybenefit.seqNumb == sn);
                    set(gca,'FontSize' , 18 , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7-13'} ,1, length(dayz)-1),...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [5 30],'YGrid' , 'on','XTickLabelRotation' , 30);
                    xlabel('Veiwing window size (W)' )
                    ylabel('% improvement')
                    title('Reduction in Inter-Press Intervals compared to Day 1' ,'FontSize' , 24)
                end
            case 'subjEffectiveHorizon'
                seqN = {[0] , [1 2]};
                allcount = 1;
                dayz = {[1] [2 3] [4 5]};
                daylab = {'Early training (day 1)' 'Mid-training (day 3)' 'Trained (days 4, 5)'};
                subjs  = unique(Dall.SN);
                poolpress = {[1:2] , [5:9] , [12:13]};
                for pp = 1:length(poolpress)
                    for sn = 1:length(subjs)
                        dcount = 1;
                        for d  = 1:length(dayz)
                            for sq = 1:length(seqN)
                                EH.Day(allcount,1) = d;
                                EH.SN(allcount,1) = sn;
                                EH.sq(allcount,1) = sq;
                                EH.poolpress(allcount,1) = pp;
                                for h = 1:6
                                    stats = se2_SigTest(Dall , 'IPI' , 'seqNumb' , seqN{sq} , 'Day' , dayz{d} , 'Horizon' , [h:13],...
                                        'PoolDays' , 1,'whatIPI','ipistoEachother','PoolSequences' , 0 ,...
                                        'PoolHorizons' , [],'ipiOfInterest' , poolpress(pp) , 'poolIPIs' , 0 , 'subjnum' , subjs(sn));
                                    pval{sq}(sn,h,dcount) = stats(1);
                                end
                                temp = squeeze(pval{sq}(sn,:,dcount));
                                if ~isempty(find(temp>0.05 ,1 , 'first'))
                                    EH.effH(allcount,1) = find(temp>0.05 ,1 , 'first');
                                else
                                    EH.effH(allcount,1) = 7;
                                end
                                allcount = allcount+1;
                            end
                            dcount = dcount+1;
                        end
                    end
                end
                E = getrow(EH , EH.sq == 1 );
                stats = anovaMixed(E.effH  , E.SN ,'within',[E.Day E.poolpress] ,{ 'Day' , 'pp'},'intercept',1) ;
                
                for sn = 1:2
                    figure('color' , 'white')
                    colorz = colz([1,end] , sn);
                    barplot([EH.poolpress] ,EH.effH , 'split' , EH.Day  , 'plotfcn' , 'mean',...
                        'facecolor' , colorz,'edgecolor' , 'none',...
                        'errorwidth' , 1 ,'leg' , daylab , 'subset' ,EH.sq == sn & EH.poolpress~=4 & ~ismember(EH.Day , 2));
                    hold on
                    ylabel('Plannig horizon')
                    set(gca , 'FontSize' , 7 , 'YLim' , [1 3.5] , 'YTick' , [1 2 3 4] , 'XtickLabels' , {})
                end
                lineplot([EH.Day] ,EH.effH , 'split' , EH.poolpress  , 'plotfcn' , 'mean',...
                        'subset' ,EH.sq == sn & EH.poolpress~=4 ,'style_thickline' , 'leg' , 'auto');
                    
        end
    case 'RT'
        Dall.RT = Dall.AllPressTimes(:,1)-1500;
        ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1:6]) & ~Dall.isError);
        ANA.seqNumb(ANA.seqNumb >=1) = 1;
        
        ANA = getrow(ANA , ANA.RT <= 9000 );

        A = [];
        for d = 1:length(dayz)
            T = getrow(ANA , ismember(ANA.Day , dayz{d}));
            T.Day = ones(size(T.Day))*d;
            A = addstruct(A , T);
        end
        ANA = A;
        
        RT = tapply(ANA , {'Horizon' , 'BN' , 'seqNumb' , 'SN' , 'Day'} , {'RT'});
        RT = normData(RT , {'RT'});
        % segments
        ANA = getrow(Dall , Dall.isgood & ~ismember(Dall.seqNumb , [0 1:6]) & ~Dall.isError);
        ANA.seqNumb(ismember(ANA.seqNumb ,[103 203 303])) = 3;
        ANA.seqNumb(ismember(ANA.seqNumb ,[104 204 304])) = 4;
        
        
        ANA = getrow(ANA , ANA.RT <= 9000 );
        
        RTseg = tapply(ANA , {'BN' , 'seqNumb' , 'SN'} , {'RT'});
        RTseg = normData(RTseg , {'RT'});
        
       
        switch nowWhat
            case 'RandvsStructCommpare'
                figure('color' , 'white');
                H = unique(RT.Horizon);
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d);
                    colorz = colz(d,1:2);
                    lineplot([RT.Horizon] , RT.normRT , 'plotfcn' , 'nanmean',...
                        'split', RT.seqNumb , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                        'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Random'  , 'Structured'}  , 'subset' , ismember(RT.Day , dayz{d}));
                    set(gca,'FontSize' , 18 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [600 900],'YTick' , [700 800] ,...
                        'YTickLabels' , [0.7 0.8] , 'YGrid' , 'on','XLim' , [1 13]);
                    xlabel('Viewing window size' )
                    ylabel('Execution time [s]')
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end
            case 'RandStructAcrossDays'
                figure('color' , 'white');
                H = unique(RT.Horizon);
%                 RT.Horizon(RT.Horizon>7) = 7;
                for sn = 0:1
                    figure('color' , 'white');
                    colorz = colz([1, 6],sn+1);
                    RT.DUMMY = RT.Day;
                    lineplot([RT.Day RT.Horizon ] , RT.normRT , 'plotfcn' , 'nanmean',...
                        'split' , RT.DUMMY , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 1 , 'markertype' , repmat({'s'} , 1  , 2) , 'markerfill' , colorz,...
                        'markersize' , 5, 'markercolor' , colorz , 'leg' , daylab  , 'subset' , ismember(RT.seqNumb , sn));
                    set(gca,'FontSize' , 7 , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7' '8' '13'} ,1,length(dayz)), ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [600 850],'YTick' ,  [700 800] ,...
                        'YTickLabels' , [0.7 0.8] , 'YGrid' , 'on');
                    xlabel('Viewing window size' )
                    ylabel('Execution time [s]')
                end
            case 'BoxFirstLastDays'
                %% THE BOX PLOTS
                dayz = {[1] [5]};
                horz = {[1] [2] [3] [4] [5] [6;13]};
                RT.Horizon(RT.Horizon>6) = 6;
                
                figure('color' , 'white');
                M = getrow(RT , ismember(RT.Day , [1 5]));
                M.Day(ismember(M.Day , [5])) = 2;
                for d = 1:2
                    subplot(2,1,d)
                    myboxplot(M.Horizon ,M.RT, 'notch',1 ,'plotall',0, 'fillcolor',{colz{dayz{d},1} colz{dayz{d},2}},...
                        'linecolor',{colz{dayz{d},1} colz{dayz{d},2}},...
                        'whiskerwidth',3,'split' , [M.seqNumb] , 'subset' , M.Day == d);
                    
                    set(gca,'FontSize' , 18 ,  ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [1500 9000],...
                        'YTick' , [3000 4000 5000 6000 7000 8000] , 'YTickLabels' , [3 4 5 6 7 8]);%,'XTick' , xs , 'XTickLabel' , xs_labs ,);
                    ylabel('Sec','FontSize' , 21)
                end

            case 'LearningEffectShade'
                H = unique(RT.Horizon);
                for sn = 0:1
                    figure('color' , 'white');
                    colorz = colz(4,sn+1);
                    hcolorz = horzcolor(3,1);
                    lineplot([RT.Day] , RT.normRT , 'plotfcn' , 'nanmean',...
                        'split', RT.Horizon , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,hcolorz,...
                        'linewidth' , 3 , 'markertype' , {'o' , 's' , '<' , '*' , 'x' , '+' , 'd' , '>' , 'p'} , 'markerfill' , colorz,...
                        'markersize' , 15, 'markercolor' , colorz , 'leg' , 'auto'  , 'subset' , ismember(RT.seqNumb , sn));
                    set(gca,'FontSize' , 18 , 'XTick' , [1:length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YLim' , [3500 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                    xlabel('Viewing window size' )
                    ylabel('Execution time [s]')
                end
            case 'compareLearning'
                figure('color' , 'white');
                H = unique(RT.Horizon);
                RT.Horizon(RT.Horizon>6) = 6;
                colorz = colz(3,:);
                lineplot([ RT.Horizon RT.Day] , RT.normRT , 'plotfcn' , 'nanmean',...
                    'split', RT.seqNumb, 'linecolor' , colorz,...
                    'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                    'markersize' , 10, 'markercolor' , colorz , 'leg' , 'auto' );
                set(gca,'FontSize' , 18  , 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [2500 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                xlabel('Viewing window size' )
                ylabel('Execution time [s]')

                
            case 'subjEffectiveHorizon'
                seqN = {[0] , [1 2]};
                allcount = 1;
                dayz = {[1] [3] [4 5]};
                daylab = {'Early training (day 1)' 'Mid-training (day 3)' 'Trained (day 5)'};
                subjs  = unique(Dall.SN);
%                 Dall.MT = Dall.AllPressTimes(3) -Dall.AllPressTimes(1);
                for sn = 1:length(subjs)
                    dcount = 1;
                    for d  = 1:length(dayz)
                        for sq = 1:length(seqN)
                            EH.Day(allcount,1) = d;
                            EH.SN(allcount,1) = sn;
                            EH.sq(allcount,1) = sq;
                            for h = 1:6
                                stats = se2_SigTest(Dall , 'RT' , 'seqNumb' , seqN{sq} , 'Day' , dayz{d} , 'Horizon' , [h:13],...
                                    'PoolDays' , 1,'whatIPI','WithBetRand','PoolSequences' , 0 ,...
                                    'PoolHorizons' , [],'ipiOfInterest' , [] , 'poolIPIs' , 0 , 'subjnum' , subjs(sn));
                                pval{sq}(sn,h,dcount) = stats(1);
                            end
                            temp = squeeze(pval{sq}(sn,:,dcount));
                            if ~isempty(find(temp>0.05 ,1 , 'first'))
                                EH.effH(allcount,1) = find(temp>0.05 ,1 , 'first');
                            else
                                EH.effH(allcount,1) = 7;
                            end
                            allcount = allcount+1;
                        end
                        dcount = dcount+1;
                    end
                end
                figure('color' , 'white')
                for sn = 1:2
                    subplot(1,2,sn)
                    colorz = colz(1:length(dayz) , sn);
                    barplot([EH.sq] ,EH.effH , 'split' ,  EH.Day , 'plotfcn' , 'mean',...
                        'facecolor' , colorz,'edgecolor' , 'none',...
                        'errorwidth' , 1 ,'leg' , daylab , 'subset' ,EH.sq == sn);
                    hold on
                    ylabel('Plannig horizon')
                    set(gca , 'FontSize' , 18 , 'YLim' , [1 4] , 'YTick' , [1 2 3 ] , 'XtickLabels' , {})
                end
        end
        out = [];
    case 'MT_asymptote'
        Hex = 0;
        out = [];
        plotcoef = 0;
%         Hex = input('What horizons to exclude? (0 = include all)');
        MT = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError );
        MT.seqNumb(MT.seqNumb == 2) = 1;
%         MT  = tapply(MT , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmedian(x)'});
        MT.MT_pred = zeros(size(MT.MT));
        % account for pooling dayz together
        for d = 1:length(dayz)
            MT.Day(ismember(MT.Day , dayz{d})) = d;
        end
        MTs = [];
        coefs = [];
        if poolDays
            sigSeq = [NaN 3 2];
            sigMT  = [3 4 5;4 5 4];
        else
            sigSeq = [NaN 3 3 2 2];
            sigMT  = [3 4 4 6 3;4 3 4 4 4];
        end
        h0 = figure;
        for subjnum = 1:length(subj_name)-1
            for d = 1:length(dayz)
                for sn = 0:1
                    [d subjnum]
                    % Structured
                    id = ismember(MT.SN , subjnum) & ismember(MT.Day , d) & ismember(MT.seqNumb , sn);
                    MTsn = getrow(MT , id);
                    exp_model1 = @(b,x) b(1) + (b(2) - b(1))*exp(-(x-1)/b(3)); % Model Function
                    %                 exp_model1 = @(b,x) b(1)*exp(-(x-1)/b(2)); % Model Function
                    x = MTsn.Horizon';
                    yx = MTsn.MT';                                    % this would be a typical MT vs Horizon vector: [5422 3548 2704 2581 2446 2592 2418 2528 2500]
                    OLS = @(b) sum((exp_model1(b,x) - yx).^2);                % Ordinary Least Squares cost function
                    opts = optimset('MaxIter', MaxIter,'TolFun',1e-5);
                    [B1 Fval] = fminsearch(OLS,[3500 7500  1], opts);        % Use ?fminsearch? to minimise the ?OLS? function
                    %                 [B1 Fval] = fminsearch(OLS,[3500  .1], opts);
                    MTsn.MT_pred = exp_model1(B1,x)';
                    B.b1 = B1(1);
                    B.b2 = B1(2);
                    B.b3 = B1(3);
                    B.seqNumb = sn;
                    B.SN  = subjnum;
                    B.Day = d;
                    B.MT_fitted(1,:) = exp_model1(B1,unique(MT.Horizon))';
                    temp = tapply(MTsn , {'Horizon' } , {'MT' , 'nanmedian(x)'});
                    B.Median_MT(1,:) = temp.MT;
                    B.Horizon(1,:) = unique(MT.Horizon);
                    MTs = addstruct(MTs,MTsn);
                    coefs = addstruct(coefs , B);
                    clear B B1
                end
            end
            if plotcoef
                MTSN = getrow(MTs , MTs.SN == subjnum);
                f1 = figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(2,length(dayz),d)
                    hold on
                    for sn = 0
                        [xx , pp0 ,ee0] = lineplot([MTSN.Horizon] , MTSN.MT , 'subset' , ismember(MTSN.seqNumb , sn)  & ismember(MTSN.Day ,d),...
                            'linecolor' , colIPI{d,sn+1} , 'plotfcn' , 'nanmean');
                    end
                    for sn = 0
                        [xx , pp1 ,ee1]=lineplot([MTSN.Horizon] , MTSN.MT_pred , 'subset' , ismember(MTSN.seqNumb , sn)  & ismember(MTSN.Day ,d),...
                            'linecolor' , colIPI{d,sn+1}, 'plotfcn' , 'nanmean');
                    end
                    xlabel('Viewing Horizon' )
                    set(gca , 'YLim' , [min([pp0 pp1]-[ee0 ee1]) max([pp0 pp1]+[ee0 ee1])])
                end
                hold off
                
                coefSN = getrow(coefs , coefs.SN == subjnum);
                d= d+1;
                subplot(2,length(dayz),d)
                hold on
                for sn = 0:1
                    coefSNch = getrow(coefSN , ismember(coefSN.seqNumb , sn));
                    plot(coefSNch.Day , coefSNch.b1 , 'o-','color' , colz{4,sn+1} )
                end
                title('b1')
                hold off
                
                d= d+1;
                subplot(2,length(dayz),d)
                hold on
                for sn = 0:1
                    coefSNch = getrow(coefSN , ismember(coefSN.seqNumb , sn));
                    plot(coefSNch.Day , coefSNch.b2 , 'o-','color' , colz{4,sn+1} )
                end
                title('b2')
                hold off
                
                d= d+1;
                subplot(2,length(dayz),d)
                hold on
                for sn = 0:1
                    coefSNch = getrow(coefSN , ismember(coefSN.seqNumb , sn));
                    plot(coefSNch.Day , coefSNch.b3 , 'o-','color' , colz{4,sn+1} )
                end
                title('b3')
                hold off
                keyboard
                close(f1)
            end
        end
      
        %%
        
        switch nowWhat
            case 'Actual&fitHorz'
                figure('color' , 'white')
                hold on
                hz = max(unique(MT.Horizon));
                xtick = [];
                for d = [1:length(dayz)]
                    xtick = [xtick ; (d*(hz+2))+coo{sn+1}(:,d)];
                    for sn = 0:1
                        hold on
                        h1 = plotshade((d*(hz+2))+coo{sn+1}(:,d)',Plot{sn+1}(:,d)',err{sn+1}(:,d)','transp' , .5 , 'patchcolor' , colz{d,sn+1} , 'linecolor' , colz{d,sn+1} , 'linewidth' , 3);
                        plot((d*(hz+2))+coo{sn+1}(:,d)',Plot{sn+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,sn+1},'MarkerFaceColor',colz{d,sn+1});
                    end
                end
                set(gca,'FontSize' , 18 , 'XTick' , xtick , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7' '8' '13'} , 1, length(dayz)) , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [min(xtick) max(xtick)], 'YLim' , [3000 7200],'YTick' ,...
                    [ 4000 5000 6000] , 'YTickLabels' , [4 5 6] , 'YGrid' , 'on');
                ylabel('Sec' ,'FontSize' , 21)
                xlabel('Viewing Horizon' , 'FontSize' , 21)
                title(['Actual Chunked MT vs. Random on Day(s) ' , num2str(dayz{d})],'FontSize' , 24)
                
                figure('color' , 'white')
                hold on
                hz = max(unique(MT.Horizon));
                xtick = [];
                for d = 1:length(dayz)
                    xtick = [xtick ; (d*(hz+2))+coo{sn+1}(:,d)];
                    for sn = 0:1
                    hold on
                    h1 = plotshade((d*(hz+2))+coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)',err_pred{sn+1}(:,d)','transp' , .5 , 'patchcolor' , colz{d,sn+1} , 'linecolor' , colz{d,sn+1} , 'linewidth' , 3);
                    plot((d*(hz+2))+coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,sn+1},'MarkerFaceColor',colz{d,sn+1});
                    end
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7' '8' '13'} , 1, length(dayz)) , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [min(xtick) max(xtick)], 'YLim' , [2000 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [4 5 6] , 'YGrid' , 'on');
                ylabel('Sec' )
                xlabel('Viewing Horizon' )
                title(['Fitted Chunked MT vs. Random on Day(s) ' , num2str(dayz{d})])
            case 'Actual&fitDayz'
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                h1 = figure;
                for d = 1:length(dayz)
                    for sn = 0:1
                        ANA = getrow(MTs , ismember(MTs.seqNumb , sn)  & ismember(MTs.Day , d));
                        ANA_pred = getrow(MTs_pred , ismember(MTs_pred.seqNumb , sn)  & ismember(MTs_pred.Day , d));
                        ANA.Horizon(ANA.Horizon>6) = 6;
                        [coo_red{sn+1}(:,d),plot_red{sn+1}(:,d),err_red{sn+1}(:,d)] = lineplot([ANA.Horizon] , ANA.MT );
                        [coo_pred_red{sn+1}(:,d),plot_pred_red{sn+1}(:,d),err_pred_red{sn+1}(:,d)] = lineplot([ANA_pred.Horizon] , ANA_pred.MT_pred );
                    end
                end
                close(h1)
                scol = [200 200 200]/255;
                ecol = [0, 0 0]/255;
                for rgb = 1:3
                    horzcolor(rgb , :) = linspace(scol(rgb),ecol(rgb) , length(unique(ANA.Horizon)));
                end
                %                 horzcolor = transpose([107, 91, 149;254, 178, 54;214, 65, 97;255, 123, 37;128, 206, 214;241, 137, 115]/255);
                
                figure('color' , 'white')
                subplot(211)
                hold on
                xtick = [];
                for i = [1:length(unique(ANA.Horizon))]
                    for sn = 0:1
                        plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_red{sn+1}(i,:),err_red{sn+1}(i,:),'transp' , .6 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , colz{6-i+1,sn+1} , 'linewidth' , 3);
                        plot((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_red{sn+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{6-i+1,sn+1},'MarkerFaceColor',colz{6-i+1,sn+1} , 'LineWidth' , 3);
                    end
                    xtick = [xtick (i-1)*(length(dayz)+1)+[1:length(dayz)]];
                    line([i*(length(dayz)+1) i*(length(dayz)+1)] , [2500 8000] , 'LineWidth' , 3 , 'LineStyle' , ':' , 'color' , [.8 .8 .8])
                end
                set(gca,'FontSize' , 18 , 'XTick' , xtick , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 i*(length(dayz)+1)], 'YLim' , [2500 8000],'YTick' ,...
                    [3000 4000 5000 6000 7000] , 'YTickLabels' , [3 4 5 6 7] , 'YGrid' , 'on',...
                    'XTickLabels' , repmat([1:5] , 1 ,length(unique(ANA.Horizon))));
                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Training Session','FontSize' , 20)
                title('Actual','FontSize' , 24)
                
                subplot(212)
                hold on
                xtick = [];
                for i = [1:length(unique(ANA.Horizon))]
                    for sn = 0:1
                        plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{sn+1}(i,:),err_pred_red{sn+1}(i,:),'transp' , .6 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , colz{6-i+1,sn+1} , 'linewidth' , 3);
                        plot((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{sn+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{6-i+1,sn+1},'MarkerFaceColor',colz{6-i+1,sn+1} , 'LineWidth' , 3);
                    end
                    xtick = [xtick (i-1)*(length(dayz)+1)+[1:length(dayz)]];
                    line([i*(length(dayz)+1) i*(length(dayz)+1)] , [2500 8000] , 'LineWidth' , 3 , 'LineStyle' , ':' , 'color' , [.8 .8 .8])
                end
                set(gca,'FontSize' , 18 , 'XTick' , xtick , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 i*(length(dayz)+1)], 'YLim' , [2500 8000],'YTick' ,...
                    [3000 4000 5000 6000 7000] , 'YTickLabels' , [3 4 5 6 7] , 'YGrid' , 'on',...
                    'XTickLabels' , repmat([1:5] , 1 ,length(unique(ANA.Horizon))));
                ylabel('Sec','FontSize' , 20 )
                xlabel('Training Session','FontSize' , 20)
                title('Fitted','FontSize' , 24)
            case 'Actual&fit%ChangeDayzTotalLearning'
                if poolDays
                    daylab = {'Session 1 to 2,3' , 'Sessions 1 to 4,5'};
                else
                    daylab = {'Session 1 to 2' , 'Session 1 to 3' , 'Session 1 to 4' , 'Session 1 to 5'};
                end
                Hz = {[1] [2] [3] , [4] [5] [6] [7:9]};
                ANA = MTs;
                ANA.Horizon(ANA.Horizon>7) = 7;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'},{'MT_pred' , 'nanmean(x)'});
                ANA.percChangeMT = zeros(size(ANA.MT));
                ANA.percChangeMT_pred = zeros(size(ANA.MT_pred));
                
                Daybenefit = [];
                for d = length(dayz)
                    for sn = 0:1
                        Db1= getrow(ANA , ANA.Day == 1 & ANA.seqNumb==sn);
                        Db_d = getrow(ANA , ANA.Day == d & ANA.seqNumb==sn);
                        Db_d.percChangeMT = 100*abs((Db1.MT - Db_d.MT)./Db1.MT);
                        Db_d.percChangeMT_pred = 100*abs((Db1.MT_pred - Db_d.MT_pred)./Db1.MT_pred);
                        Daybenefit = addstruct(Daybenefit , Db_d);
                    end
                end
                Daybenefit = normData(Daybenefit , {'percChangeMT' , 'percChangeMT_pred'});
                
                figure('color' , 'white')
                for sn = 0:1
                    subplot(1,2,sn+1)
                    colorz = colz([length(dayz)],sn+1);
%                     barplot([Daybenefit.Horizon] , Daybenefit.normpercChangeMT , 'plotfcn' , 'nanmean',...
%                         'split', Daybenefit.Day , 'facecolor' , colorz,...
%                         'edgecolor' , 'none',...
%                         'errorwidth' , 1 ,'leg' , {'Day 1 to 2'  , 'Day 1 to 3' , 'Day 1 to 4' , 'Day 1 to 5'} , 'subset' ,Daybenefit.seqNumb == sn);% & ismember(Daybenefit.Day , [2 5]));
                    
                    lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeMT , 'plotfcn' , 'nanmean',...
                        'split', Daybenefit.Day , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 3 , 'markertype' , {'o' , 's' , '<' , '*'}  , 'markerfill' , colorz,...
                        'markersize' , 15, 'markercolor' , colorz , 'leg' , {'Day 2'  , 'Day 3' , 'Day 4' , 'Day 5'} , 'subset' ,Daybenefit.seqNumb == sn & ismember(Daybenefit.Day , [5]));
                    
                    
                    set(gca,'FontSize' , 18 , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [5 25],'YTick' ,...
                        [10 20 30] , 'YGrid' , 'on');
                    ylabel('% improvement' ,'FontSize' , 20)
                    xlabel('Viewing Window size','FontSize' , 20)
                    title('Reduction in Sequence Execution Time From First to Last Day (Actual)' ,'FontSize' , 24)
                end
                
                figure('color' , 'white')
                colorz = colz(4,1:2);
                    lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeMT_pred , 'plotfcn' , 'nanmean',...
                        'split', Daybenefit.seqNumb , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                        'markersize' , 10, 'markercolor' , colorz , 'leg' , {'Random'  , 'Structured'} );
              
                set(gca,'FontSize' , 18 , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on');
                ylabel('%' ,'FontSize' , 20)
                xlabel('Viewing Window size','FontSize' , 20)
                title('Reduction in Sequence Execution Time From First to Last Day (Fitted)' ,'FontSize' , 24)
            case 'Actual&fit%ChangeDay2Day'
                dayz = {[1] [3]  [5]};
%                 dayz = {[1] [2] [3] [4] [5]};
                Hz = {[1] [2] [3] , [4] [5] [6] [7:13]};
                ANA = MTs;
                ANA.Horizon(ANA.Horizon>7) = 7;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'},{'MT_pred' , 'nanmean(x)'});
                ANA.percChangeMT = zeros(length(ANA.MT) , length(dayz)-1);
                ANA.percChangeMT_pred = zeros(length(ANA.MT) , length(dayz)-1);
                
                Daybenefit = [];
                
                % day to day
                for d = 1:length(dayz)-1
                    Dbd= getrow(ANA , ismember(ANA.Day , dayz{d}) );
                    Db = getrow(ANA , ismember(ANA.Day , dayz{d+1}));
                    Db.percChangeMT = 100*abs((Dbd.MT - Db.MT)./Dbd.MT);
                    Db.percChangeMT_pred = 100*abs((Dbd.MT_pred - Db.MT_pred)./Dbd.MT_pred);
                    Daybenefit = addstruct(Daybenefit , Db);
                end
                % overal 1 -> 5
                for d = length(dayz)
                    Db1= getrow(ANA , ismember(ANA.Day , 1) );
                    Db = getrow(ANA , ismember(ANA.Day , dayz{d}));
                    Db.percChangeMT = 100*abs((Db1.MT - Db.MT)./Db1.MT);
                    Db.percChangeMT_pred = 100*abs((Db1.MT_pred - Db.MT_pred)./Db1.MT_pred);
                    Db.Day = 10+Db.Day;
                    Daybenefit = addstruct(Daybenefit , Db);
                end
                
                Daybenefit = normData(Daybenefit , {'percChangeMT' , 'percChangeMT_pred'});
                if poolDays
                    daylab = {'Session 1 to 3' , 'Sessions 3 to 5' , 'Overal (session 1 to 5)'};
                else
                    daylab = {'Session 1 to 2' , 'Session 2 to 3' , 'Session 3 to 4' , 'Session 4 to 5' , 'Overal ()'};
                end
                xtick = 1:length(unique(ANA.Horizon));
                xtl = repmat({'W=1' , 'W=2' , 'W=3' , 'W=4' , 'W=5' ,'W=6' , 'W=7-13'} , 1, length(dayz));
                xt = xtick;
                for dd = 1:length(dayz)-1
                    xt = [xt , xtick + dd*length(unique(ANA.Horizon))+dd*.5];
                end
                
                for sn = 0%:1
                    colorz = colz(end,sn+1);
                    figure('color' , 'white')
                     barplot([Daybenefit.Day Daybenefit.Horizon] , Daybenefit.normpercChangeMT , 'plotfcn' , 'nanmean',...
                         'subset' ,Daybenefit.seqNumb == sn, 'facecolor' , colorz,'edgecolor' , 'white',...
                        'errorwidth' , 1 ,'leg' , daylab,'XTickLabel' ,{'W=1' , 'W=2' , 'W=3' , 'W=4' , 'W=5' ,'W=6' , 'W=7-13'})

                    
                    
%                     barplot([Daybenefit.Day] , Daybenefit.normpercChangeMT , 'plotfcn' , 'nanmean',...
%                         'split', Daybenefit.Horizon, 'subset' ,Daybenefit.seqNumb == sn, 'facecolor' , colorz,'edgecolor' , 'none',...
%                         'errorwidth' , 1 )
                 
         
%                     hold on
%                     lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeMT , 'plotfcn' , 'nanmean',...
%                         'split', Daybenefit.Day , 'linecolor' , colorz,...
%                         'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
%                         'linewidth' , 1 , 'markertype' , {'o' , 's' , '<' , '*'}  , 'markerfill' , colorz,...
%                         'markersize' , 5, 'markercolor' , colorz , 'subset' ,Daybenefit.seqNumb == sn , 'leg' , 'auto');
                    
                    
                    set(gca,'FontSize' , 7 , 'XTick' , xt , 'XTickLabels' , xtl,...
                        'XTickLabelRotation' , 30, 'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [0 25],'YTick' ,...
                        [10 20 30]);
                    ylabel('%' ,'FontSize' , 9)
                    xlabel('Viewing Window size','FontSize' , 9)
                    title('Day-to-day Improvement in Performance in Different Viewing Window Sizes (W)' ,'FontSize' , 10)
                end
                subplot(212)
                hold on
                Nh = length(unique(Daybenefit.Horizon));
                for sn = 0:1
                    xtick = [];
                    for h = 1:length(unique(Daybenefit.Horizon))
                        Xcoor = (h-1)*(length(dayz))+[1:length(dayz)-1];
                        plotshade(Xcoor,plot_pred_red{sn+1}(:,h)',err_pred_red{sn+1}(:,h)','transp' , .6 , 'patchcolor' , colz{4,sn+1} , 'linecolor' , colz{4,sn+1} , 'linewidth' , 3);
                        plot(Xcoor,plot_pred_red{sn+1}(:,h), '-o' , 'MarkerSize' , 15 , 'color' , colz{4,sn+1},'MarkerFaceColor',colz{4,sn+1} , 'LineWidth' , 4);
                        xtick = [xtick  Xcoor];
                    end
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabels' , repmat(daylab , length(dayz)-1),...
                    'XTickLabelRotation' , 30, 'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [0 max(xtick)], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on');
                ylabel('%' ,'FontSize' , 22)
                xlabel('Viewing Window size','FontSize' , 22)
                title('Day-to-day Improvement in Fitted Performance in Different Viewing Window Sizes (W)' ,'FontSize' , 28)
            case 'Actual&fit%ChangeSeqType'
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                ANA = MTs;
                ANA.Horizon(ANA.Horizon>=6) = 6;
                dayz = {[1] [2 3] [4 5]};
                % pool dayz for better visual
                %                 ANA.Day(ANA.Day==3) = 2;
                %                 ANA.Day(ismember(ANA.Day , [4 5])) = 3;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'},{'MT_pred' , 'nanmean(x)'});
                ANA.percChangeMT = zeros(length(ANA.MT),length(dayz)-1);
                ANA.percChangeMT_pred = zeros(length(ANA.MT_pred),length(dayz)-1);
                Seqbenefit = [];
                for d = 1:length(dayz)
                    Db1= getrow(ANA , ismember(ANA.Day , dayz{d}) & ANA.seqNumb == 1);
                    Db = getrow(ANA , ismember(ANA.Day , dayz{d}) & ANA.seqNumb == 0);
                    Db1.percChangeMT = 100*abs((Db.MT - Db1.MT)./Db.MT);
                    Db1.percChangeMT_pred = 100*abs((Db.MT_pred - Db1.MT_pred)./Db.MT_pred);
                    Seqbenefit = addstruct(Seqbenefit , Db1);
                end
                h1 = figure;
                hold on
                for d = 1:length(dayz)
                    [coo_red{d},plot_red{d},err_red{d}] = lineplot([Seqbenefit.Horizon] , Seqbenefit.percChangeMT, 'subset', ismember(Seqbenefit.Day , dayz{d}));
                    [coo_pred_red{d},plot_pred_red{d},err_pred_red{d}] = lineplot([Seqbenefit.Horizon] , Seqbenefit.percChangeMT_pred, 'subset', ismember(Seqbenefit.Day , dayz{d}));
                end
                close(h1)
                
                figure('color' , 'white')
                hold on
                for d = 1:length(dayz)
                    plotshade([1:length(unique(Seqbenefit.Horizon))],plot_red{d},err_red{d},'transp' , .6 , 'patchcolor' , avgCol{2*d} , 'linecolor' , avgCol{2*d} , 'linewidth' , 3);
                    plot([1:length(unique(Seqbenefit.Horizon))],plot_red{d}, '-o' , 'MarkerSize' , 15 , 'color' , avgCol{2*d},'MarkerFaceColor',avgCol{2*d} , 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Seqbenefit.Horizon))], 'YLim' , [0 15],'YTick' ,...
                    [5 10 15] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Seqbenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4', '5' , '6-13'});
                ylabel('%','FontSize' , 21 )
                xlabel('Viewing Window Size','FontSize' , 21)
                title('Relative Percent Performance Improvement in Structured Sequencs Compared to Random','FontSize' , 24)
                
                figure('color' , 'white')
                hold on
                for d = 1:length(dayz)
                    plotshade([1:length(unique(Seqbenefit.Horizon))],plot_pred_red{d},err_pred_red{d},'transp' , .6 , 'patchcolor' , avgCol{2*d} , 'linecolor' , avgCol{2*d} , 'linewidth' , 3);
                    plot([1:length(unique(Seqbenefit.Horizon))],plot_pred_red{d}, '-o' , 'MarkerSize' , 15 , 'color' , avgCol{2*d},'MarkerFaceColor',avgCol{2*d} , 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Seqbenefit.Horizon))], 'YLim' , [0 15],'YTick' ,...
                    [5 10 15] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Seqbenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5-13'});
                ylabel('%','FontSize' , 22)
                xlabel('Viewing Window Size' ,'FontSize' , 22)
                title('Fitted Percent Relative Performance Improvement in Structured Sequencs Compared to Random','FontSize' , 28)
            case 'Actual&fit%ChangeSeqType'
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                ANA = MT;
                ANA.Horizon(ANA.Horizon>=6) = 6;
                dayz = {[1] [2 3] [4 5]};
                % pool dayz for better visual
                %                 ANA.Day(ANA.Day==3) = 2;
                %                 ANA.Day(ismember(ANA.Day , [4 5])) = 3;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'},{'MT_pred' , 'nanmean(x)'});
                ANA.percChangeMT = zeros(length(ANA.MT),length(dayz)-1);
                ANA.percChangeMT_pred = zeros(length(ANA.MT_pred),length(dayz)-1);
                Seqbenefit = [];
                for d = 1:length(dayz)
                    Db1= getrow(ANA , ismember(ANA.Day , dayz{d}) & ANA.seqNumb == 1);
                    Db = getrow(ANA , ismember(ANA.Day , dayz{d}) & ANA.seqNumb == 0);
                    Db1.percChangeMT = 100*abs((Db.MT - Db1.MT)./Db.MT);
                    Db1.percChangeMT_pred = 100*abs((Db.MT_pred - Db1.MT_pred)./Db.MT_pred);
                    Seqbenefit = addstruct(Seqbenefit , Db1);
                end
                h1 = figure;
                hold on
                for d = 1:length(dayz)
                    [coo_red{d},plot_red{d},err_red{d}] = lineplot([Seqbenefit.Horizon] , Seqbenefit.percChangeMT, 'subset', ismember(Seqbenefit.Day , dayz{d}));
                    [coo_pred_red{d},plot_pred_red{d},err_pred_red{d}] = lineplot([Seqbenefit.Horizon] , Seqbenefit.percChangeMT_pred, 'subset', ismember(Seqbenefit.Day , dayz{d}));
                end
                close(h1)
                
                figure('color' , 'white')
                subplot(211)
                hold on
                for d = 1:length(dayz)
                    plotshade([1:length(unique(Seqbenefit.Horizon))],plot_red{d},err_red{d},'transp' , .6 , 'patchcolor' , avgCol{2*d} , 'linecolor' , avgCol{2*d} , 'linewidth' , 3);
                    plot([1:length(unique(Seqbenefit.Horizon))],plot_red{d}, '-o' , 'MarkerSize' , 15 , 'color' , avgCol{2*d},'MarkerFaceColor',avgCol{2*d} , 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Seqbenefit.Horizon))], 'YLim' , [0 15],'YTick' ,...
                    [5 10 15] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Seqbenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5-13'});
                ylabel('%','FontSize' , 22 )
                xlabel('Viewing Window Size','FontSize' , 22)
                title('Relative Percent Performance Improvement in Structured Sequencs Compared to Random','FontSize' , 28)
                
                subplot(212)
                hold on
                for d = 1:length(dayz)
                    plotshade([1:length(unique(Seqbenefit.Horizon))],plot_pred_red{d},err_pred_red{d},'transp' , .6 , 'patchcolor' , avgCol{2*d} , 'linecolor' , avgCol{2*d} , 'linewidth' , 3);
                    plot([1:length(unique(Seqbenefit.Horizon))],plot_pred_red{d}, '-o' , 'MarkerSize' , 15 , 'color' , avgCol{2*d},'MarkerFaceColor',avgCol{2*d} , 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Seqbenefit.Horizon))], 'YLim' , [0 15],'YTick' ,...
                    [5 10 15] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Seqbenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5-13'});
                ylabel('%','FontSize' , 22)
                xlabel('Viewing Window Size' ,'FontSize' , 22)
                title('Fitted Percent Relative Performance Improvement in Structured Sequencs Compared to Random','FontSize' , 28)
            case 'ActualvsfitHorz'
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(2,length(dayz),d)
                    hold on
                    for sn = 1
                        h1 = plotshade(coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)',err_pred{sn+1}(:,d)','transp' , .5 , 'patchcolor' , .5*colz{d,sn+1} , 'linecolor' , .5*colz{d,sn+1} , 'linewidth' , 3);
                        plot(coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , .5*colz{d,sn+1},'MarkerFaceColor',.5*colz{d,sn+1});
                        
                        h1 = plotshade(coo{sn+1}(:,d)',Plot{sn+1}(:,d)',err{sn+1}(:,d)','transp' , .5 , 'patchcolor' , colz{d,sn+1} , 'linecolor' , colz{d,sn+1} , 'linewidth' , 3);
                        plot(coo{sn+1}(:,d)',Plot{sn+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,sn+1},'MarkerFaceColor',colz{d,sn+1});
                     end
                    set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [2000 8000],'YTick' ,...
                        [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                    ylabel('Sec' )
                    xlabel('Viewing Horizon' )
                    title(['Chunked vs. fitted - Day ' , num2str(dayz{d})])
                end
                legend({'Chunked' , 'Fitted Chunked'})
                
                for d = 1:length(dayz)
                    subplot(2,length(dayz),length(dayz)+d)
                    hold on
                    for sn = 0
                        h1 = plotshade(coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)',err_pred{sn+1}(:,d)','transp' , .5 , 'patchcolor' , .5*colz{d,sn+1} , 'linecolor' , .5*colz{d,sn+1} , 'linewidth' , 3);
                        plot(coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , .5*colz{d,sn+1},'MarkerFaceColor',.5*colz{d,sn+1});
                        
                        h1 = plotshade(coo{sn+1}(:,d)',Plot{sn+1}(:,d)',err{sn+1}(:,d)','transp' , .5 , 'patchcolor' , colz{d,sn+1} , 'linecolor' , colz{d,sn+1} , 'linewidth' , 3);
                        plot(coo{sn+1}(:,d)',Plot{sn+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,sn+1},'MarkerFaceColor',colz{d,sn+1});
                    end
                    set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [2000 8000],'YTick' ,...
                        [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                    ylabel('Sec' )
                    xlabel('Viewing Horizon' )
                    title(['Random vs. fitted Random- Day ' , num2str(dayz{d})])
                end
            case 'Actual&fitSeqType'
                figure('color' , 'white')
                subplot(221)
                sn = 1;
                for d = 1:length(dayz)
                    errorbar(coo{sn+1}(:,d)',Plot{sn+1}(:,d)',err{sn+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
                hold on
                xlabel('Horizon')
                ylabel('msec')
                title('Chunked MT vs. Horizon')
                grid on
                
                subplot(222)
                sn = 0;
                for d = 1:length(dayz)
                    errorbar(coo{sn+1}(:,d)',Plot{sn+1}(:,d)',err{sn+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
                xlabel('Horizon')
                ylabel('msec')
                title('Random MT vs. Horizon')
                legend({'Day1' , 'Day2' ,'Day3','Day4','Day5'}, 'Box' , 'off')
                grid on
                
                subplot(223)
                sn = 1;
                for d = 1:length(dayz)
                    errorbar(coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)',err_pred{sn+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
                hold on
                xlabel('Horizon')
                ylabel('msec')
                title('Fitted Chunked MT vs. Horizon')
                grid on
                
                subplot(224)
                sn = 0;
                for d = 1:length(dayz)
                    errorbar(coo_pred{sn+1}(:,d)',plot_pred{sn+1}(:,d)',err_pred{sn+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20 , 'GridAlpha' , 1, 'Box' , 'off')
                xlabel('Horizon')
                title('Fitted Random MT vs. Horizon')
                grid on
            case 'plotCoef'
                h0 = figure;
               
                for sn = 0:1
                    [coo_b1{sn+1},plot_b1{sn+1},err_b1{sn+1}] = lineplot([coefs.Day] , coefs.b1 , 'subset' , ismember(coefs.seqNumb , [sn]));
                    [coo_b2{sn+1},plot_b2{sn+1},err_b2{sn+1}] = lineplot([coefs.Day] , coefs.b2 , 'subset' , ismember(coefs.seqNumb , [sn]));
                    [coo_b3{sn+1},plot_b3{sn+1},err_b3{sn+1}] = lineplot([coefs.Day] , coefs.b3 , 'subset' , ismember(coefs.seqNumb , [sn]));
                    [coo_invb3{sn+1},plot_invb3{sn+1},err_invb3{sn+1}] = lineplot([coefs.Day] , (coefs.b3).^-1 , 'subset' , ismember(coefs.seqNumb , [sn]));
                end
                close(h0)
                figure('color' , 'white')
                subplot(131)
                hold on
                for sn = 0:1
                    h1 = plotshade(coo_b1{sn+1}',plot_b1{sn+1},err_b1{sn+1},'transp' , .5 , 'patchcolor' , colz{3,sn+1} , 'linecolor' , colz{3,sn+1} , 'linewidth' , 3);
                    plot(coo_b1{sn+1},plot_b1{sn+1}, 'o' , 'MarkerSize' , 10 , 'color' , colz{3,sn+1},'MarkerFaceColor',colz{3,sn+1});
                end
                set(gca,'FontSize' , 20 ,'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YGrid' , 'on');
                hold on
                xlabel('Training session')
                title('b1 Coefficient')
                
                subplot(132)
                hold on
                for sn = 0:1
                    h1 = plotshade(coo_b2{sn+1}',plot_b2{sn+1},err_b2{sn+1},'transp' , .5 , 'patchcolor' , colz{3,sn+1} , 'linecolor' , colz{3,sn+1} , 'linewidth' , 3);
                    plot(coo_b2{sn+1},plot_b2{sn+1}, 'o' , 'MarkerSize' , 10 , 'color' , colz{3,sn+1},'MarkerFaceColor',colz{3,sn+1});
                end
                set(gca,'FontSize' , 20 ,'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YGrid' , 'on');
                hold on
                xlabel('Training session')
                title('b2 Coefficient')
                
                subplot(133)
                hold on
                for sn = 0:1
                    h1 = plotshade(coo_b3{sn+1}',plot_b3{sn+1},err_b3{sn+1},'transp' , .5 , 'patchcolor' , colz{3,sn+1} , 'linecolor' , colz{3,sn+1} , 'linewidth' , 3);
                    plot(coo_b3{sn+1},plot_b3{sn+1}, 'o' , 'MarkerSize' , 10 , 'color' , colz{3,sn+1},'MarkerFaceColor',colz{3,sn+1});
                end
                set(gca,'FontSize' , 20 ,'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YGrid' , 'on');
                hold on
                xlabel('Training session')
                title('b3 Coefficient')
            case 'testCoef'
                x = [1:10];
                b1 = 6;
                b2 = 10;
                
                figure('color' , 'white')
                subplot(121)
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
                subplot(122)
                for b2 = 0:20000:10^5
                    plot(x , b1 + (b2 - b1)*exp(-(x-1)/b3) ,'LineWidth' , 3)
                    hold on
                end
                legend({'b2  = 0' , 'b2  = 20000 ' ,'b2  = 40000' ,'b2  = 60000' ,'b2  = 80000' ,'b2  = 100000'}, 'Box' , 'off')
                title(['b1 = ' ,num2str(b1)  , ',     b3 = ' ,num2str(b3)])
                set(gca , 'FontSize' , 20 , 'Box' , 'off', 'GridAlpha' , 1)
                grid on
        end    
    case 'IPI_asymptote'
        structNumb = [1 2];
        out = [];
        plotcoef = 0;
        %         plotfcn = input('nanmean or nanmean?' , 's');
        %% IPIs vs horizon
        % this is the output of the case: 'transitions_All' that is saved to disc
        
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 , structNumb])  & ~Dall.isError);
        ANA.seqNumb(ANA.seqNumb==2) = 1;
        for tn = 1:length(ANA.TN)
            n = (ANA.AllPressIdx(tn , sum(~isnan(ANA.AllPressIdx(tn , :))))  - ANA.AllPressIdx(tn , 1)) / 1000;
            nIdx(tn , :) = (ANA.AllPressIdx(tn , :) - ANA.AllPressIdx(tn , 1))/n;
            ANA.IPI_norm(tn , :) = diff(nIdx(tn ,:) , 1 , 2);
        end
        for tn  = 1:length(ANA.TN)
            ANA.IPI_Horizon(tn , :) = ANA.Horizon(tn)*ones(1,13);
            ANA.IPI_SN(tn , :) = ANA.SN(tn)*ones(1,13);
            ANA.IPI_Day(tn , :) = ANA.Day(tn)*ones(1,13);
            ANA.IPI_prsnumb(tn , :) = [1 :13];
            ANA.IPI_seqNumb(tn , :) = ANA.seqNumb(tn)*ones(1,13);
        end
        IPItable.IPI = reshape(ANA.IPI , numel(ANA.IPI) , 1);
        IPItable.ChunkBndry = reshape(ANA.IPIarrangement , numel(ANA.IPI) , 1);
        IPItable.Horizon = reshape(ANA.IPI_Horizon , numel(ANA.IPI) , 1);
        IPItable.SN  = reshape(ANA.IPI_SN , numel(ANA.IPI) , 1);
        IPItable.Day = reshape(ANA.IPI_Day , numel(ANA.IPI) , 1);
        IPItable.prsnumb = reshape(ANA.IPI_prsnumb , numel(ANA.IPI) , 1);
        IPItable.seqNumb = reshape(ANA.IPI_seqNumb , numel(ANA.IPI) , 1);
        
        
        IPItable = getrow(IPItable , ismember(IPItable.prsnumb , [4:10]));
        IPItable.IPI_pred = zeros(size(IPItable.IPI));
        IPItable  = tapply(IPItable , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI' , 'nanmedian(x)'});
        for d = 1:length(dayz)
            IPItable.Day(ismember(IPItable.Day , dayz{d})) = d;
        end
        coefs = [];
        IPIs = [];
        for subjnum = 1:length(subj_name)-1
            for chp = 0:2
                for d = 1:length(dayz)
                    [d subjnum]
                    id = ismember(IPItable.SN , subjnum) & ismember(IPItable.Day , d) & ismember(IPItable.ChunkBndry , chp);
                    IPIch = getrow(IPItable , id);
                    exp_model1 = @(b,x) b(1) + (b(2) - b(1))*exp(-(x-1)/b(3)); % Model Function
                    %                 exp_model1 = @(b,x) b(1)*exp(-(x-1)/b(2)); % Model Function
                    x = IPIch.Horizon';
                    yx = IPIch.IPI';                                    % this would be a typical MT vs Horizon vector: [5422 3548 2704 2581 2446 2592 2418 2528 2500]
                    OLS = @(b) sum((exp_model1(b,x) - yx).^2);                % Ordinary Least Squares cost function
                    opts = optimset('MaxIter', MaxIter,'TolFun',1e-5);
                    [B1 Fval] = fminsearch(OLS,[3500 7500  1], opts);        % Use ?fminsearch? to minimise the ?OLS? function
                    IPIch.IPI_pred = exp_model1(B1,x)';
                    B.b1 = B1(1);
                    B.b2 = B1(2);
                    B.b3 = B1(3);
                    B.SN = subjnum;
                    B.Day = d;
                    B.ChunkBndry = chp;
                    B.IPI_fitted(1,:) = exp_model1(B1,unique(IPItable.Horizon))';
                    temp = tapply(IPIch , {'Horizon' } , {'IPI' , 'nanmedian(x)'});
                    B.Median_IPI(1,:) = temp.IPI;
                    B.Horizon(1,:) = unique(IPItable.Horizon);
                    coefs = addstruct(coefs , B);
                    IPIs = addstruct(IPIs , IPIch);
                    
                end
            end
            if plotcoef
                IPISN = getrow(IPIs , IPIs.SN == subjnum);
                f1 = figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(2,length(dayz),d)
                    hold on
                    for chp = 0
                        [xx , pp0 ,ee0] = lineplot([IPISN.Horizon] , IPISN.IPI , 'subset' , ismember(IPISN.ChunkBndry , chp)  & ismember(IPISN.Day ,d),...
                            'linecolor' , colIPI{d,chp+1} , 'plotfcn' , 'nanmean');
                    end
                    for chp = 0
                        [xx , pp1 ,ee1]=lineplot([IPISN.Horizon] , IPISN.IPI_pred , 'subset' , ismember(IPISN.ChunkBndry , chp)  & ismember(IPISN.Day ,d),...
                            'linecolor' , colIPI{d,chp+1}, 'plotfcn' , 'nanmean');
                    end
                    xlabel('Viewing Horizon' )
                    set(gca , 'YLim' , [min([pp0 pp1]-[ee0 ee1]) max([pp0 pp1]+[ee0 ee1])])
                end
                hold off
                
                coefSN = getrow(coefs , coefs.SN == subjnum);
                d= d+1;
                subplot(2,length(dayz),d)
                hold on
                for chp = 0:2
                    coefSNch = getrow(coefSN , ismember(coefSN.ChunkBndry , chp));
                    plot(coefSNch.Day , coefSNch.b1 , 'o-','color' , colIPI{4,chp+1} )
                end
                title('b1')
                hold off
                
                d= d+1;
                subplot(2,length(dayz),d)
                hold on
                for chp = 0:2
                    coefSNch = getrow(coefSN , ismember(coefSN.ChunkBndry , chp));
                    plot(coefSNch.Day , coefSNch.b2 , 'o-','color' , colIPI{4,chp+1} )
                end
                title('b2')
                hold off
                
                d= d+1;
                subplot(2,length(dayz),d)
                hold on
                for chp = 0:2
                    coefSNch = getrow(coefSN , ismember(coefSN.ChunkBndry , chp));
                    plot(coefSNch.Day , coefSNch.b3 , 'o-','color' , colIPI{4,chp+1} )
                end
                title('b3')
                hold off
                keyboard
                close(f1)
            end
        end
        IPIs_Pred = tapply(IPIs , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI_pred' , 'nanmean(x)'} , {'IPI' , 'nanmean(x)'});
        IPIs_Pred = normData(IPIs_Pred , {'IPI' , 'IPI_pred'});
        h0 = figure;
        for d = 1:length(dayz)
            for subjnum = 1:length(subj_name)-1
                for chp = 0:2
                    [coo{chp+1}(:,d),Plot{chp+1}(:,d),err{chp+1}(:,d)] = lineplot([IPIs_Pred.Horizon] , IPIs_Pred.normIPI , 'subset' , ismember(IPIs_Pred.ChunkBndry , chp)  & ismember(IPIs_Pred.Day ,d));
                    [coo_pred{chp+1}(:,d),plot_pred{chp+1}(:,d),err_pred{chp+1}(:,d)] = lineplot([IPIs_Pred.Horizon] , IPIs_Pred.normIPI_pred , 'subset' , ismember(IPIs_Pred.ChunkBndry , chp) & ismember(IPIs_Pred.Day , d));
                end
            end
        end
        close(h0)
        %%
        switch nowWhat
            case 'Actual&fitHorz'
                figure('color' , 'white')
                hz = length(unique(IPIs.Horizon));
                if hz == 9
                    hz = 13;
                end
                subplot(211)
                hold on
                for d = 1:length(dayz)
                    xtick = [];
                    for chp = 0:2
                        h1 = plotshade(chp*(hz+1)+coo{chp+1}(:,d)',Plot{chp+1}(:,d)',err{chp+1}(:,d)','transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                        plot(chp*(hz+1)+coo{chp+1}(:,d)',Plot{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                        xtick = [xtick ; chp*(hz+1)+coo{chp+1}(:,d)];
                    end
                    
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7' '8' '13'} ,1,chp+1), ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 max(xtick)], 'YLim' , [210 600],'YTick' ,...
                    [300 400 500 ] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on');
                ylabel('Sec' )
                xlabel('Viewing Horizon' )
                title(['IPIs - Day(s) ' , num2str(dayz{d})])
                
                
                subplot(212)
                hold on
                for d = 1:length(dayz)
                    xtick = [];
                    for chp = 0:2
                        h1 = plotshade(chp*(hz+1)+coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)','transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                        plot(chp*(hz+1)+coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                        xtick = [xtick ; chp*(hz+1)+coo{chp+1}(:,d)];
                    end
                    
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7' '8' '13'} ,1,chp+1), ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 max(xtick)], 'YLim' , [210 600],'YTick' ,...
                    [300 400 500 ] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on');
                ylabel('Sec' )
                xlabel('Viewing Horizon' )
                title(['Fitted IPIs - Day(s) ' , num2str(dayz{d})])
            case 'Actual&fitDayz'
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                h1 = figure;
                for d = 1:length(dayz)
                    for chp = 0:2
                        ANA = getrow(IPIs , ismember(IPIs.ChunkBndry , chp)  & ismember(IPIs.Day , d));
                        ANA_pred = getrow(IPIs_Pred , ismember(IPIs_Pred.ChunkBndry , chp)  & ismember(IPIs_Pred.Day , d));
                        ANA.Horizon(ANA.Horizon>6) = 6;
                        ANA_pred.Horizon(ANA_pred.Horizon>6) = 6;
                        [coo_red{chp+1}(:,d),Plot_red{chp+1}(:,d),err_red{chp+1}(:,d)] = lineplot([ANA.Horizon] , ANA.normIPI );
                        [coo_pred_red{chp+1}(:,d),plot_pred_red{chp+1}(:,d),err_pred_red{chp+1}(:,d)] = lineplot([ANA_pred.Horizon] , ANA_pred.normIPI_pred );
                    end
                end
                close(h1)
                scol = [200 200 200]/255;
                ecol = [0, 0 0]/255;
                for rgb = 1:3
                    horzcolor(rgb , :) = linspace(scol(rgb),ecol(rgb) , length(unique(ANA.Horizon)));
                end
                %                 horzcolor = transpose([107, 91, 149;254, 178, 54;214, 65, 97;255, 123, 37;128, 206, 214;241, 137, 115]/255);
                
                figure('color' , 'white')
%                 subplot(211)
                hold on
                xtick = [];
                for i = [1:length(unique(ANA.Horizon))]
                    for chp = 0:2
                        plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],Plot_red{chp+1}(i,:),err_red{chp+1}(i,:),'transp' , .6 , 'patchcolor' , colIPI{6-i+1,chp+1} , 'linecolor' , colIPI{6-i+1,chp+1} , 'linewidth' , 3);
                        plot((i-1)*(length(dayz)+1)+[1:length(dayz)],Plot_red{chp+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colIPI{6-i+1,chp+1},'MarkerFaceColor',colIPI{6-i+1,chp+1} , 'LineWidth' , 3);
                    end
                    xtick = [xtick (i-1)*(length(dayz)+1)+[1:length(dayz)]];
                    line([i*(length(dayz)+1) i*(length(dayz)+1)] , [210 600] , 'LineWidth' , 3 , 'LineStyle' , ':' , 'color' , [.8 .8 .8])
                end
                set(gca,'FontSize' , 18 , 'XTick' , xtick , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 i*(length(dayz)+1)], 'YLim' , [210 600],'YTick' ,...
                    [300 400 500] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on',...
                    'XTickLabels' , repmat([1:5] , 1 ,length(unique(ANA.Horizon))));
                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Training Session','FontSize' , 20)
                title('Actual','FontSize' , 24)
                
                figure('color' , 'white')
                hold on
                xtick = [];
                for i = [1:length(unique(ANA.Horizon))]
                    for chp = 0:2
                        plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{chp+1}(i,:),err_pred_red{chp+1}(i,:),'transp' , .6 , 'patchcolor' , colIPI{6-i+1,chp+1} , 'linecolor' , colIPI{6-i+1,chp+1} , 'linewidth' , 3);
                        plot((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{chp+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colIPI{6-i+1,chp+1},'MarkerFaceColor',colIPI{6-i+1,chp+1} , 'LineWidth' , 3);
                    end
                    xtick = [xtick (i-1)*(length(dayz)+1)+[1:length(dayz)]];
                    line([i*(length(dayz)+1) i*(length(dayz)+1)] , [210 600] , 'LineWidth' , 3 , 'LineStyle' , ':' , 'color' , [.8 .8 .8])
                end
                set(gca,'FontSize' , 18 , 'XTick' , xtick , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 i*(length(dayz)+1)], 'YLim' , [210 600],'YTick' ,...
                    [300 400 500] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on',...
                    'XTickLabels' , repmat([1:5] , 1 ,length(unique(ANA.Horizon))));
                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Training Session','FontSize' , 20)
                title('Fittted','FontSize' , 24)
            case 'Actual&fit%ChangeDayzTotalLearning'
                Hz = {[1] [2] [3] , [4] [5] [6] [7:9]};
                ANA = IPIs;
                ANA.Horizon(ANA.Horizon>Hz{end}(1)) = Hz{end}(1);
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI' , 'nanmean(x)'},{'IPI_pred' , 'nanmean(x)'});
                ANA.percChangeIPI = zeros(size(ANA.IPI));
                ANA.percChangeIPI_pred = zeros(size(ANA.IPI_pred));
                Daybenefit = [];
                for d = 2:length(dayz)
                    for chp = 0:2
                        Db1= getrow(ANA , ANA.Day == 1 & ANA.ChunkBndry==chp);
                        Db_d = getrow(ANA , ANA.Day == d & ANA.ChunkBndry==chp);
                        Db_d.percChangeIPI = 100*abs((Db_d.IPI - Db1.IPI)./Db1.IPI);
                        Db_d.percChangeIPI_pred = 100*abs((Db_d.IPI_pred - Db1.IPI_pred)./Db1.IPI_pred);
                        Daybenefit = addstruct(Daybenefit , Db_d);
                    end
                end
      
                Daybenefit  = tapply(Daybenefit , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'percChangeIPI' , 'nanmean(x)'},{'percChangeIPI_pred' ,  'nanmean(x)'});
                Daybenefit = normData(Daybenefit , {'percChangeIPI' , 'percChangeIPI_pred'});
                figure('color' , 'white')
                colorz = colIPI(d,1:4);
                lineplot([Daybenefit.Day Daybenefit.Horizon ] , Daybenefit.percChangeIPI , 'plotfcn' , 'nanmean',...
                    'split', Daybenefit.ChunkBndry , 'linecolor' , colorz,...
                    'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                    'linewidth' , 2 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'Random'  , 'Within' , 'Between'} );
                set(gca,'FontSize' , 18 , 'XTickLabel' , repmat({'1' '2' '3' '4' '5' '6' '7-13'} ,1, length(dayz)-1),...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [5 30],'YGrid' , 'on','XTickLabelRotation' , 30);
                xlabel('Veiwing window size (W)' )
                ylabel('% improvement')
                title('Reduction in Inter-Press Intervals compared to Day 1' ,'FontSize' , 24)
            case 'Actual&fit%ChangeDay2Day'
                Hz = {[1] [2] [3] , [4] [5] [6] [7:9]};
                ANA = IPIs;
                ANA.Horizon(ANA.Horizon>Hz{end}(1)) = Hz{end}(1);
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI' , 'nanmean(x)'},{'IPI_pred' , 'nanmean(x)'});
                ANA.percChangeIPI = zeros(length(ANA.IPI) , length(dayz)-1);
                ANA.percChangeIPI_pred = zeros(length(ANA.IPI) , length(dayz)-1);
                
                Daybenefit = [];
                for chp = 0:2
                    for d = 1:length(dayz)-1
                        Db1= getrow(ANA , ANA.Day == d & ANA.ChunkBndry==chp);
                        Db = getrow(ANA , ANA.Day == d+1 & ANA.ChunkBndry==chp);
                        Db.percChangeIPI(:,d) = 100*abs((Db1.IPI - Db.IPI)./Db1.IPI);
                        Db.percChangeIPI_pred(:,d) = 100*abs((Db1.IPI_pred - Db.IPI_pred)./Db1.IPI_pred);
                        Daybenefit = addstruct(Daybenefit , Db);
                    end
                end
                if poolDays
                    daylab = {'Session 1 to 2,3' , 'Sessions 2,3 to 4,5'};
                else
                    daylab = {'Session 1 to 2' , 'Session 2 to 3' , 'Session 3 to 4' , 'Session 4 to 5'};
                end
                Daybenefit  = tapply(Daybenefit , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'percChangeIPI' , 'nanmean(x)'},{'percChangeIPI_pred' ,  'nanmean(x)'});
                Daybenefit = normData(Daybenefit , {'percChangeIPI' , 'percChangeIPI_pred'});
                figure('color' , 'white')
                colorz = colIPI(d,1:3);
                lineplot([ Daybenefit.Horizon Daybenefit.Day] , Daybenefit.percChangeIPI , 'plotfcn' , 'nanmean',...
                    'split', Daybenefit.ChunkBndry , 'linecolor' , colorz,...
                    'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colorz,...
                    'linewidth' , 2 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colorz,...
                    'markersize' , 5, 'markercolor' , colorz , 'leg' , {'Random'  , 'Within' , 'Between'} );
                set(gca,'FontSize' , 18 , 'XTickLabels' , repmat(daylab , 1,length(Hz)), ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [0 6],'YGrid' , 'on','XTickLabelRotation' , 30);
                xlabel('Day compared to day 1' )
                ylabel('% improvement')
                title(['Day-to-day Improvement within IPI type (4-10 included) , W = 1, 2, 3, 4, 5, 6, 7:13'])
            case 'ActualvsfitHorz'
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(2,length(dayz),d)
                    hold on
                    for chp = 0:2
                        h1 = plotshade(coo{chp+1}(:,d)',Plot{chp+1}(:,d)',err{chp+1}(:,d)','transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                        plot(coo{chp+1}(:,d)',Plot{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                    end
                    title('Actual')
                    set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [200 600],'YTick' ,...
                        [300 400 500] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on');
                    subplot(2,length(dayz),length(dayz)+d)
                    hold on
                    for chp = 0:2
                        h1 = plotshade(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)','transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                        plot(coo{chp+1}(:,d)',plot_pred{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                    end
                    title('Fitted')
                    set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [200 600],'YTick' ,...
                        [300 400 500] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on');
                    ylabel('Sec' )
                    xlabel('Viewing Horizon' )
                end
            case 'plotCoef'
                h0 = figure;
                if poolDays
                    IPIs.Day(IPIs.Day==3) = 2;
                    IPIs.Day(ismember(IPIs.Day , [4 5])) = 3;
                end
                
                IPIs  = tapply(IPIs , {'Day' ,'SN' , 'ChunkBndry'} , {'b1' , 'nanmean(x)'},{'b2' ,  'nanmean(x)'},{'b3' ,  'nanmean(x)'});
                IPIs = normData(IPIs , {'b1' , 'b2' , 'b3'});
                for chp = 0:2
                    [coo_b1{chp+1},plot_b1{chp+1},err_b1{chp+1}] = lineplot([coefs.Day] , coefs.b1 , 'subset' , ismember(coefs.ChunkBndry , [chp]));
                    [coo_b2{chp+1},plot_b2{chp+1},err_b2{chp+1}] = lineplot([coefs.Day] , coefs.b2 , 'subset' , ismember(coefs.ChunkBndry , [chp]));
                    [coo_b3{chp+1},plot_b3{chp+1},err_b3{chp+1}] = lineplot([coefs.Day] , coefs.b3 , 'subset' , ismember(coefs.ChunkBndry , [chp]));
                    [coo_invb3{chp+1},plot_invb3{chp+1},err_invb3{chp+1}] = lineplot([coefs.Day] , (coefs.b3).^-1 , 'subset' , ismember(coefs.ChunkBndry , [chp]));
                end
                close(h0)
                figure('color' , 'white')
                subplot(121)
                hold on
                for chp = 0:2
                    h1 = plotshade(coo_b1{chp+1}',plot_b1{chp+1},err_b1{chp+1},'transp' , .4 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{4,chp+1} , 'linewidth' , 3);
                    plot(coo_b1{chp+1},plot_b1{chp+1}, 'o' , 'MarkerSize' , 10 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{4,chp+1});
                end
                set(gca,'FontSize' , 20 ,'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YGrid' , 'on');
                hold on
                xlabel('Training session')
                title('b1 Coefficient')
                
                subplot(122)
                hold on
                for chp = 0:2
                    h1 = plotshade(coo_b2{chp+1}',plot_b2{chp+1},err_b2{chp+1},'transp' , .4 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{4,chp+1} , 'linewidth' , 3);
                    plot(coo_b2{chp+1},plot_b2{chp+1}, 'o' , 'MarkerSize' , 10 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{4,chp+1});
                end
                set(gca,'FontSize' , 20 ,'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YGrid' , 'on');
                hold on
                xlabel('Training session')
                title('b2 Coefficient')
                
                figure('color' , 'white')
                hold on
                for chp = 0:2
                    h1 = plotshade(coo_b3{chp+1}',plot_b3{chp+1},err_b3{chp+1},'transp' , .4 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{4,chp+1} , 'linewidth' , 3);
                    plot(coo_b3{chp+1},plot_b3{chp+1}, 'o' , 'MarkerSize' , 10 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{4,chp+1});
                end
                set(gca,'FontSize' , 20 ,'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], ...
                    'XTick' , [1:5] , 'YGrid' , 'on');
                hold on
                xlabel('Training Sessions')
                title('Exponential Decay Fit Time Constant')
            case 'testCoef'
                x = [1:10];
                b1 = 6;
                b2 = 10;
                
                figure('color' , 'white')
                subplot(121)
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
                subplot(122)
                for b2 = 0:20000:10^5
                    plot(x , b1 + (b2 - b1)*exp(-(x-1)/b3) ,'LineWidth' , 3)
                    hold on
                end
                legend({'b2  = 0' , 'b2  = 20000 ' ,'b2  = 40000' ,'b2  = 60000' ,'b2  = 80000' ,'b2  = 100000'}, 'Box' , 'off')
                title(['b1 = ' ,num2str(b1)  , ',     b3 = ' ,num2str(b3)])
                set(gca , 'FontSize' , 20 , 'Box' , 'off', 'GridAlpha' , 1)
                grid on
        end
    case 'Eye'
        calc = 0;
        if isSymmetric 
            filename = 'se2_eyeInfo.mat';
        else
            filename = 'se2_eyeInfo_asym.mat';
        end
        if calc
            eyeinfo.PB         = [];   % preview benefit
            eyeinfo.CB         = [];   % chunk boundry
            eyeinfo.Horizon    = [];
            eyeinfo.sn         = [];
            eyeinfo.Day        = [];
            eyeinfo.BN         = [];
            eyeinfo.TN         = [];
            eyeinfo.sacPerSec  = [];
            eyeinfo.sacDur     = [];
            eyeinfo.sacPeakVel = [];
            eyeinfo.sacAmp     = [];
            eyeinfo.seqNumb    = [];
            eyeinfo.DigFixDur  = [];
            eyeinfo.prsnumb    = [];
            eyeinfo.prstimepos = [];
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:2]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & cellfun(@length , Dall.xEyePosDigit)>1);
           
            for tn  = 1:length(ANA.TN)
                if ismember(ANA.seqNumb(tn) , [1:2])
                    ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                    a = find(ANA.ChunkBndry(tn , :));
                    ANA.ChunkBndry(tn , a(2:end)-1) = 3;
                    ANA.ChunkBndry(tn ,end) = 3;
                    ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                else
                    ANA.ChunkBndry(tn , :) = zeros(size(ANA.ChnkArrang(tn , :)));
                end
                ANA.DigFixWeighted(tn , :) = zeros(1 ,14);
                window = 74;
                if isSymmetric
                    for p = 1:14
                        id = ANA.xEyePosDigit{tn , 1}<=p+.5 & ANA.xEyePosDigit{tn , 1}>p-.5;
                        if sum(id)
                            ANA.DigFixWeighted(tn , p) = mean(abs(ANA.xEyePosDigit{tn , 1}(id) - p))*(sum(id)/500)*1000;
                        else
                            ANA.DigFixWeighted(tn , p) = 0;
                        end
                    end
                else
                    for p = 1:14
                        id = ANA.xEyePosDigit{tn , 1}<=p+.3 & ANA.xEyePosDigit{tn , 1}>p-.7;
                        if sum(id)
                            ANA.DigFixWeighted(tn , p) = mean(abs(ANA.xEyePosDigit{tn , 1}(id) - p))*(sum(id)/500)*1000;
                        else
                            ANA.DigFixWeighted(tn , p) = 0;
                        end
                    end
                    
                end
                
                for p = 1:14
                    id = [ANA.AllPressIdx(tn , p) - window :ANA.AllPressIdx(tn , p) + 1];
                    if id(1) > length(ANA.xEyePosDigit{tn}) | sum(id<0)>0
                        ANA.EyePressTimePos(tn , p) = NaN;
                    elseif id(end)>length(ANA.xEyePosDigit{tn})
                        ANA.EyePressTimePos(tn , p) = nanmedian(ANA.xEyePosDigit{tn}(id(1):end));
                    else
                        ANA.EyePressTimePos(tn , p) = nanmedian(ANA.xEyePosDigit{tn}(id));
                    end
                end
                perv_Ben           = [1:14] - ANA.EyePressTimePos(tn , :);
                goodid             = ~(abs(perv_Ben)>=3.5);
                prsnumb            = find(goodid);
                count              = sum(goodid);
                eyeinfo.PB         = [eyeinfo.PB ;perv_Ben(goodid)'];
                eyeinfo.prsnumb    = [eyeinfo.prsnumb ;find(goodid')];
                eyeinfo.CB         = [eyeinfo.CB ;ANA.ChunkBndry(tn ,goodid)'];
                eyeinfo.Horizon    = [eyeinfo.Horizon ; ANA.Horizon(tn)*ones(count , 1)];
                eyeinfo.sn         = [eyeinfo.sn ; ANA.SN(tn)*ones(count , 1)];
                eyeinfo.Day        = [eyeinfo.Day ; ANA.Day(tn)*ones(count , 1)];
                eyeinfo.BN         = [eyeinfo.BN ; ANA.BN(tn)*ones(count , 1)];
                eyeinfo.TN         = [eyeinfo.TN ; ANA.TN(tn)*ones(count , 1)];
                eyeinfo.sacPerSec  = [eyeinfo.sacPerSec ; ANA.SaccPerSec(tn)*ones(count , 1)];
                eyeinfo.sacDur     = [eyeinfo.sacDur ; mean(ANA.SaccDuration{tn})*ones(count , 1)];
                eyeinfo.sacPeakVel = [eyeinfo.sacPeakVel ; mean(ANA.SaccPeakVel{tn})*ones(count , 1)];
                eyeinfo.sacAmp     = [eyeinfo.sacAmp ; mean(ANA.SaccAmplitude{tn})*ones(count , 1)];
                eyeinfo.seqNumb    = [eyeinfo.seqNumb ; ANA.seqNumb(tn)*ones(count , 1)];
                eyeinfo.DigFixDur  = [eyeinfo.DigFixDur ;ANA.DigFixWeighted(tn ,goodid)'];
                eyeinfo.prstimepos = [eyeinfo.prstimepos ;ANA.EyePressTimePos(tn ,goodid)'];
            end

            save([baseDir , '/' , filename] , 'eyeinfo','-v7.3')
        else
            load([baseDir , '/', filename])
        end
        if poolDays
            eyeinfo = getrow(eyeinfo , ismember(eyeinfo.Day , [1 3 4 5]));
            eyeinfo.Day(eyeinfo.Day==3) = 2;
            eyeinfo.Day(ismember(eyeinfo.Day , [4 5])) = 3;
        end
        eyeinfo = getrow(eyeinfo , ismember(eyeinfo.sn , subjnum));
        out = [];
        switch nowWhat
            case 'sacDurSplitDay'
                Ho  = unique(eyeinfo.Horizon);
                
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb'} , {'sacDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacDur'});
                figure('color' , 'white')
                subplot(121)
                lineplot([K.Horizon] , K.normsacDur , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 0) ,'split', K.Day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),1),...
                    'errorcolor' , colz(1:length(dayz),1) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),1),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),1),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),1));
                ylabel('Mean Saccade duration per trial [ms]' )
                xlabel('Viewing window Size' )
                title('Random Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [20  70],'YTick' , [30:10:70] ,...
                    'YGrid' , 'on');
                
                subplot(122)
                lineplot([K.Horizon] , K.normsacDur , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 1) ,'split', K.Day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),2),...
                    'errorcolor' , colz(1:length(dayz),2) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),2),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),2),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),2));
                ylabel('Mean Saccade duration per trial [ms]' )
                xlabel('Viewing window Size' )
                title('Structured Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [20  70],'YTick' , [30:10:70] ,...
                    'YGrid' , 'on','YColor' , 'none');
            case 'sacDurSplitseqType'
                Ho  = unique(eyeinfo.Horizon);
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb'} , {'sacDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacDur'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Horizon] , K.normsacDur , 'plotfcn' , 'nanmean',...
                        'split', K.seqNumb , 'subset' , ismember(K.Day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colz(d,1:2),...
                        'errorcolor' , colz(d,1:2) , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colz(d,1:2),...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colz(d,1:2),...
                        'markersize' , 10, 'markercolor' , colz(d,1:2));
                    ylabel('Mean Saccade duration per trial [ms]' )
                    xlabel('Viewing window Size' )
                    title('Random Sequences')
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [20  70],'YTick' , [30:10:70] ,...
                        'YGrid' , 'on');
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end  
            case 'sacAmpSplitDay'
                Ho  = unique(eyeinfo.Horizon);
                
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb'} , {'sacAmp' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacAmp'});
                figure('color' , 'white')
                subplot(121)
                lineplot([K.Horizon] , K.normsacAmp , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 0) ,'split', K.Day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),1),...
                    'errorcolor' , colz(1:length(dayz),1) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),1),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),1),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),1));
                ylabel('Mean Saccade amplitude per trial [deg]' )
                xlabel('Viewing window Size' )
                title('Random Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [-4 20],'YTick' , [0:5:20] ,...
                    'YGrid' , 'on');
                
                subplot(122)
                lineplot([K.Horizon] , K.normsacAmp , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 1) ,'split', K.Day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),2),...
                    'errorcolor' , colz(1:length(dayz),2) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),2),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),2),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),2));
                ylabel('Mean Saccade amplitude per trial [deg]' )
                xlabel('Viewing window Size' )
                title('Structured Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [-4 20],'YTick' , [0:5:20] ,...
                    'YGrid' , 'on','YColor' , 'none');
            case 'sacAmpSplitseqType'
                Ho  = unique(eyeinfo.Horizon);
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb'} , {'sacAmp' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacAmp'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Horizon] , K.normsacAmp , 'plotfcn' , 'nanmean',...
                        'split', K.seqNumb , 'subset' , ismember(K.Day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colz(d,1:2),...
                        'errorcolor' , colz(d,1:2) , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colz(d,1:2),...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colz(d,1:2),...
                        'markersize' , 10, 'markercolor' , colz(d,1:2));
                    ylabel('Mean Saccade amplitude per trial [deg]' )
                    xlabel('Viewing window Size' )
                    title('Random Sequences')
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [-4 20],'YTick' , [0:5:20] ,...
                        'YGrid' , 'on');
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end
            case 'sacFreqSplitDay'
                Ho  = unique(eyeinfo.Horizon);
                
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb'} , {'sacPerSec' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacPerSec'});
                figure('color' , 'white')
                subplot(121)
                lineplot([K.Horizon] , K.normsacPerSec , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 0) ,'split', K.Day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),1),...
                    'errorcolor' , colz(1:length(dayz),1) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),1),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),1),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),1));
                ylabel('Mean Saccade frequency [Hz]' )
                xlabel('Viewing window Size' )
                title('Random Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.9 2.3],'YTick' , [1:.2:2.2] ,...
                        'YGrid' , 'on');
                
                subplot(122)
                lineplot([K.Horizon] , K.normsacPerSec , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 1) ,'split', K.Day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),2),...
                    'errorcolor' , colz(1:length(dayz),2) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),2),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),2),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),2));
                ylabel('Mean Saccade frequency [Hz]' )
                xlabel('Viewing window Size' )
                title('Structured Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.9 2.3],'YTick' , [1:.2:2.2] ,...
                        'YGrid' , 'on','YColor' , 'none');
            case 'sacFreqSplitseqType'
                Ho  = unique(eyeinfo.Horizon);
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb'} , {'sacPerSec' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacPerSec'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Horizon] , K.normsacPerSec , 'plotfcn' , 'nanmean',...
                        'split', K.seqNumb , 'subset' , ismember(K.Day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colz(d,1:2),...
                        'errorcolor' , colz(d,1:2) , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colz(d,1:2),...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colz(d,1:2),...
                        'markersize' , 10, 'markercolor' , colz(d,1:2));
                    ylabel('Mean Saccade frequency [Hz]' )
                    xlabel('Viewing window Size' )
                    title('Saccade rate')
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.9 2.3],'YTick' , [1:.2:2.2] ,...
                        'YGrid' , 'on');
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end
            case 'FixDurSplitipitype'
                Ho  = unique(eyeinfo.Horizon);
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','CB'} , {'DigFixDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                
                K = normData(K , {'DigFixDur'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Horizon] , K.normDigFixDur , 'plotfcn' , 'nanmean',...
                        'split', K.CB , 'subset' , ismember(K.Day , dayz{d}) , 'leg' , {'Random' , 'First' , 'Middle' , 'Last'} , 'linecolor' , colIPI(d,1:4),...
                        'errorcolor' , colIPI(d,1:4) , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colIPI(d,1:4),...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colIPI(d,1:4),...
                        'markersize' , 10, 'markercolor' , colIPI(d,1:4));
                    ylabel('Mean digit fixation Duration [ms]' )
                    xlabel('Viewing window Size' )
                    title('Digit fixation duration')
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [40 140],'YTick' , [50:10:140] ,...
                        'YGrid' , 'on');
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end    
            case 'FixDurSplitwindow'
                Ho  = unique(eyeinfo.Horizon);
                scol = [200 200 200]/255;
                ecol = [0, 0 0]/255;
                for rgb = 1:3
                    temp(: , rgb) = linspace(scol(rgb),ecol(rgb) , length(dayz));
                end
                for d = 1:length(dayz)
                    horzcolor{d} = temp(d,:);
                end
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','CB'} , {'DigFixDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K = normData(K , {'DigFixDur'});
                figure('color' , 'white')
                for h = 1:length(Ho)
                    subplot(1,length(Ho) , h)
                    colorz  = horzcolor;
                    lineplot([K.CB] , K.normDigFixDur , 'plotfcn' , 'nanmean',...
                        'split', K.Day , 'subset' , ismember(K.Horizon , Ho(h)) , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colorz,...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colorz,...
                        'markersize' , 10, 'markercolor' , colorz , 'leg' , {'day1' , 'day2' , 'day3' , 'day4' , 'day5'});
                    ylabel('Mean digit fixation Duration [ms]' )
                    xlabel('IPI type' )
                    title(['W = ' , num2str(Ho(h))])
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [40 140],'YTick' , [50:10:140] ,...
                        'YGrid' , 'on');
                    if h>1
                        set(gca,'YColor' , 'none');
                    end
                end
            case 'previewSplitipitype'
                Ho  = unique(eyeinfo.Horizon);
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','CB'} , {'PB' , 'nanmedian'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                
                K = normData(K , {'PB'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Horizon] , -K.normPB , 'plotfcn' , 'nanmean',...
                        'split', K.CB , 'subset' , ismember(K.Day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colIPI(d,1:4),...
                        'errorcolor' , colIPI(d,1:4) , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colIPI(d,1:4),...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colIPI(d,1:4),...
                        'markersize' , 10, 'markercolor' , colIPI(d,1:4));
                    ylabel('Mean preview [digits]' )
                    xlabel('Viewing window Size' )
                    title(['Look-ahead - day ' , num2str(dayz{d})])
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.4 1.4],'YTick' , [.5:.1:1.3] ,...
                        'YGrid' , 'on');
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end    
            case 'previewSplitwindow'
                Ho  = unique(eyeinfo.Horizon);
                
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','CB'} , {'PB' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K = normData(K , {'PB'});
                scol = [200 200 200]/255;
                ecol = [0, 0 0]/255;
                for rgb = 1:3
                    temp(: , rgb) = linspace(scol(rgb),ecol(rgb) , length(dayz));
                end
                for d = 1:length(dayz)
                    horzcolor{d} = temp(d,:);
                end
                figure('color' , 'white')
                for h = 1:length(Ho)
                    subplot(1,length(Ho) , h)
                    colorz  = horzcolor;
                    lineplot([K.CB] , -K.normPB , 'plotfcn' , 'nanmean',...
                        'split', K.Day , 'subset' , ismember(K.Horizon , Ho(h)) & ismember(K.Day , [1 5]) , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colorz,...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colorz,...
                        'markersize' , 10, 'markercolor' , colorz);
                    ylabel('Mean preview [digits]' )
                    xlabel('IPI type' )
                    title(['Look-ahead  W = ' , num2str(Ho(h))])
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.4 1.4],'YTick' , [.5:.1:1.3] ,...
                        'YGrid' , 'on');
                    if h>1
                        set(gca,'YColor' , 'none');
                    end
                end
                
            case 'previewSplitDays'
                Ho  = unique(eyeinfo.Horizon);
                
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','CB'} , {'PB' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K = normData(K , {'PB'});
                K = getrow(K , ismember(K.Day , [1 , length(dayz)]));
                figure('color' , 'white')
                for cb = {[0] [1 2 3]}
                    subplot(1,2,cb{1}(1)+1)
                    colorz  = colz(1:length(dayz) , cb{1}(1)+1);
                    
                    lineplot([K.Horizon] , -K.normPB , 'plotfcn' , 'nanmean',...
                        'split', K.Day , 'subset' ,  ismember(K.CB , cb{1})  , 'linecolor' , colorz,...
                        'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colorz,...
                        'linewidth' , 3 , 'markertype' , {'o' , '>' , '<' , 'x' , 's'}  , 'markerfill' , colorz,...
                        'markersize' , 15, 'markercolor' , colorz , 'leg' , 'auto');
                    ylabel('Mean preview [digits]' )
                    xlabel('Viewing window size (W)' )
                    title(['Look-ahead  CB = ' , num2str(cb{1})])
%                     set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.4 1.4],'YTick' , [.5:.1:1.3] ,...
%                         'YGrid' , 'on');
                end
                
            case 'presspositionlook_ahead'
                K = eyeinfo;
                dayz = {[1] , [4 5]};
                K = getrow(K , ismember(K.seqNumb , [0]) & ~isnan(K.PB) );
                %  K.prsnumb(ismember(K.prsnumb ,[5:9])) = 5;
                 K.Horizon(K.Horizon>8) = 9;
                A = [];
                for dd = 1:length(dayz)
                    temp  = getrow(K , ismember(K.Day , dayz{dd}));
                    temp.Day(1:end)  = dd;
                    A = addstruct(A , temp);
                end
                K = A;
                K = tapply(K , {'Day' , 'Horizon' , 'sn', 'prsnumb' , 'seqNumb'} , {'PB' , 'nanmedian'} );
                K = normData(K , {'PB'});

                figure('color' , 'white')
                colorz  = colz([1,end] , unique(K.seqNumb)+1);
                lineplot([K.Horizon ] , -K.normPB , 'plotfcn' , 'nanmean',...
                    'split', K.Day , 'linecolor' , colorz,...
                    'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colorz,...
                    'linewidth' , 1 , 'markertype' , {'^', '*'}  , 'markerfill' , colorz,...
                    'markersize' , 5, 'markercolor' , colorz , 'leg' , 'auto', 'subset' , K.prsnumb==1);
                ylabel('Mean preview [digits]' )
                xlabel('Viewing window size (W)' )
                title(['Look-ahead  CB = ' , num2str(cb{1})])
                set(gca,'FontSize' , 7 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [-1 2],'XLim' , [.5 66.5],'YTick' , [-.5:.5:2] ,...
                    'XTickLabel' , repmat({'1','','','','','','','','13'} , 1 , length(unique(K.prsnumb))),'XTickLabelRotation' , 0)
            case 'startlookahead'
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','CB' , 'prsnumb'} , {'PB' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [1]));
                K = normData(K , {'PB'});
                K = getrow(K , ismember(K.Day , [1 , length(dayz)]));
                figure('color' , 'white')
                    for cb = {[0] [1 2 3]}
                        subplot(1,2,cb{1}(1)+1)
                        colorz  = colz(1:length(dayz) , cb{1}(1)+1);
                        
                        lineplot([K.prsnumb K.Horizon ] , -K.normPB , 'plotfcn' , 'nanmean',...
                            'split', K.Day , 'subset' ,  ismember(K.CB , cb{1}) , 'linecolor' , colorz,...
                            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colorz,...
                            'linewidth' , 3 , 'markertype' , {'+' , '*'}  , 'markerfill' , colorz,...
                            'markersize' , 15, 'markercolor' , colorz , 'leg' , 'auto');
                        ylabel('Mean preview [digits]' )
                        xlabel('Viewing window size (W)' )
                        title(['Look-ahead  CB = ' , num2str(cb{1})])
                        set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [.6 2],'YTick' , [.6:.2:2] ,...
                            'YGrid' , 'on');
                    end
            case 'EyePrsTimePos'
                Ho  = unique(eyeinfo.Horizon);
                scol = [200 200 200]/255;
                ecol = [0, 0 0]/255;
                for rgb = 1:3
                    temp(: , rgb) = linspace(scol(rgb),ecol(rgb) , length(dayz));
                end
                for d = 1:length(dayz)
                    horzcolor{d} = temp(d,:);
                end
                K = tapply(eyeinfo , {'Day' , 'Horizon' , 'sn','seqNumb' , 'prsnumb'} , {'prstimepos' , 'nanmean'});
                K = normData(K , {'prstimepos'});
                figure('color' , 'white')
                fcount = 1;
                for sqn = 0:2
                    for h = 1:length(Ho)
                        subplot(3,length(Ho) , fcount)
                        colorz  = horzcolor;
                        lineplot([K.prsnumb] , K.normprstimepos , 'plotfcn' , 'nanmean',...
                            'split', K.Day , 'subset' , ismember(K.Horizon , Ho(h)) & ismember(K.seqNumb , sqn) , 'linecolor' , colorz,...
                            'errorcolor' , colorz , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colorz,...
                            'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colorz,...
                            'markersize' , 10, 'markercolor' , colorz);
%                         ylabel('Eye position [digits]' )
%                         xlabel('Press position' )
%                         title(['W = ' , num2str(Ho(h)) , 'seqNum = ' , num2str(sqn) ])
                        axis square
                        set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [1 14],'YTick' , [1 : 14] ,...
                            'XLim' , [1 14],'XTick' , [1 : 14],'YGrid' , 'on' , 'XGrid' , 'on');
                        if h>1
                            set(gca,'YColor' , 'none');
                        end
                        fcount = fcount + 1;
                    end
                end
        end
end


