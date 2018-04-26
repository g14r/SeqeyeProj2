
function out  = se2_pubFigs(Dall , what, nowWhat , varargin)

subj_name = {'AT1' , 'CG1' , 'HB1' , 'JT1' , 'CB1' , 'YM1' , 'NL1' , 'SR1' , 'IB1' , 'MZ1' , 'DW1', 'RA1' ,'CC1', 'All'};
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
        otherwise
            error('Unknown option: %s',varargin{c});
    end
end

%%
prefix = 'se1_';
baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
% baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye2/analyze';          %iMac



if subjnum == length(subj_name)
    %     subjnum = [1 3:length(subj_name)-1];
    subjnum = 1:length(subj_name)-1;
end


days  = {1 ,2 ,3 ,4 ,5,[1:5] ,[2:5] [2:3] [4:5],[3:5]};
if poolDays
    dayz = {1 [2 3] [4 5]};
else
    dayz = {[1] [2] [3] [4] [5]};
end

clear tempcol
c1 = [255, 153, 179]/255; % Random Red Tones
ce = [153, 0, 51]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
end
for i = 1:length(tempcol)
    colz{i,1} = tempcol(: , i)';
end

clear tempcol
c1 = [153, 194, 255]/255; % structured blue tones
ce = [0, 0, 102]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
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
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
end
for i = 1:length(tempcol)
    colIPI{i,2} = tempcol(: , i)';
end
clear tempcol
c1 = [179, 255, 153]/255; % middle chunk green tones
ce = [19, 77, 0]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
end
for i = 1:length(tempcol)
    colIPI{i,3} = tempcol(: , i)';
end

clear tempcol
c1 = [255, 179, 255]/255; % last chunk fuscia tones
ce = [102, 0, 102]/255;
for rgb = 1:3
    tempcol(rgb , :) = linspace(c1(rgb),ce(rgb) , 6);
end
for i = 1:length(tempcol)
    colIPI{i,4} = tempcol(: , i)';
end




switch what
    
    case 'MT'
        ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError);
        ANA = getrow(Dall , Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError );
        ANA.seqNumb(ANA.seqNumb == 2) = 1;
        MT  = ANA;%tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb','BN'} , {'MT' , 'nanmean(x)'});
        
        MT = getrow(MT , MT.MT <= 9000 );
        
        for d = 1:length(dayz)
            MT.Day(ismember(MT.Day , dayz{d})) = d;
        end
        h1 = figure('color' , 'white');
        MT = tapply(MT , {'Horizon' , 'Day' , 'seqNumb' , 'SN'} , {'MT'});
        MT = normData(MT , {'MT'});
        for d=  1:length(dayz)
            hc = 1;
            [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(MT.Horizon,MT.normMT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 1 & ismember(MT.Day , d));
            hold on
            [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(MT.Horizon,MT.normMT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 0 & ismember(MT.Day , d));
        end
        close(h1);
        if poolDays
            sigSeq = [NaN 3 2];
            sigMT  = [1 4 5;4 5 4];
        else
            sigSeq = [NaN 3 3 2 2];
            sigMT  = [3 NaN NaN NaN 5;4 NaN NaN NaN 4];
        end
        
        switch nowWhat
            case 'RandvsStructCommpare'
                figure('color' , 'white');
                for d=  1:length(dayz)
                    subplot(1,length(dayz),d)
                    h1 = plotshade(xcoords{d}',PLOTs{d} , ERRORs{d},'transp' , .5 , 'patchcolor' , colz{d,2} ,...
                        'linecolor' ,colz{d,2} , 'linewidth' , 3 );
                    hold on
                    h2 = plotshade(xcoordr{d}',PLOTr{d} , ERRORr{d},'transp' , .5 , 'patchcolor' , colz{d,1} , 'linecolor' , colz{d,1} , 'linewidth' , 3 );
                    set(gca,'FontSize' , 40 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [2400 7000],'YTick' , [4000 5000 6000] ,...
                        'YTickLabels' , [4 5 6] , 'YGrid' , 'on','XLim' , [1 13]);
                    %             title(['Execution time - Training Session(s) ' , num2str(dayz{d})])
                    ylabel('Sec' )
                    xlabel('Viewing Horizon Size' )
                    hold on
                    plot(xcoords{d},PLOTs{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,2},'MarkerFaceColor',colz{d,2})
                    plot(xcoords{d},PLOTr{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,1},'MarkerFaceColor',colz{d,1})
                    area([sigSeq(d) 13 13 sigSeq(d)],[2400 2400 2600 2600] , 'FaceColor',avgCol{d},'EdgeColor' , 'none','FaceAlpha',1)
                    line([sigSeq(d) sigSeq(d)] , [2400 6800] , 'color' , avgCol{d}, 'LineWidth' , 3 , 'LineStyle' , ':')
                    %             legend([h1 h2] ,{'Structured Sequences' , 'Random Sequences'})
                end
            case 'RandStructAcrossDays'
                figure('color' , 'white');
                subplot(1,2,1);hold on
                for d=1:length(dayz)
                    h1 = plotshade(xcoords{d}',PLOTs{d} , ERRORs{d},'transp' , .5 , 'patchcolor' , colz{d,2} , 'linecolor' , colz{d,2} , 'linewidth' , 3);
                    plot(xcoords{d},PLOTs{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,2},'MarkerFaceColor',colz{d,2});
                    if ~isnan(sigMT(2,d))
                        area([1 sigMT(2,d) sigMT(2,d) 1],[6800+(d-1)*200 6800+d*200 6800+d*200 7000+(d-1)*200] , 'FaceColor',colz{d,2},'EdgeColor' , 'none','FaceAlpha',.6)
                        line([sigMT(2,d) sigMT(2,d)] , [PLOTs{d}(sigMT(2,d)) 7000+(d-1)*200] , 'color' , colz{d,2} , 'LineWidth' , 3 , 'LineStyle' , ':')
                    end
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [2600 8000],'YTick' ,...
                    [ 3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6], 'YGrid' , 'on');
                title(['Execution time for Structured Sequences'],'FontSize' , 24)
                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Viewing Horizon' ,'FontSize' , 20)
                %                 legend([h1 h2 h3] ,{'Training Session 1' , 'Training Sessions 2,3' , 'Training Sessions 4,5'})
                
                subplot(1,2,2);hold on
                for d=1:length(dayz)
                    h1 = plotshade(xcoordr{d}',PLOTr{d} , ERRORr{d},'transp' , .5 , 'patchcolor' , colz{d,1} , 'linecolor' , colz{d,1} , 'linewidth' , 3);
                    plot(xcoordr{d},PLOTr{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,1},'MarkerFaceColor',colz{d,1});
                    if ~isnan(sigMT(1,d))
                        area([1 sigMT(1,d) sigMT(1,d) 1],[6800+(d-1)*200 6800+d*200 6800+d*200 7000+(d-1)*200] , 'FaceColor', colz{d,1},'EdgeColor' , 'none','FaceAlpha',.6)
                        line([sigMT(1,d) sigMT(1,d)] , [PLOTr{d}(sigMT(1,d)) 7000+(d-1)*200] , 'color' , colz{d,1} , 'LineWidth' , 3 , 'LineStyle' , ':')
                    end
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [2600 8000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Viewing Horizon' ,'FontSize' , 20)
                %                 legend([h1 h2 h3] ,{'Session 1' , 'Sessions 2,3' , 'Sessions 4,5'})
                title(['Execution time for Random Sequences'],'FontSize' , 24)
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
            case 'BoxDaySeqType'
                %% BOX PLOT 2
                xs = zeros(2,9);
                % text x locations
                xs(1,:) = 1:1.7:1.7*9;
                for xx = 1:9
                    xs(2,xx) = xs(1,xx)+.7;
                end
                xs = xs(:);
                % text colors
                d_cols = repmat([1 2]' , 1 , 9);
                d_cols = d_cols(:);
                % star or no star
                sigstars(1,:) = {'' , '' , '' ,'','','','','','','' , '' , '' ,'','','','','',''};
                sigstars(2,:) = {'' , '' , '' ,'','*','*','*','*','*','*' , '*' , '*' ,'*','*','*','*','*','*'};
                sigstars(3,:) = {'' , '' , '*' ,'*','*','*','*','*','*','*' , '*' , '*' ,'*','*','*','*','*','*'};
                
                figure('color' , 'white');
                
                for d = 1:length(dayz)
                    cols = {colz{d,1} colz{d,2}};
                    subplot(3,1,d);hold on
                    M = getrow(MT , ismember(MT.Day , dayz{d}));
                    myboxplot(M.Horizon ,M.MT, 'notch',1 ,'plotall',0, 'fillcolor',cols,'linecolor',cols,...
                        'whiskerwidth',3,'split' , M.seqNumb,'xtickoff');
                    hold on
                    ystar = zeros(2,9);
                    c = 1;
                    for h = [1:8,13]
                        for sn = 0:1
                            M = getrow(MT , MT.seqNumb==0 & ismember(MT.Day , dayz{d}) & MT.Horizon==h);
                            ystar(sn+1,c) = (std(M.MT)/sqrt(length(M.MT)))+max(M.MT);
                            M = getrow(MT , MT.seqNumb~=0 & ismember(MT.Day , dayz{d}) & MT.Horizon==h);
                            ystar(sn+1,c) = (std(M.MT)/sqrt(length(M.MT)))+max(M.MT);
                        end
                        c = c+1;
                    end
                    ystar = ystar(:);
                    for xl = 1:length(xs)
                        text(xs(xl) , ystar(xl) , sigstars{d,xl} ,'FontSize' , 40,'color',cols{d_cols(xl)})
                    end
                    xs_labs = repmat({'Random' , 'Structured'} , 1, 9);
                    % title('Structured')
                    set(gca,'FontSize' , 30 ,  ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [0 9000],...
                        'YTick' , [1000 2000 3000 4000 5000 6000 7000 8000] , 'YTickLabels' , [1 2 3 4 5 6 7 8],...
                        'XTickLabelRotation' , 30,'YGrid' , 'on');%'XTick' , xs , 'XTickLabel' , xs_labs )
                    ylabel('Sec')
                end
            case 'LearningEffectHeat'
                h1 = figure('color' , 'white');
                hold on
                M.Day = MT.Day;
                h = 1;
                for hc  = [1:8 13]
                    [xcoords_med{h},PLOTs_med{h},ERRORs_med{h}] = lineplot([MT.seqNumb M.Day],MT.normMT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(MT.Horizon , [hc]) & ismember(MT.seqNumb , 1));
                    [xcoordr_med{h},PLOTr_med{h},ERRORr_med{h}] = lineplot([MT.seqNumb M.Day],MT.normMT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(MT.Horizon , [hc]) & ismember(MT.seqNumb , 0));
                    h  = h + 1;
                end
                close(h1)
                subplot(121)
                imagesc(cell2mat(PLOTs_med') , [2500 7000])
                set(gca , 'XTick' , [1:length(unique(M.Day))] , 'YTick', [1:9] , 'YTickLabels' , [1 :8 , 13],'FontSize' , 20)
                title('Execution Time in  Structured Sequences')
                ylabel('Viewing Horizon Size')
                xlabel('Training Session')
                
                subplot(121)
                imagesc(cell2mat(PLOTr_med'), [2500 7000])
                colorbar
                set(gca , 'XTick' , [1:length(unique(M.Day))] , 'YTick', [1:9] , 'YTickLabels' , [1 :8 , 13],'FontSize' , 20)
                title('Execution Time in  Random Sequences')
                ylabel('Viewing Horizon Size')
                xlabel('Training Session')
            case 'LearningEffectShade'
                % lump h = 6 and  up together
                MT.Horizon(MT.Horizon>5) = 5;
%                 MT  = tapply(MT , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'});
                h1 = figure('color' , 'white');
                for d=  1:length(dayz)
                    hc = 1;
                    [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(MT.Horizon,MT.normMT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 1 & ismember(MT.Day , dayz{d}));
                    hold on
                    [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(MT.Horizon,MT.normMT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 0 & ismember(MT.Day , dayz{d}));
                end
                close(h1);
                figure('color' , 'white');
                subplot(121)
                hold on
                P = cell2mat(PLOTs')';
                E = cell2mat(ERRORs')';
                X = cell2mat(xcoords')';
                % the patch color represents the horizon size
                scol = [200 200 200]/255;
                ecol = [0, 0 0]/255;
                for rgb = 1:3
                    horzcolor(rgb , :) = linspace(scol(rgb),ecol(rgb) , length(unique(ANA.Horizon)));
                end
                for i = 1:length(unique(MT.Horizon))
                    h1 = plotshade([1:length(dayz)],P(i,:) , E(i,:),'transp' , .5 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , horzcolor(:,i)  , 'linewidth' , 3);
                    plot([1:length(dayz)],P(i,:) , '-o' , 'MarkerSize' , 14 , 'color' , colz{4,2},'MarkerFaceColor',horzcolor(:,i) ,'MarkerEdgeColor',colz{4,2} , 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YLim' , [2500 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');
                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Training Session','FontSize' , 20)
                title('Effect of Practice on Performance in Structured Sequences in Different Viewing Window Sizes','FontSize' , 24)
                
                subplot(122)
                hold on
                P = cell2mat(PLOTr')';
                E = cell2mat(ERRORr')';
                X = cell2mat(xcoordr')';
                % the patch color represents the horizon size
                for i = 1:length(unique(MT.Horizon))
                    h1 = plotshade([1:length(dayz)],P(i,:) , E(i,:),'transp' , .5 , 'patchcolor' , horzcolor(:,i) , 'linecolor' ,  horzcolor(:,i) , 'linewidth' , 3);
                    plot([1:length(dayz)],P(i,:) , '-o' , 'MarkerSize' , 14 , 'color' , colz{4,1},'MarkerFaceColor',horzcolor(:,i) ,'MarkerEdgeColor',colz{4,1}, 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YLim' , [2500 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on');

                ylabel('Sec' ,'FontSize' , 20)
                xlabel('Training Session','FontSize' , 20)
                title('Effect of Practice on Performance in Random Sequences in Different Viewing Window Sizes','FontSize' , 24)
            case 'compareLearning'
                MT.Horizon(MT.Horizon>5) = 5;
                h1 = figure('color' , 'white');
                for d=  1:length(dayz)
                    hc = 1;
                    [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(MT.Horizon,MT.normMT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 1 & ismember(MT.Day , dayz{d}));
                    hold on
                    [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(MT.Horizon,MT.normMT , 'plotfcn' , 'nanmean' , 'subset' , MT.seqNumb == 0 & ismember(MT.Day , dayz{d}));
                end
                close(h1);
                horzcolor = linspace(0.2,.8 , length(unique(MT.Horizon)));
                As = cell2mat(PLOTs');
                Bs = cell2mat(ERRORs');
                Ar = cell2mat(PLOTr');
                Br = cell2mat(ERRORr');
                horz = {[1] [2] [3] [4] [5:9]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5 - 13'} , 1 , length(dayz));
                figure('color' , 'white')
                hold on
                allX = [];
                for i = [1: length(dayz)]
                    X = ((i-1)*(length(horz)))+[1:length(horz)];
                    plotshade(X , As(i,:) , Bs(i,:),'transp' , .5 , 'patchcolor' , colz{i,2}  , 'linecolor' , colz{i,2} , 'linewidth' , 3);
                    plot(X , As(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{i,2} ,'MarkerFaceColor',colz{i,2}  , 'LineWidth' , 3);
                    plotshade(X , Ar(i,:) , Br(i,:),'transp' , .5 , 'patchcolor' , colz{i,1}  , 'linecolor' , colz{i,1} , 'linewidth' , 3);
                    plot(X , Ar(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{i,1} ,'MarkerFaceColor',colz{i,1}  , 'LineWidth' , 3);
                    allX = [allX X];
                end
                set(gca,'FontSize' , 20 , 'XTick' , allX , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 max(X)], 'YLim' , [2500 7000],'YTick' ,...
                    [3000 4000 5000 6000] , 'YTickLabels' , [3 4 5 6] , 'YGrid' , 'on','XTickLabels',hlab);
                xlabel('Viewing Window Size')
                ylabel('sec')
        end
        out = [];
    case 'IPI'
        structNumb = [1 2];
        out = [];
        %         plotfcn = input('nanmean or nanmean?' , 's');
        %% IPIs vs horizon
        % this is the output of the case: 'transitions_All' that is saved to disc
        calc = 1;
        if ~ calc
            load([baseDir , '/se2_TranProb.mat'] , 'All');
            IPItable = All;
            IPItable = getrow(IPItable , ismember(IPItable.Day , dayz{day}));
        else
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
        end
        IPItable =  getrow(IPItable,ismember(IPItable.prsnumb , [4:10]));
        IPItable  = tapply(IPItable , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI' , 'nanmean(x)'});
        
        IPItable = normData(IPItable , {'IPI'});
        
        
        h0 = figure;
        [x1 , plot1 , error1] = lineplot(IPItable.Horizon , IPItable.normIPI , 'subset' , IPItable.ChunkBndry == 1 ,'plotfcn' , 'nanmean');
        hold on
        [x2 , plot2 , error2] = lineplot(IPItable.Horizon , IPItable.normIPI , 'subset' , IPItable.ChunkBndry == 2 ,'plotfcn' , 'nanmean');
        
        close(h0)
        ChunkedIPI  = pivottable(IPItable.Horizon, IPItable.ChunkBndry, IPItable.IPI ,'nanmean(x)' );
        switch nowWhat
            case 'IPIFullDispHeat'
                h0 = figure;
                hold on
                hc = 1;
                horz = {[1] [2] [3] [4] [5:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5 - 13'} , 1 , length(dayz));
                horzcolor = linspace(0.2,.8 , length(unique(MT.Horizon)));
                ipitable = IPItable;
                ipitable.Horizon(ipitable.Horizon>5) = 5;
                for h = 1:length(unique(ipitable.Horizon))
                    for d = 1:length(dayz)
                        for sn = 0:1
                            [IPI_x{sn+1,d}(hc,:),IPI_plot{sn+1,d}(hc,:) , IPI_plot_error{sn+1,d}(hc,:)] = lineplot(ipitable.prsnumb ,ipitable.normIPI , 'subset' , ipitable.Horizon == h & ismember(ipitable.seqNumb , sn) & ...
                                ismember(ipitable.Day , dayz{d}),'plotfcn' , 'nanmean');
                            hold on
                        end
                    end
                    hc = hc+1;
                end
                close(h0)
                ttl = {'Random' ['Structures ' , num2str(structNumb)]};
                figure('color' , 'white');
                for d = 1:length(dayz)
                    for sn = 0:1
                        subplot(2,length(dayz),sn*length(dayz)+d)
                        imagesc(IPI_plot{sn+1,d}, [100 600]);
                        colorbar
                        hold on
                        title([ttl{sn+1} , ' IPIs - Day(s) ' , num2str(dayz{d})])
                        xlabel('Press Number')
                        set(gca,'FontSize' , 16, 'XTick' , [1:13] ,'YTick' , [1:5] , 'YTickLabel' ,...
                            fliplr({'H = 5-13' 'H = 4' 'H = 3' 'H = 2' 'H = 1'}));
                        axis square
                    end
                end
            case 'IPIFullDispShade'
                h0 = figure;
                hold on
                hc = 1;
                
                ipitable = IPItable;
                ipitable.Horizon(ipitable.Horizon>5) = 5;
                horzcolor = repmat(linspace(0.2,.8 , length(unique(ipitable.Horizon))) , 3,1);
                for h = 1:length(unique(ipitable.Horizon))
                    for d = 1:length(dayz)
                        for sn = 0:1
                            [IPI_x{sn+1,d}(hc,:),IPI_plot{sn+1,d}(hc,:) , IPI_plot_error{sn+1,d}(hc,:)] = lineplot(ipitable.prsnumb ,ipitable.normIPI , 'subset' , ipitable.Horizon == h & ismember(ipitable.seqNumb , sn) & ...
                                ismember(ipitable.Day , dayz{d}),'plotfcn' , 'nanmean');
                            hold on
                        end
                    end
                    hc = hc+1;
                end
                close(h0)
                ttl = {'Random' ['Structures ' , num2str(structNumb)]};
                figure('color' , 'white');
                for sn = 0:1
                    subplot(2,1,sn+1)
                    hold on
                    for horzz = 1:length(unique(ipitable.Horizon))
                        for d = 1:length(dayz)
                            h1 = plotshade((horzz-1)*13+IPI_x{sn+1,d}(horzz,:),IPI_plot{sn+1,d}(horzz,:) , IPI_plot_error{sn+1,d}(horzz,:),'transp' , .5 , 'patchcolor' , horzcolor(:,horzz) , 'linecolor' , horzcolor(:,horzz) , 'linewidth' , 3);
                            plot((horzz-1)*13+IPI_x{sn+1,d}(horzz,:),IPI_plot{sn+1,d}(horzz,:) , 'o' , 'MarkerSize' , 10 , 'color' , colz{horzz,sn+1},'MarkerFaceColor',colz{horzz,sn+1});
                            title([ttl{sn+1} , ' IPIs (left to right horizons)'])
                        end
                    end
                    set(gca,'FontSize' , 20, 'XTick' , [1:(horzz-1)*13+13],'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTickLabel' , repmat([1:13] , 1,(horzz-1)));
                    xlabel('Press Number')
                end
                
                out = [];
            case 'compareLearning'
                horz = {[1] [2] [3] [4] [5] [6:13]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5' , '6 - 13'} , 1 , length(dayz));
                IPIs = IPItable;
                IPIs.Horizon(IPIs.Horizon>6) = 6;
                h1 = figure('color' , 'white');
                hold on
                for hz=  1:length(unique(IPIs.Horizon))
                    for chp=  0:2
                        hc = 1;
                        [xcoordshz{hz ,chp+1},PLOTshz{hz,chp+1},ERRORshz{hz,chp+1}] = lineplot(IPIs.Day,IPIs.normIPI , 'plotfcn' , 'nanmean' , 'subset' , IPIs.ChunkBndry == chp &...
                            ismember(IPIs.Horizon , horz{hz}));
                    end
                end
                for d=  1:length(dayz)
                    for chp=  0:2
                        [xcoords{d ,chp+1},PLOTs{d,chp+1},ERRORs{d,chp+1}] = lineplot(IPIs.Horizon,IPIs.normIPI , 'plotfcn' , 'nanmean' , 'subset' , IPIs.ChunkBndry == chp &...
                            ismember(IPIs.Day , dayz{d}));
                    end
                end
                close(h1)
                hoz = max(IPIs.Horizon);
                figure('color' , 'white');
                for d = 1:length(dayz)
                    hold on
                    xtick = [];
                    for chp = 0:2
                        h1 = plotshade(chp*(hz+1)+xcoords{d ,chp+1}',PLOTs{d,chp+1},ERRORs{d,chp+1},'transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                        plot(chp*(hz+1)+xcoords{d ,chp+1},PLOTs{d,chp+1} , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                        xtick = [xtick ; chp*(hz+1)+xcoords{d ,chp+1}];
                    end
                end
                set(gca,'FontSize' , 18, 'XTick' , xtick,'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTickLabel' , ...
                    repmat({'1' '2' '3' '4' '5-13'} , 1,chp) , 'YLim' , [250 550] , 'YTick' , [300 400 500],...
                    'YTickLabel' , [0.3 0.4 0.5]);
                xlabel('IPI Type (sub-catagory Viewing Window 1,2,3,4,5-13)' ,'FontSize' , 20)
                ylabel('Sec','FontSize' , 20)
                title('Learning Effect as a Function of Inter-Press-Interval Placement Within the Sequence' ,'FontSize' , 24)
                
                figure('color' , 'white');
                hold on
                d = length(dayz);
                for  chp = 0:2
                    xtick = [];
                    for hz=  1:length(unique(IPIs.Horizon))
                        h1 = plotshade((hz-1)*(d+1)+xcoordshz{hz ,chp+1}',PLOTshz{hz,chp+1},ERRORshz{hz,chp+1},'transp' , .5 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{4,chp+1} , 'linewidth' , 3);
                        plot((hz-1)*(d+1)+xcoordshz{hz ,chp+1},PLOTshz{hz,chp+1} , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{4,chp+1});
                        xtick = [xtick ; (hz-1)*(d+1)+xcoordshz{hz ,chp+1}];
                    end
                end
                set(gca,'FontSize' , 20, 'XTick' , xtick,'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTickLabel' ,...
                    repmat([1:d] , 1 , chp+1), 'YLim' , [210 550] , 'YTick' , [300 400 500],...
                    'YTickLabel' , [0.3 0.4 0.5]);
                title('Learning Effect as a Function of Inter-Press-Interval Placement Within the Sequence' ,'FontSize' , 24)
                ylabel('Sec','FontSize' , 20)
                xlabel('Viewing Window 1,2,3,4,5-13 (sub-catagory days)','FontSize' , 20)
                
%                 figure('color' , 'white');
%                 d = length(dayz);
%                 hold on
%                 for d=  1:length(dayz)
%                     xtick = [];
%                     for  hz=  1:length(unique(ipitable.Horizon))
%                         h1 = plotshade((hz-1)*3+xcoordsch{hz ,d}',PLOTsch{hz,d},ERRORsch{hz,d},'transp' , .5 , 'patchcolor' , horzcolor(:,d) , 'linecolor' , horzcolor(:,d) , 'linewidth' , 3);
%                         plot((hz-1)*3+xcoordsch{hz ,d}',PLOTsch{hz,d} , 'o' , 'MarkerSize' , 10 , 'color' , horzcolor(:,d),'MarkerFaceColor',horzcolor(:,d));
%                         xtick = [xtick ; (hz-1)*3+xcoordsch{hz ,d}];
%                     end
%                 end
%                 set(gca,'FontSize' , 20, 'XLim' , [-1 (hz-1)*3+2] , 'XTick' , xtick,'GridAlpha' , .2 , 'Box' , 'off','YGrid' , 'on','XTickLabel' , ...
%                     repmat({'R' , 'S_W' , 'S_B'} , 1,hz-1),'YLim' , [250 550] , 'YTick' , [300 400 500],...
%                     'YTickLabel' , [0.3 0.4 0.5]);
%                 title('Learning Effect as a Function of Inter-Press-Interval Placement Within the Sequence' ,'FontSize' , 24)
%                 ylabel('Sec','FontSize' , 20)
%                 xlabel('Viewing Window 1,2,3,4,5-13 ((sub-catagory IPI Type)','FontSize' , 20)
                
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
        end
    case 'RT'
        ANA = getrow(Dall , ismember(Dall.SN , subjnum) & Dall.isgood & ismember(Dall.seqNumb , [0 1 2]) & ~Dall.isError);
        ANA.seqNumb(ANA.seqNumb == 2) = 1;
        RT = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'RT' , 'nanmean(x)'});
        RT.RT = RT.AllPressTimes(:,1)-1500;        
        RT = normData(RT , {'RT'});
        h1 = figure('color' , 'white');
        for d=  1:length(dayz)
            hc = 1;
            [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(RT.Horizon,RT.normRT , 'plotfcn' , 'nanmean' , 'subset' , RT.seqNumb == 1 & ismember(RT.Day , dayz{d}));
            hold on
            [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(RT.Horizon,RT.normRT , 'plotfcn' , 'nanmean' , 'subset' , RT.seqNumb == 0 & ismember(RT.Day , dayz{d}));
            
        end
        close(h1);
        if poolDays
            sigSeq = [NaN 3 2];
            sigRT  = [3 4 5;4 5 4];
        else
            sigSeq = [NaN 3 3 2 2];
            sigRT  = [3 4 4 6 3;4 3 4 4 4];
        end
        
        switch nowWhat
            case 'RandvsStructCommpare'
                figure('color' , 'white');
                for d=  1:length(dayz)
                    subplot(1,length(dayz),d)
                    h1 = plotshade(xcoords{d}',PLOTs{d} , ERRORs{d},'transp' , .5 , 'patchcolor' , colz{d,2} ,...
                        'linecolor' ,colz{d,2} , 'linewidth' , 3 );
                    hold on
                    h2 = plotshade(xcoordr{d}',PLOTr{d} , ERRORr{d},'transp' , .5 , 'patchcolor' , colz{d,1} , 'linecolor' , colz{d,1} , 'linewidth' , 3 );
                    set(gca,'FontSize' , 40 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YGrid' , 'on','XLim' , [1 13],'YLim' , [500 850],'YTick' , [ 600 700 800]);
                    %             title(['Execution time - Training Session(s) ' , num2str(dayz{d})])
                    ylabel('Sec' )
                    xlabel('Viewing Horizon Size' )
                    hold on
                    plot(xcoords{d},PLOTs{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,2},'MarkerFaceColor',colz{d,2})
                    plot(xcoords{d},PLOTr{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,1},'MarkerFaceColor',colz{d,1})
%                     patch([sigSeq(d) 13 13 sigSeq(d)],[2400 2400 2600 2600] ,[NaN NaN NaN NaN],'FaceColor', avgCol{d},'EdgeColor' , 'none','FaceAlpha',1)
%                     line([sigSeq(d) sigSeq(d)] , [2400 6800] , 'color' , avgCol{d}, 'LineWidth' , 3 , 'LineStyle' , ':')
                    %             legend([h1 h2] ,{'Structured Sequences' , 'Random Sequences'})
                end
            case 'RandStructAcrossDays'
                figure('color' , 'white');
                subplot(1,2,1);hold on
                for d=1:length(dayz)
                    h1 = plotshade(xcoords{d}',PLOTs{d} , ERRORs{d},'transp' , .5 , 'patchcolor' , colz{d,2} , 'linecolor' , colz{d,2} , 'linewidth' , 3);
                    plot(xcoords{d},PLOTs{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,2},'MarkerFaceColor',colz{d,2});
%                     patch([1 sigRT(2,d) sigRT(2,d) 1],[6800+(d-1)*200 6800+d*200 6800+d*200 7000+(d-1)*200] , colz{d,2},'EdgeColor' , 'none','FaceAlpha',.6)
%                     line([sigRT(2,d) sigRT(2,d)] , [PLOTs{d}(sigRT(2,d)) 7000+(d-1)*200] , 'color' , colz{d,2} , 'LineWidth' , 3 , 'LineStyle' , ':')
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [500 850],'YTick' , [ 600 700 800] , 'YGrid' , 'on');
                title(['Execution time for Structured Sequences'])
                ylabel('Sec' )
                xlabel('Viewing Horizon' )
                %                 legend([h1 h2 h3] ,{'Training Session 1' , 'Training Sessions 2,3' , 'Training Sessions 4,5'})
                
                subplot(1,2,2);hold on
                for d=1:length(dayz)
                    h1 = plotshade(xcoordr{d}',PLOTr{d} , ERRORr{d},'transp' , .5 , 'patchcolor' , colz{d,1} , 'linecolor' , colz{d,1} , 'linewidth' , 3);
                    plot(xcoordr{d},PLOTr{d} , 'o' , 'MarkerSize' , 10 , 'color' , colz{d,1},'MarkerFaceColor',colz{d,1});
%                     patch([1 sigRT(1,d) sigRT(1,d) 1],[6800+(d-1)*200 6800+d*200 6800+d*200 7000+(d-1)*200] , colz{d,1},'EdgeColor' , 'none','FaceAlpha',.6)
%                     line([sigRT(1,d) sigRT(1,d)] , [PLOTr{d}(sigRT(1,d)) 7000+(d-1)*200] , 'color' , colz{d,1} , 'LineWidth' , 3 , 'LineStyle' , ':')
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [500 850],'YTick' , [ 600 700 800], 'YGrid' , 'on');
                ylabel('Sec' )
                xlabel('Viewing Horizon' )
                %                 legend([h1 h2 h3] ,{'Session 1' , 'Sessions 2,3' , 'Sessions 4,5'})
                title(['Execution time for Random Sequences'])
            case 'BoxAcrossDays'
                %% THE BOX PLOTS
                
                xs = zeros(3,9);
                % text x locations
                xs(1,:) = 1:2.4:2.4*9;
                for xx = 1:9
                    for d = 2:length(dayz)
                        xs(d,xx) = xs(d-1,xx)+.7;
                    end
                end
                xs = xs(:);
                % text colors
                d_cols = repmat([1:3]' , 1 , 9);
                d_cols = d_cols(:);
                % star or no star
                sigstars_r = cell(3,9);
                sigstars_r(1,:) = {'*' , '*' , '*' ,'','','','','',''};
                sigstars_r(2,:) = {'*' , '*' , '*' ,'','','','','',''};
                sigstars_r(3,:) = {'*' , '*' , '*' ,'*','','','','',''};
                sigstars_r = sigstars_r(:);
                sigstars_s = cell(3,9);
                sigstars_s(1,:) = {'*' , '*' , '*' ,'','','','','',''};
                sigstars_s(2,:) = {'*' , '*' , '*' ,'','','','','',''};
                sigstars_s(3,:) = {'*' , '*' , '*' ,'','','','','',''};
                sigstars_s = sigstars_s(:);
                ystar_r = zeros(3,9);
                ystar_s = zeros(3,9);
                c = 1;
                for h = [1:8,13]
                    for d = 1:3
                        M = getrow(RT , RT.seqNumb==0 & ismember(RT.Day , dayz{d}) & RT.Horizon==h);
                        ystar_r(d,c) = (std(M.RT)/sqrt(length(M.RT)))+max(M.RT);
                        M = getrow(RT , RT.seqNumb~=0 & ismember(RT.Day , dayz{d}) & RT.Horizon==h);
                        ystar_s(d,c) = (std(M.RT)/sqrt(length(M.RT)))+max(M.RT);
                    end
                    c = c+1;
                end
                ystar_r = ystar_r(:);
                ystar_s = ystar_s(:);
                
                
                figure('color' , 'white');
                xs_labs = repmat({'S1' , 'S2,3' , 'S4,5'} , 1, 9);
                
                subplot(211);hold on
                M = getrow(RT , RT.seqNumb~=0);
                M.Day(ismember(M.Day , [2 3])) = 2;
                M.Day(ismember(M.Day , [4 5])) = 3;
                myboxplot(M.Horizon ,M.RT, 'notch',1 ,'plotall',0, 'fillcolor',colz{1,2},'linecolor',colz{1,2},...
                    'whiskerwidth',3,'split' , M.Day,'xtickoff' );
                hold on
                for xl = 1:length(xs)
                    text(xs(xl) , ystar_s(xl) , sigstars_s{xl} ,'FontSize' , 40,'color',colz{d_cols(xl),2})
                end
                
                % title('Structured')
                set(gca,'FontSize' , 40 ,  ...
                    'GridAlpha' , .2 , 'Box' , 'off' ,...
                    'XTickLabelRotation' , 30, 'YGrid' , 'on');%,'XTick' , xs , 'XTickLabel' , xs_labs ,);
                ylabel('Sec')
                subplot(212);hold on
                M = getrow(RT , RT.seqNumb==0);
                M.Day(ismember(M.Day , [2 3])) = 2;
                M.Day(ismember(M.Day , [4 5])) = 3;
                myboxplot(M.Horizon ,M.RT, 'notch',1 ,'plotall',0, 'fillcolor',colz{1,1},'linecolor',colz{1,1},...
                    'whiskerwidth',3,'split' , M.Day,'xtickoff');
                ylabel('Sec')
                % title('Random')
                
                set(gca,'FontSize' , 40 ,  ...
                    'GridAlpha' , .2 , 'Box' , 'off' ,...
                    'XTickLabelRotation' , 30, 'YGrid' , 'on');%,'XTick' , xs , 'XTickLabel' , xs_labs ,);
                
                hold on
                for xl = 1:length(xs)
                    text(xs(xl) , ystar_r(xl) , sigstars_r{xl} ,'FontSize' , 40,'color',colz{d_cols(xl),1})
                end
            case 'BoxAcrosSeqType'
                %% BOX PLOT 2
                xs = zeros(2,9);
                % text x locations
                xs(1,:) = 1:1.7:1.7*9;
                for xx = 1:9
                    xs(2,xx) = xs(1,xx)+.7;
                end
                xs = xs(:);
                % text colors
                d_cols = repmat([1 2]' , 1 , 9);
                d_cols = d_cols(:);
                % star or no star
                sigstars(1,:) = {'' , '' , '' ,'','','','','','','' , '' , '' ,'','','','','',''};
                sigstars(2,:) = {'' , '' , '' ,'','*','*','*','*','*','*' , '*' , '*' ,'*','*','*','*','*','*'};
                sigstars(3,:) = {'' , '' , '*' ,'*','*','*','*','*','*','*' , '*' , '*' ,'*','*','*','*','*','*'};
                
                figure('color' , 'white');
                
                for d = 1:length(dayz)
                    cols = {colz{d,1} colz{d,2}};
                    subplot(3,1,d);hold on
                    M = getrow(RT , ismember(RT.Day , dayz{d}));
                    myboxplot(M.Horizon ,M.RT, 'notch',1 ,'plotall',0, 'fillcolor',cols,'linecolor',cols,...
                        'whiskerwidth',3,'split' , M.seqNumb,'xtickoff');
                    hold on
                    ystar = zeros(2,9);
                    c = 1;
                    for h = [1:8,13]
                        for sn = 0:1
                            M = getrow(RT , RT.seqNumb==0 & ismember(RT.Day , dayz{d}) & RT.Horizon==h);
                            ystar(sn+1,c) = (std(M.RT)/sqrt(length(M.RT)))+max(M.RT);
                            M = getrow(RT , RT.seqNumb~=0 & ismember(RT.Day , dayz{d}) & RT.Horizon==h);
                            ystar(sn+1,c) = (std(M.RT)/sqrt(length(M.RT)))+max(M.RT);
                        end
                        c = c+1;
                    end
                    ystar = ystar(:);
                    for xl = 1:length(xs)
                        text(xs(xl) , ystar(xl) , sigstars{d,xl} ,'FontSize' , 40,'color',cols{d_cols(xl)})
                    end
                    xs_labs = repmat({'Random' , 'Structured'} , 1, 9);
                    % title('Structured')
                    set(gca,'FontSize' , 30 ,  ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [0 9000],...
                        'YTick' , [1000 2000 3000 4000 5000 6000 7000 8000] , 'YTickLabels' , [1 2 3 4 5 6 7 8],...
                        'XTickLabelRotation' , 30,'YGrid' , 'on');%'XTick' , xs , 'XTickLabel' , xs_labs )
                    ylabel('Sec')
                end
            case 'LearningEffectHeat'
                %%
                
                h1 = figure('color' , 'white');
                hold on
                M.Day = RT.Day;
                %         M.Day(ismember(M.Day ,[2,3])) = 2;
                %         M.Day(ismember(M.Day,[4,5])) = 3;
                h = 1;
                for hc  = [1:8 13]
                    [xcoords_med{h},PLOTs_med{h},ERRORs_med{h}] = lineplot([RT.seqNumb M.Day],RT.normRT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(RT.Horizon , [hc]) & ismember(RT.seqNumb , 1));
                    [xcoordr_med{h},PLOTr_med{h},ERRORr_med{h}] = lineplot([RT.seqNumb M.Day],RT.normRT , 'plotfcn' , 'nanmean' , 'subset' ,  ismember(RT.Horizon , [hc]) & ismember(RT.seqNumb , 0));
                    h  = h + 1;
                end
                close(h1)
                
                
                subplot(121)
                imagesc(cell2mat(PLOTs_med') , [500 800])
                set(gca , 'XTick' , [1:length(unique(M.Day))] , 'YTick', [1:9] , 'YTickLabels' , [1 :8 , 13],'FontSize' , 20)
                title('Reaction Time in  Structured Sequences')
                ylabel('Viewing Horizon Size')
                xlabel('Training Session')
                axis square
                
                subplot(122)
                imagesc(cell2mat(PLOTr_med'), [500 800])
                colorbar
                set(gca , 'XTick' , [1:length(unique(M.Day))] , 'YTick', [1:9] , 'YTickLabels' , [1 :8 , 13],'FontSize' , 20)
                title('Reaction Time in  Random Sequences')
                ylabel('Viewing Horizon Size')
                xlabel('Training Session')
                axis square
            case 'LearningEffectShade'
                % lump h = 6 and  up together
                RT.Horizon(RT.Horizon>5) = 5;
                h1 = figure('color' , 'white');
                for d=  1:length(dayz)
                    hc = 1;
                    [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(RT.Horizon,RT.normRT , 'plotfcn' , 'nanmean' , 'subset' , RT.seqNumb == 1 & ismember(RT.Day , dayz{d}));
                    hold on
                    [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(RT.Horizon,RT.normRT , 'plotfcn' , 'nanmean' , 'subset' , RT.seqNumb == 0 & ismember(RT.Day , dayz{d}));
                end
                close(h1);
                figure('color' , 'white');
                subplot(121)
                hold on
                P = cell2mat(PLOTs')';
                E = cell2mat(ERRORs')';
                X = cell2mat(xcoords')';
                % the patch color represents the horizon size
                horzcolor = linspace(0.2,.8 , length(unique(RT.Horizon)));
                for i = 1:length(unique(RT.Horizon))
                    h1 = plotshade([1:length(dayz)],P(i,:) , E(i,:),'transp' , .5 , 'patchcolor' , repmat(horzcolor(i) , 1,3) , 'linecolor' , colz{3,2} , 'linewidth' , 3);
                    plot([1:length(dayz)],P(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{3,2},'MarkerFaceColor',repmat(horzcolor(i) , 1,3) , 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YLim' , [500 850],'YTick' ,...
                    [600 700 800] , 'YGrid' , 'on');
                ylabel('msec' )
                xlabel('Training Session')
                
                subplot(122)
                hold on
                P = cell2mat(PLOTr')';
                E = cell2mat(ERRORr')';
                X = cell2mat(xcoordr')';
                % the patch color represents the horizon size
                horzcolor = linspace(0.2,.8 , length(unique(RT.Horizon)));
                for i = 1:length(unique(RT.Horizon))
                    h1 = plotshade([1:length(dayz)],P(i,:) , E(i,:),'transp' , .5 , 'patchcolor' , repmat(horzcolor(i) , 1,3) , 'linecolor' , colz{3,1} , 'linewidth' , 3);
                    plot([1:length(dayz)],P(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{3,1},'MarkerFaceColor',repmat(horzcolor(i) , 1,3), 'LineWidth' , 4);
                end
                set(gca,'FontSize' , 20 , 'XTick' , [1:length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(dayz)], 'YLim' , [500 850],'YTick' ,...
                    [600 700 800] , 'YGrid' , 'on');
                ylabel('msec' )
                xlabel('Training Session')
            case 'compareLearning'
                RT.Horizon(RT.Horizon>5) = 5;
                h1 = figure('color' , 'white');
                for d=  1:length(dayz)
                    hc = 1;
                    [xcoords{d},PLOTs{d},ERRORs{d}] = lineplot(RT.Horizon,RT.normRT , 'plotfcn' , 'nanmean' , 'subset' , RT.seqNumb == 1 & ismember(RT.Day , dayz{d}));
                    hold on
                    [xcoordr{d},PLOTr{d},ERRORr{d}] = lineplot(RT.Horizon,RT.normRT , 'plotfcn' , 'nanmean' , 'subset' , RT.seqNumb == 0 & ismember(RT.Day , dayz{d}));
                end
                close(h1);
                horzcolor = linspace(0.2,.8 , length(unique(RT.Horizon)));
                As = cell2mat(PLOTs');
                Bs = cell2mat(ERRORs');
                Ar = cell2mat(PLOTr');
                Br = cell2mat(ERRORr');
                horz = {[1] [2] [3] [4] [5:9]};
                hlab = repmat({'1' , '2' , '3' , '4'  , '5 - 13'} , 1 , length(dayz));
                figure('color' , 'white')
                hold on
                allX = [];
                for i = [1: length(dayz)]
                    X = ((i-1)*(length(horz)))+[1:length(horz)];
                    plotshade(X , As(i,:) , Bs(i,:),'transp' , .5 , 'patchcolor' , colz{i,2}  , 'linecolor' , colz{i,2} , 'linewidth' , 3);
                    plot(X , As(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{i,2} ,'MarkerFaceColor',colz{i,2}  , 'LineWidth' , 3);
                    plotshade(X , Ar(i,:) , Br(i,:),'transp' , .5 , 'patchcolor' , colz{i,1}  , 'linecolor' , colz{i,1} , 'linewidth' , 3);
                    plot(X , Ar(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{i,1} ,'MarkerFaceColor',colz{i,1}  , 'LineWidth' , 3);
                    allX = [allX X];
                end
                set(gca,'FontSize' , 20 , 'XTick' , allX , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 max(X)], 'YLim' , [500 850],'YTick' ,...
                    [600 700 800] , 'YGrid' , 'on','XTickLabels',hlab);
                xlabel('Viewing Window Size')
                ylabel('msec')
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
                    id = ismember(MT.SN , subjnum) & ismember(MT.Day , dayz{d}) & ismember(MT.seqNumb , sn);
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
        MTs_pred = tapply(MTs , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT_pred' , 'nanmean(x)'});
        for d = 1:length(dayz)
            for subjnum = 1:length(subj_name)-1
                for sn = 0:1
                    [coo{sn+1}(:,d),Plot{sn+1}(:,d),err{sn+1}(:,d)] = lineplot([MTs.Horizon] , MTs.MT , 'subset' , ismember(MTs.seqNumb , sn)  & ismember(MTs.Day , dayz{d}));
                    [coo_pred{sn+1}(:,d),plot_pred{sn+1}(:,d),err_pred{sn+1}(:,d)] = lineplot([MTs_pred.Horizon] , MTs_pred.MT_pred , 'subset' , ismember(MTs_pred.seqNumb , sn) & ismember(MTs_pred.Day , dayz{d}));
                end
            end
        end
        close(h0)
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
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                ANA = MTs;
                ANA.Horizon(ANA.Horizon>6) = 6;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'},{'MT_pred' , 'nanmean(x)'});
                ANA.percChangeMT = zeros(size(ANA.MT));
                ANA.percChangeMT_pred = zeros(size(ANA.MT_pred));
                
                Daybenefit = [];
                for sn = 0:1
                    Db1= getrow(ANA , ANA.Day == 1 & ANA.seqNumb==sn);
                    Db5 = getrow(ANA , ANA.Day == length(dayz) & ANA.seqNumb==sn);
                    Db1.percChangeMT = 100*abs((Db1.MT - Db5.MT)./Db1.MT);
                    Db1.percChangeMT_pred = 100*abs((Db1.MT_pred - Db5.MT_pred)./Db1.MT_pred);
                    Daybenefit = addstruct(Daybenefit , Db1);
                end
                h1 = figure;
                hold on
                for sn = 0:1
                    [coo_red{sn+1},plot_red{sn+1},err_red{sn+1}] = lineplot([Daybenefit.Horizon] , Daybenefit.percChangeMT, 'subset', ismember(Daybenefit.seqNumb , sn));
                    [coo_pred_red{sn+1},plot_pred_red{sn+1},err_pred_red{sn+1}] = lineplot([Daybenefit.Horizon] , Daybenefit.percChangeMT_pred, 'subset', ismember(Daybenefit.seqNumb , sn));
                end
                close(h1)
                
                figure('color' , 'white')
                subplot(211)
                hold on
                for sn = 0:1
                    plotshade([1:length(unique(Daybenefit.Horizon))],plot_red{sn+1},err_red{sn+1},'transp' , .6 , 'patchcolor' , colz{3,sn+1} , 'linecolor' , colz{3,sn+1} , 'linewidth' , 3);
                    plot([1:length(unique(Daybenefit.Horizon))],plot_red{sn+1}, '-o' , 'MarkerSize' , 10 , 'color' , colz{3,sn+1},'MarkerFaceColor',colz{3,sn+1} , 'LineWidth' , 3);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Daybenefit.Horizon))], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Daybenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5' , '6-13'});
                ylabel('%' ,'FontSize' , 20)
                xlabel('Viewing Window size','FontSize' , 20)
                title('Reduction in Sequence Execution Time From First to Last Day (Actual)' ,'FontSize' , 24)
                
                subplot(212)
                hold on
                for sn = 0:1
                    plotshade([1:length(unique(Daybenefit.Horizon))],plot_pred_red{sn+1},err_pred_red{sn+1},'transp' , .6 , 'patchcolor' , colz{3,sn+1} , 'linecolor' , colz{3,sn+1} , 'linewidth' , 3);
                    plot([1:length(unique(Daybenefit.Horizon))],plot_pred_red{sn+1}, '-o' , 'MarkerSize' , 10 , 'color' , colz{3,sn+1},'MarkerFaceColor',colz{3,sn+1} , 'LineWidth' , 3);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Daybenefit.Horizon))], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Daybenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5' , '6-13'});
                ylabel('%' )
                xlabel('Viewing Window Size')
                title('Reduction in Sequence Execution Time From First to Last Day (Fitted)' ,'FontSize' , 24)
            case 'Actual&fit%ChangeDay2Day'
                Hz = {[1] [2] [3] , [4] [5:9]};
                ANA = MTs;
                ANA.Horizon(ANA.Horizon>5) = 5;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'seqNumb'} , {'MT' , 'nanmean(x)'},{'MT_pred' , 'nanmean(x)'});
                ANA.percChangeMT = zeros(length(ANA.MT) , length(dayz)-1);
                ANA.percChangeMT_pred = zeros(length(ANA.MT) , length(dayz)-1);
                
                Daybenefit = [];
                for sn = 0:1
                    for d = 1:length(dayz)-1
                        Db1= getrow(ANA , ANA.Day == d & ANA.seqNumb==sn);
                        Db = getrow(ANA , ANA.Day == d+1 & ANA.seqNumb==sn);
                        Db.percChangeMT(:,d) = 100*abs((Db1.MT - Db.MT)./Db1.MT);
                        Db.percChangeMT_pred(:,d) = 100*abs((Db1.MT_pred - Db.MT_pred)./Db1.MT_pred);
                        Daybenefit = addstruct(Daybenefit , Db);
                    end
                end
                h1 = figure;
                hold on
                for sn = 0:1
                    for d = 1:length(dayz)-1
                        [~,plot_red{sn+1}(d,:),err_red{sn+1}(d,:)] = lineplot([Daybenefit.Horizon] , Daybenefit.percChangeMT(:,d), 'subset', ismember(Daybenefit.seqNumb , sn) & ismember(Daybenefit.Day , [d+1]));
                        [~,plot_pred_red{sn+1}(d,:),err_pred_red{sn+1}(d,:)] = lineplot([Daybenefit.Horizon] , Daybenefit.percChangeMT_pred(:,d), 'subset', ismember(Daybenefit.seqNumb , sn) & ismember(Daybenefit.Day , [d+1]));
                    end
                end
                if poolDays
                    daylab = {'Session 1 to 2,3' , 'Sessions 2,3 to 4,5'};
                else
                    daylab = {'Session 1 to 2' , 'Session 2 to 3' , 'Session 3 to 4' , 'Session 4 to 5'};
                end
                close(h1)
                figure('color' , 'white')
                subplot(211)
                hold on
                Nh = length(unique(Daybenefit.Horizon));
                for sn = 0:1
                    xtick = [];
                    for h = 1:length(unique(Daybenefit.Horizon))
                        Xcoor = (h-1)*(length(dayz))+[1:length(dayz)-1];
                        plotshade(Xcoor,plot_red{sn+1}(:,h)',err_red{sn+1}(:,h)','transp' , .6 , 'patchcolor' , colz{4,sn+1} , 'linecolor' , colz{4,sn+1} , 'linewidth' , 3);
                        plot(Xcoor,plot_red{sn+1}(:,h), '-o' , 'MarkerSize' , 15 , 'color' , colz{4,sn+1},'MarkerFaceColor',colz{4,sn+1} , 'LineWidth' , 4);
                        xtick = [xtick  Xcoor];
                    end
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabels' , repmat(daylab , length(dayz)-1),...
                    'XTickLabelRotation' , 30, 'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [0 max(xtick)], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on');
                ylabel('%' ,'FontSize' , 22)
                xlabel('Viewing Window size','FontSize' , 22)
                title('Day-to-day Improvement in Performance in Different Viewing Window Sizes (W)' ,'FontSize' , 28)

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
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                ANA = IPIs;
                ANA.Horizon(ANA.Horizon>6) = 6;
                ANA  = tapply(ANA , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'IPI' , 'nanmean(x)'},{'IPI_pred' , 'nanmean(x)'});
                ANA.percChangeIPI = zeros(size(ANA.IPI));
                ANA.percChangeIPI_pred = zeros(size(ANA.IPI_pred));
                
                Daybenefit = [];
                for chp = 0:2
                    Db1= getrow(ANA , ANA.Day == 1 & ANA.ChunkBndry==chp);
                    Db5 = getrow(ANA , ANA.Day == length(dayz) & ANA.ChunkBndry==chp);
                    Db1.percChangeIPI = 100*abs((Db1.IPI - Db5.IPI)./Db1.IPI);
                    Db1.percChangeIPI_pred = 100*abs((Db1.IPI_pred - Db5.IPI_pred)./Db1.IPI_pred);
                    Daybenefit = addstruct(Daybenefit , Db1);
                end
                Daybenefit  = tapply(Daybenefit , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'percChangeIPI' , 'nanmean(x)'},{'percChangeIPI_pred' ,  'nanmean(x)'});
                Daybenefit = normData(Daybenefit , {'percChangeIPI' , 'percChangeIPI_pred'});
                h1 = figure;
                hold on
                for chp = 0:2
                    [coo_red{chp+1},plot_red{chp+1},err_red{chp+1}] = lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeIPI, 'subset', ismember(Daybenefit.ChunkBndry , chp));
                    [coo_pred_red{chp+1},plot_pred_red{chp+1},err_pred_red{chp+1}] = lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeIPI_pred, 'subset', ismember(Daybenefit.ChunkBndry , chp));
                end
                close(h1)
                
                figure('color' , 'white')
                hold on
                for chp = 0:2
                    plotshade([1:length(unique(Daybenefit.Horizon))],plot_red{chp+1},err_red{chp+1},'transp' , .5 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{3,chp+1} , 'linewidth' , 3);
                    plot([1:length(unique(Daybenefit.Horizon))],plot_red{chp+1}, '-o' , 'MarkerSize' , 10 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{3,chp+1} , 'LineWidth' , 3);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Daybenefit.Horizon))], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Daybenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5' , '6-13'});
                ylabel('%' ,'FontSize' , 20)
                xlabel('Viewing Window size','FontSize' , 20)
                title('Reduction in Inter-Press Intervals From First to Last Day (Actual)' ,'FontSize' , 24)
                
                figure('color' , 'white')
                hold on
                for chp = 0:2
                    plotshade([1:length(unique(Daybenefit.Horizon))],plot_pred_red{chp+1},err_pred_red{chp+1},'transp' , .6 , 'patchcolor' , colIPI{3,chp+1} , 'linecolor' , colIPI{3,chp+1} , 'linewidth' , 3);
                    plot([1:length(unique(Daybenefit.Horizon))],plot_pred_red{chp+1}, '-o' , 'MarkerSize' , 10 , 'color' , colIPI{3,chp+1},'MarkerFaceColor',colIPI{3,chp+1} , 'LineWidth' , 3);
                end
                set(gca,'FontSize' , 18 , 'XTick' , [1:i*length(dayz)] , ...
                    'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 length(unique(Daybenefit.Horizon))], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on',...
                    'XTick' , [1:length(unique(Daybenefit.Horizon))] , 'XTickLabels' , {'1' , '2' , '3' , '4',  '5' , '6-13'});
                ylabel('%' )
                xlabel('Viewing Window Size')
                title('Reduction in Inter-Press Intervals From First to Last Day (Fitted)' ,'FontSize' , 24)
            case 'Actual&fit%ChangeDay2Day'
                Hz = {[1] [2] [3] , [4] [5] [6:9]};
                ANA = IPIs;
                ANA.Horizon(ANA.Horizon>6) = 6;
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
                
                Daybenefit  = tapply(Daybenefit , {'Horizon' , 'Day' ,'SN' , 'ChunkBndry'} , {'percChangeIPI' , 'nanmean(x)'},{'percChangeIPI_pred' ,  'nanmean(x)'});
                Daybenefit = normData(Daybenefit , {'percChangeIPI' , 'percChangeIPI_pred'});
                h1 = figure;
                hold on
                for chp = 0:2
                    for d = 1:length(dayz)-1
                        [~,plot_red{chp+1}(d,:),err_red{chp+1}(d,:)] = lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeIPI(:,d), 'subset', ismember(Daybenefit.ChunkBndry , chp) & ismember(Daybenefit.Day , [d+1]));
                        [~,plot_pred_red{chp+1}(d,:),err_pred_red{chp+1}(d,:)] = lineplot([Daybenefit.Horizon] , Daybenefit.normpercChangeIPI_pred(:,d), 'subset', ismember(Daybenefit.ChunkBndry , chp) & ismember(Daybenefit.Day , [d+1]));
                    end
                end
                if poolDays
                    daylab = {'Session 1 to 2,3' , 'Sessions 2,3 to 4,5'};
                else
                    daylab = {'Session 1 to 2' , 'Session 2 to 3' , 'Session 3 to 4' , 'Session 4 to 5'};
                end
                close(h1)
                figure('color' , 'white')
                subplot(211)
                hold on
                Nh = length(unique(Daybenefit.Horizon));
                for chp = 0:2
                    xtick = [];
                    for h = 1:length(unique(Daybenefit.Horizon))
                        Xcoor = (h-1)*(length(dayz))+[1:length(dayz)-1];
                        plotshade(Xcoor,plot_red{chp+1}(:,h)',err_red{chp+1}(:,h)','transp' , .6 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{4,chp+1} , 'linewidth' , 3);
                        plot(Xcoor,plot_red{chp+1}(:,h), '-o' , 'MarkerSize' , 15 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{4,chp+1} , 'LineWidth' , 4);
                        xtick = [xtick  Xcoor];
                    end
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabels' , repmat(daylab , length(dayz)-1),...
                    'XTickLabelRotation' , 30, 'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [0 max(xtick)], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on');
                ylabel('%' ,'FontSize' , 22)
                xlabel('Viewing Window size','FontSize' , 22)
                title('Day-to-day Improvement in Performance in Different Viewing Window Sizes (W)' ,'FontSize' , 28)

                subplot(212)
                hold on
                Nh = length(unique(Daybenefit.Horizon));
                for chp = 0:2
                    xtick = [];
                    for h = 1:length(unique(Daybenefit.Horizon))
                        Xcoor = (h-1)*(length(dayz))+[1:length(dayz)-1];
                        plotshade(Xcoor,plot_pred_red{chp+1}(:,h)',err_pred_red{chp+1}(:,h)','transp' , .6 , 'patchcolor' , colIPI{4,chp+1} , 'linecolor' , colIPI{4,chp+1} , 'linewidth' , 3);
                        plot(Xcoor,plot_pred_red{chp+1}(:,h), '-o' , 'MarkerSize' , 15 , 'color' , colIPI{4,chp+1},'MarkerFaceColor',colIPI{4,chp+1} , 'LineWidth' , 4);
                        xtick = [xtick  Xcoor];
                    end
                end
                set(gca,'FontSize' , 20 , 'XTick' , xtick , 'XTickLabels' , repmat(daylab , length(dayz)-1),...
                    'XTickLabelRotation' , 30, 'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [0 max(xtick)], 'YLim' , [0 35],'YTick' ,...
                    [10 20 30] , 'YGrid' , 'on');
                ylabel('%' ,'FontSize' , 22)
                xlabel('Viewing Window size','FontSize' , 22)
                title('Day-to-day Improvement in Fitted Performance in Different Viewing Window Sizes (W)' ,'FontSize' , 28)
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
            eyeinfo.Hor        = [];
            eyeinfo.sn         = [];
            eyeinfo.day        = [];
            eyeinfo.BN         = [];
            eyeinfo.TN         = [];
            eyeinfo.sacPerSec  = [];
            eyeinfo.sacDur     = [];
            eyeinfo.sacPeakVel = [];
            eyeinfo.sacAmp     = [];
            eyeinfo.seqNumb    = [];
            eyeinfo.DigFixDur  = [];
            eyeinfo.prsnumb    = [];
            ANA = getrow(Dall ,ismember(Dall.seqNumb , [0:2]) & ismember(Dall.SN , subjnum) & Dall.isgood & ~Dall.isError & cellfun(@length , Dall.xEyePosDigit)>1);
           
            for tn  = 1:length(ANA.TN)
                if ismember(ANA.seqNumb , [1:2])
                    ANA.ChunkBndry(tn , :) = [1 diff(ANA.ChnkArrang(tn,:))];
                    a = find(ANA.ChunkBndry(tn , :));
                    ANA.ChunkBndry(tn , a(ignorDig:end)-1) = 3;
                    ANA.ChunkBndry(tn , ANA.ChunkBndry(tn , :) == 0) = 2;
                else
                    ANA.ChunkBndry(tn , :) = ANA.ChnkArrang(tn,:);
                end
                ANA.DigFixWeighted(tn , :) = zeros(1 ,14);
                window = 50; %samples = 100ms
                if isSymmetric
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
                else
                    for p = 1:14
                        id = ANA.xEyePosDigit{tn , 1}<=p+.3 & ANA.xEyePosDigit{tn , 1}>p-.7;
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
                goodid             = ~(abs(perv_Ben)>=3.5);
                prsnumb            = find(goodid);
                count              = sum(goodid);
                eyeinfo.PB         = [eyeinfo.PB ;perv_Ben(goodid)'];
                eyeinfo.prsnumb    = [eyeinfo.prsnumb ;find(goodid')];
                eyeinfo.CB         = [eyeinfo.CB ;ANA.ChunkBndry(tn ,goodid)'];
                eyeinfo.Hor        = [eyeinfo.Hor ; ANA.Horizon(tn)*ones(count , 1)];
                eyeinfo.sn         = [eyeinfo.sn ; ANA.SN(tn)*ones(count , 1)];
                eyeinfo.day        = [eyeinfo.day ; ANA.Day(tn)*ones(count , 1)];
                eyeinfo.BN         = [eyeinfo.BN ; ANA.BN(tn)*ones(count , 1)];
                eyeinfo.TN         = [eyeinfo.TN ; ANA.TN(tn)*ones(count , 1)];
                eyeinfo.sacPerSec  = [eyeinfo.sacPerSec ; ANA.SaccPerSec(tn)*ones(count , 1)];
                eyeinfo.sacDur     = [eyeinfo.sacDur ; mean(ANA.SaccDuration{tn})*ones(count , 1)];
                eyeinfo.sacPeakVel = [eyeinfo.sacPeakVel ; mean(ANA.SaccPeakVel{tn})*ones(count , 1)];
                eyeinfo.sacAmp     = [eyeinfo.sacAmp ; mean(ANA.SaccAmplitude{tn})*ones(count , 1)];
                eyeinfo.seqNumb    = [eyeinfo.seqNumb ; ANA.seqNumb(tn)*ones(count , 1)];
                eyeinfo.DigFixDur  = [eyeinfo.DigFixDur ;ANA.DigFixWeighted(tn ,goodid)'];
            end

            save([baseDir , '/' , filename] , 'eyeinfo','-v7.3')
        else
            load([baseDir , '/', filename])
        end
        out = [];
        switch nowWhat
            case 'sacDurSplitDay'
                Ho  = unique(eyeinfo.Hor);
                
                K = tapply(eyeinfo , {'day' , 'Hor' , 'sn','seqNumb'} , {'sacDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacDur'});
                figure('color' , 'white')
                subplot(121)
                lineplot([K.Hor] , K.normsacDur , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 0) ,'split', K.day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),1),...
                    'errorcolor' , colz(1:length(dayz),1) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),1),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),1),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),1));
                ylabel('Mean Saccade duration per trial [ms]' )
                xlabel('Viewing window Size' )
                title('Random Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [20  70],'YTick' , [30:10:70] ,...
                    'YGrid' , 'on');
                
                subplot(122)
                lineplot([K.Hor] , K.normsacDur , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 1) ,'split', K.day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),2),...
                    'errorcolor' , colz(1:length(dayz),2) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),2),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),2),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),2));
                ylabel('Mean Saccade duration per trial [ms]' )
                xlabel('Viewing window Size' )
                title('Structured Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [20  70],'YTick' , [30:10:70] ,...
                    'YGrid' , 'on');
            case 'sacDurSplitseqType'
                Ho  = unique(eyeinfo.Hor);
                K = tapply(eyeinfo , {'day' , 'Hor' , 'sn','seqNumb'} , {'sacDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacDur'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Hor] , K.normsacDur , 'plotfcn' , 'nanmean',...
                        'split', K.seqNumb , 'subset' , ismember(K.day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colz(d,1:2),...
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
                Ho  = unique(eyeinfo.Hor);
                
                K = tapply(eyeinfo , {'day' , 'Hor' , 'sn','seqNumb'} , {'sacAmp' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacAmp'});
                figure('color' , 'white')
                subplot(121)
                lineplot([K.Hor] , K.normsacAmp , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 0) ,'split', K.day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),1),...
                    'errorcolor' , colz(1:length(dayz),1) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),1),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),1),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),1));
                ylabel('Mean Saccade amplitude per trial [deg]' )
                xlabel('Viewing window Size' )
                title('Random Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [-4 20],'YTick' , [0:5:20] ,...
                    'YGrid' , 'on');
                
                subplot(122)
                lineplot([K.Hor] , K.normsacAmp , 'plotfcn' , 'nanmean',...
                    'subset' , ismember(K.seqNumb , 1) ,'split', K.day , 'leg' , 'auto' , 'linecolor' , colz(1:length(dayz),2),...
                    'errorcolor' , colz(1:length(dayz),2) , 'errorbars' , repmat({'shade'} , 1 , length(dayz)) , 'shadecolor' ,colz(1:length(dayz),2),...
                    'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , length(dayz)) , 'markerfill' , colz(1:length(dayz),2),...
                    'markersize' , 10, 'markercolor' , colz(1:length(dayz),2));
                ylabel('Mean Saccade amplitude per trial [deg]' )
                xlabel('Viewing window Size' )
                title('Structured Sequences')
                set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [-4 20],'YTick' , [0:5:20] ,...
                    'YGrid' , 'on');
            case 'sacAmpSplitseqType'
                Ho  = unique(eyeinfo.Hor);
                K = tapply(eyeinfo , {'day' , 'Hor' , 'sn','seqNumb'} , {'sacAmp' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                K.seqNumb(K.seqNumb>1) = 1;
                K = normData(K , {'sacAmp'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Hor] , K.normsacAmp , 'plotfcn' , 'nanmean',...
                        'split', K.seqNumb , 'subset' , ismember(K.day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colz(d,1:2),...
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
            case 'FixDurSplitipitype'
                Ho  = unique(eyeinfo.Hor);
                K = tapply(eyeinfo , {'day' , 'Hor' , 'sn','CB'} , {'DigFixDur' , 'nanmean'} , ...
                    'subset' , ismember(eyeinfo.prsnumb , [4:10]));
                
                K = normData(K , {'DigFixDur'});
                figure('color' , 'white')
                for d = 1:length(dayz)
                    subplot(1,length(dayz) , d)
                    lineplot([K.Hor] , K.normDigFixDur , 'plotfcn' , 'nanmean',...
                        'split', K.CB , 'subset' , ismember(K.day , dayz{d}) , 'leg' , 'auto' , 'linecolor' , colIPI(d,1:4),...
                        'errorcolor' , colIPI(d,1:4) , 'errorbars' , repmat({'shade'} , 1 , 2) , 'shadecolor' ,colIPI(d,1:4),...
                        'linewidth' , 3 , 'markertype' , repmat({'o'} , 1  , 2) , 'markerfill' , colIPI(d,1:4),...
                        'markersize' , 10, 'markercolor' , colIPI(d,1:4));
                    ylabel('Mean Saccade amplitude per trial [deg]' )
                    xlabel('Viewing window Size' )
                    title('Random Sequences')
                    set(gca,'FontSize' , 18 ,'GridAlpha' , .2 , 'Box' , 'off' , 'YLim' , [20 140],'YTick' , [30:10:140] ,...
                        'YGrid' , 'on');
                    if d>1
                        set(gca,'YColor' , 'none');
                    end
                end    
               
 
        end
        
    
end


