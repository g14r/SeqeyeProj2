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
                
                
                
                for d = 1:length(dayz)
                    subplot(2,length(dayz),length(dayz)+d)
                    hold on
                    for chp = 0:2
                        h1 = plotshade(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)' ,'transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                        plot(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                    end
                    set(gca,'FontSize' , 20 , 'XTick' , [1:8,13] , 'XTickLabel' , {'1' '2' '3' '4' '5' '6' '7' '8' '13'} , ...
                        'GridAlpha' , .2 , 'Box' , 'off' , 'XLim' , [1 13], 'YLim' , [210 600],'YTick' ,...
                        [300 400 500 ] , 'YTickLabels' , [0.3 0.4 0.5] , 'YGrid' , 'on');
                    ylabel('Sec' )
                    xlabel('Viewing Horizon' )
                    title(['fitted IPIs - Day(s) ' , num2str(dayz{d})])
                end
            case 'Actual&fitDayz'
                Hz = {[1] [2] [3] , [4] [5:9]};
                ANA = MT;
                
                h1 = figure;
                for d = 1:length(dayz)
                    ANA = getrow(MT , ismember(MT.seqNumb , [1])  & ismember(MT.Day , d));
                    ANA.Horizon(ANA.Horizon>5) = 5;
                    [coo_red{chp+1}(:,d),Plot_red{chp+1}(:,d),err_red{chp+1}(:,d)] = lineplot([ANA.Horizon] , ANA.MT );
                    [coo_pred_red{chp+1}(:,d),plot_pred_red{chp+1}(:,d),err_pred_red{chp+1}(:,d)] = lineplot([ANA.Horizon] , ANA.MT_pred );
                    
                    
                    ANA = getrow(MT , ismember(MT.seqNumb , [0])  & ismember(MT.Day , d));
                    ANA.Horizon(ANA.Horizon>5) = 5;
                    [coo_red{chp+1}(:,d),Plot_red{chp+1}(:,d),err_red{chp+1}(:,d)] = lineplot([ANA.Horizon] , ANA.MT );
                    [coo_pred_red{chp+1}(:,d),plot_pred_red{chp+1}(:,d),err_pred_red{chp+1}(:,d)] = lineplot([ANA.Horizon] , ANA.MT_pred);
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
                    plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],Plot_red{chp+1}(i,:),err_red{chp+1}(i,:),'transp' , .6 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , colz{6-i+1,2} , 'linewidth' , 3);
                    plot((i-1)*(length(dayz)+1)+[1:length(dayz)],Plot_red{chp+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{6-i+1,2},'MarkerFaceColor',colz{6-i+1,2} , 'LineWidth' , 3);
                    
                    plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],Plot_red{chp+1}(i,:),err_red{chp+1}(i,:),'transp' , .6 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , colz{6-i+1,1} , 'linewidth' , 3);
                    plot((i-1)*(length(dayz)+1)+[1:length(dayz)],Plot_red{chp+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{6-i+1,1},'MarkerFaceColor',colz{6-i+1,1} , 'LineWidth' , 3);
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
                    plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{chp+1}(i,:),err_pred_red{chp+1}(i,:),'transp' , .6 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , colz{6-i+1,2} , 'linewidth' , 3);
                    plot((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{chp+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{6-i+1,2},'MarkerFaceColor',colz{6-i+1,2} , 'LineWidth' , 3);
                    
                    plotshade((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{chp+1}(i,:),err_pred_red{chp+1}(i,:),'transp' , .6 , 'patchcolor' , horzcolor(:,i) , 'linecolor' , colz{6-i+1,1} , 'linewidth' , 3);
                    plot((i-1)*(length(dayz)+1)+[1:length(dayz)],plot_pred_red{chp+1}(i,:) , '-o' , 'MarkerSize' , 10 , 'color' , colz{6-i+1,1},'MarkerFaceColor',colz{6-i+1,1} , 'LineWidth' , 3);
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
                Hz = {[1] [2] [3] , [4] [5:9]};
                ANA = MT;
                ANA.Horizon(ANA.Horizon>5) = 5;
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
                ANA = MT;
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
                Hz = {[1] [2] [3] , [4] [5:9]};
                ANA = MT;
                ANA.Horizon(ANA.Horizon>=5) = 5;
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
            case 'Actual&fit%ChangeSeqType'
                Hz = {[1] [2] [3] , [4] [5:9]};
                ANA = MT;
                ANA.Horizon(ANA.Horizon>=5) = 5;
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
                    h1 = plotshade(coo{chp+1}(:,d)',Plot{chp+1}(:,d)',err{chp+1}(:,d)','transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                    plot(coo{chp+1}(:,d)',Plot{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                    h1 = plotshade(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)','transp' , .3 , 'patchcolor' , [.5 .5 .5] , 'linecolor' , [.5 .5 .5] , 'linewidth' , 3);
                    plot(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)' , 'o' , 'MarkerSize' , 8 , 'color' , [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
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
                    h1 = plotshade(coo{chp+1}(:,d)',Plot{chp+1}(:,d)',err{chp+1}(:,d)','transp' , .5 , 'patchcolor' , colIPI{d,chp+1} , 'linecolor' , colIPI{d,chp+1} , 'linewidth' , 3);
                    plot(coo{chp+1}(:,d)',Plot{chp+1}(:,d)' , 'o' , 'MarkerSize' , 10 , 'color' , colIPI{d,chp+1},'MarkerFaceColor',colIPI{d,chp+1});
                    h1 = plotshade(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)','transp' , .3 , 'patchcolor' , [.5 .5 .5] , 'linecolor' , [.5 .5 .5] , 'linewidth' , 3);
                    plot(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)' , 'o' , 'MarkerSize' , 8 , 'color' , [.5 .5 .5],'MarkerFaceColor',[.5 .5 .5]);
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
                for d = 1:length(dayz)
                    errorbar(coo{chp+1}(:,d)',Plot{chp+1}(:,d)',err{chp+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13], 'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
                hold on
                xlabel('Horizon')
                ylabel('msec')
                title('Chunked MT vs. Horizon')
                grid on
                
                subplot(222)
                for d = 1:length(dayz)
                    errorbar(coo{chp+1}(:,d)',Plot{chp+1}(:,d)',err{chp+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
                xlabel('Horizon')
                ylabel('msec')
                title('Random MT vs. Horizon')
                legend({'Day1' , 'Day2' ,'Day3','Day4','Day5'}, 'Box' , 'off')
                grid on
                
                subplot(223)
                for d = 1:length(dayz)
                    errorbar(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20, 'GridAlpha' , 1, 'Box' , 'off')
                hold on
                xlabel('Horizon')
                ylabel('msec')
                title('Fitted Chunked MT vs. Horizon')
                grid on
                
                subplot(224)
                for d = 1:length(dayz)
                    errorbar(coo_pred{chp+1}(:,d)',plot_pred{chp+1}(:,d)',err_pred{chp+1}(:,d)' , 'LineWidth' , 3);
                    hold on
                end
                set(gca, 'YLim' , [3000 , 8000] , 'XTick' , [1:8 , 13],'FontSize' , 20 , 'GridAlpha' , 1, 'Box' , 'off')
                xlabel('Horizon')
                title('Fitted Random MT vs. Horizon')
                grid on
            case 'plotCoef'
                h0 = figure;
                if poolDays
                    MT.Day(MT.Day==3) = 2;
                    MT.Day(ismember(MT.Day , [4 5])) = 3;
                end
                for sn = 0:1
                    [coo_b1{sn+1},plot_b1{sn+1},err_b1{sn+1}] = lineplot([MT.Day] , MT.b1 , 'subset' , ismember(MT.seqNumb , [sn]));
                    [coo_b2{sn+1},plot_b2{sn+1},err_b2{sn+1}] = lineplot([MT.Day] , MT.b2 , 'subset' , ismember(MT.seqNumb , [sn]));
                    [coo_b3{sn+1},plot_b3{sn+1},err_b3{sn+1}] = lineplot([MT.Day] , MT.b3 , 'subset' , ismember(MT.seqNumb , [sn]));
                    [coo_invb3{sn+1},plot_invb3{sn+1},err_invb3{sn+1}] = lineplot([MT.Day] , (MT.b3).^-1 , 'subset' , ismember(MT.seqNumb , [sn]));
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