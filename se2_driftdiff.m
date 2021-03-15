function se2_driftdiff(what)

N = 14;   % number of presses in a sequence
% colors = rand(N , 3);
%load('/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze/se2_colors.mat')
load('/Volumes/MotorControl/data/SeqEye2/analyze/se2_colors.mat');

switch what
    case 'Scenario1'
        % assumption 1: you can only execute a digit that has been a 100% planned
        % assumption 2: execution pushes your planning percentages one down the exponential ladder
        % assumtion  3: you only worry about the xx digits ahead that you can plen/see. Even if you see more, or can plan more does not matter 
        Tp = 100; % time you need to plan the 100% of 1 press
        Te = 150; % time of executing onepress when you have one press ahead planned
        beta =  0.1 ;
        gamma = 1.1; % the factor by which the planning time increases if sonr in synchrony with execution
        Tp_plan_while_exe = gamma*Tp;
        for Horizon = 1:14
            for B = 1:13
                
                xx = min(B , Horizon);
                syms Digit_start Digit_mid PlanPerc plantime planNum percPlanning PlanPerUnit plantimePort ExePerUnit DIFPLAN
                Plan_Decay_start = 100*exp((Digit_start/B)*(1-Digit_start));   % happen at the beggining, so no pressing
                Plan_Decay_mid   = 70*exp((Digit_mid/B)*(1-Digit_mid));     % happen in the middle while pressing
                Te_X1_planned_ahead = ExePerUnit + ExePerUnit*exp(2*(2-xx));
                Te_exe_while_plan   = Te + Te*(DIFPLAN/100);
                % Careful that plantime can be relating to different percentages of a planning a press
                plantimePortion         =  (2.7-exp(.01*(100-PlanPerc)))/1.7 ;   % the symbolic PlanPerc is the portion of 100 that you are planning a digit
                % plantime*PlanPerUnit is the time you need to plan any given percentage of a number
                % PlanPerUnit is either Tp (if not executing) or Tp_plan_while_exe (if executing)
                % plantimePort is the output of plantimePortion symbolic function
                
                Tp_moreThanOne   = (plantimePort * PlanPerUnit)*(1+beta*(planNum-1));   
                
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                
                alpha = .6; % the percentage of presses in the horizon that you fully plan before execution
                full_plan_lim = ceil(alpha*xx);
                prep_perc = zeros(1,N);
                prep_perc(1:full_plan_lim) =  100;
                rem_plan = xx - full_plan_lim; 
                Digit_start = 1:rem_plan;
                prep_perc(full_plan_lim+1:full_plan_lim+rem_plan) = eval(Plan_Decay_start);
                
                MT(B , Horizon).com_plan = zeros(N , xx);
                MT(B , Horizon).prs_numb = zeros(N , xx);
                MT(B , Horizon).dif_plan = zeros(N , xx);
                MT(B , Horizon).PlanTimePerPress = zeros(N , xx);
                
                for press = 1:N
                    press
                    id = min(xx , N-press+1);
                    switch press
                        case 1  % this is essentially reaction time
                            ToBePlanned = prep_perc(press:press+(id-1));
                            MT(B , Horizon).prs_numb(press ,1 :id) = press : press+id-1; % the presses you're planning
                            planNum = id;
                            PlanPerUnit = Tp;
                            tempid =  [];
%                             midl_plan_perc = 100;
                            for horz = 1:xx
                                Digit_start = horz+1;
                                midl_plan_perc(horz) = eval(Plan_Decay_start); % the planning top-up that you can add to the previous planned ones
                            end
                        otherwise
                            end_id = 2:size(MT(B , Horizon).com_plan , 2);
                            MT(B , Horizon).prs_numb(press ,1 :id) = press : press+id-1;
                            if ~isempty(end_id)
                                ToBePlanned = [MT(B , Horizon).com_plan(press-1 , end_id) prep_perc(press + id-1)];
                            else
                                ToBePlanned = [prep_perc(press)];
                            end
                            
                            planNum = id;
                            PlanPerUnit = Tp_plan_while_exe;
                            tempid =  ToBePlanned == 100;
%                             midl_plan_perc = 50;
                            for horz = 1:xx
                                Digit_mid = horz+1;
                                midl_plan_perc(horz) = eval(Plan_Decay_mid); % the planning top-up that you can add to the previous planned ones
                            end
                            
                    end
                    MT(B , Horizon).plan_Done(press,1:id)  = midl_plan_perc(1:id);
                    [B , Horizon]
                    MT(B , Horizon).plan_Done(press,:) = [MT(B , Horizon).plan_Done(press,1:id) zeros(1,xx-id)] + ToBePlanned; % MT(B , Horizon).plan_Done(press,:) is the planned movements after the planning limits while pressing have been added
                    % if after adding the execc plannig, a certian press exceeds 100, divide the axxecc planning between the other
                    temp3 = (MT(B , Horizon).plan_Done(press,:) > 100); % 
                    if sum(temp3)
                       temp4 = sum(MT(B , Horizon).plan_Done(press , temp3))  - 100*sum(temp3);
                       MT(B , Horizon).plan_Done(press,temp3) = 100;
                       temp5 = (MT(B , Horizon).plan_Done(press,:) < 100 & MT(B , Horizon).plan_Done(press,:)~=0);
                    end
                    temp3 = (MT(B , Horizon).plan_Done(press,:) > 100); % 
                    if sum(temp3)
                       temp4 = sum(MT(B , Horizon).plan_Done(press,temp3))  - 100*sum(temp3);
                       MT(B , Horizon).plan_Done(press,temp3) = 100;
                       temp5 = (MT(B , Horizon).plan_Done(press,:) < 100 & MT(B , Horizon).plan_Done(press,:)~=0);
                    end
                    % if the immediate press that you have to d o right after the press, is not fully prepared, you have prepare it fully then
                    
                    MT(B , Horizon).com_plan(press , :)        = MT(B , Horizon).plan_Done(press,:);
                    if press>1
                        MT(B , Horizon).dif_plan(press , :)           = MT(B , Horizon).plan_Done(press,:) - ToBePlanned;
                    else
                        MT(B , Horizon).dif_plan(press , :) = MT(B , Horizon).plan_Done(press,:); %right before the first press you are planning everything from scratch
                    end
                    PlanPerc                      = MT(B , Horizon).dif_plan(press , :);
%                   plantimePort                  = eval(plantimePortion);
                    plantimePort                = PlanPerc/100;
                    plantimePort(MT(B , Horizon).dif_plan(press , :) == 0) = 0;
                    MT(B , Horizon).PlanTimePerPress(press , :) =  eval(Tp_moreThanOne); % this is as if your planning each of them this way
                    
                    if MT(B , Horizon).plan_Done(press,1) < 100
                        PlanPerUnit = Tp; % because the execution has to wait for the  prep, so you
                        PlanPerc = 100*ones(size(MT(B , Horizon).plan_Done(press,:))) - MT(B , Horizon).plan_Done(press,:);
                        plantimePort     = eval(plantimePortion);
                        MT(B , Horizon).PlanTimePerPress(press , :) = MT(B , Horizon).PlanTimePerPress(press , :) + eval(Tp_moreThanOne); % this is as if your planning each of them this way
                    end
                    
                end
                for press = 1:N-1
                    if sum(MT(B , Horizon).PlanTimePerPress(press+1,:) , 2)
                        DIFPLAN  = sum(MT(B , Horizon).dif_plan(press+1 , :));
                        ExePerUnit = eval(Te_exe_while_plan);
                    else
                        ExePerUnit = Te;
                        DIFPLAN  = sum(MT(B , Horizon).dif_plan(press+1 , :));
                    end
                    MT(B , Horizon).Exec_Time(press) = eval(Te_X1_planned_ahead);
                end
                press = N;
                ExePerUnit = Te;
                DIFPLAN  = 0;
                MT(B , Horizon).Exec_Time(press) = eval(Te_X1_planned_ahead);
                
                MT(B , Horizon).Reac_Time = sum(MT(B , Horizon).PlanTimePerPress(1,:));
                MT(B , Horizon).Plan_Time = sum(MT(B , Horizon).PlanTimePerPress(2:end,:) , 2);
                MT(B , Horizon).Exec_Time(MT(B , Horizon).Plan_Time == 0) = Te;
                a = [0 ;cumsum(MT(B , Horizon).Plan_Time)]';
                b = cumsum(MT(B , Horizon).Exec_Time);
                for press = 1:N
                    temp(press) = b(press);%max(a(press) , b(press));
                end
                MT(B , Horizon).IPI = diff(temp);
                MT(B , Horizon).Move_Time = sum(MT(B , Horizon).IPI);
            end
        end
        
    case 'Scenario2' % negative planning
        % assumption 1: you can only execute a digit that has been a 100% planned
        % assumption 2: execution pushes your planning percentages one down the exponential ladder
        % assumtion  3: you only worry about the xx digits ahead that you can plen/see. Even if you see more, or can plan more does not matter 
        Tp = 100; % time you need to plan the 100% of 1 press
        Te = 150; % time of executing onepress when you have one press ahead planned
        beta =  0.1 ;
        gamma = 1.1; % the factor by which the planning time increases if sonr in synchrony with execution
        Tp_plan_while_exe = gamma*Tp;
        for Horizon = 1:14
            for B = 1:13
                prePlan = zeros(1,5);
                xx = min(B , Horizon);
                syms Digit_start Digit_mid PlanPerc plantime planNum percPlanning PlanPerUnit plantimePort ExePerUnit DIFPLAN
                Plan_Decay_start = 100*exp((Digit_start/B)*(1-Digit_start));   % happen at the beggining, so no pressing
                Plan_Decay_mid   = 70*exp((Digit_mid/B)*(1-Digit_mid));     % happen in the middle while pressing
                Te_X1_planned_ahead = ExePerUnit + ExePerUnit*exp(2*(2-xx));
                Te_exe_while_plan   = Te + Te*(DIFPLAN/100);
                % Careful that plantime can be relating to different percentages of a planning a press
                plantimePortion         =  (2.7-exp(.01*(100-PlanPerc)))/1.7 ;   % the symbolic PlanPerc is the portion of 100 that you are planning a digit
                % plantime*PlanPerUnit is the time you need to plan any given percentage of a number
                % PlanPerUnit is either Tp (if not executing) or Tp_plan_while_exe (if executing)
                % plantimePort is the output of plantimePortion symbolic function
                
                Tp_moreThanOne   = (plantimePort * PlanPerUnit)*(1+beta*(planNum-1));   
                
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                
                alpha = .6; % the percentage of presses in the horizon that you fully plan before execution
                full_plan_lim = ceil(alpha*xx);
                prep_perc = zeros(1,N);
                prep_perc(1:full_plan_lim) =  100;
                rem_plan = xx - full_plan_lim; 
                Digit_start = 1:rem_plan;
                prep_perc(full_plan_lim+1:full_plan_lim+rem_plan) = eval(Plan_Decay_start);
                
                MT(B , Horizon).com_plan = zeros(N , xx);
                MT(B , Horizon).prs_numb = zeros(N , xx);
                MT(B , Horizon).dif_plan = zeros(N , xx);
                MT(B , Horizon).PlanTimePerPress = zeros(N , xx);
                
                for press = 1:N
                    press
                    id = min(xx , N-press+1);
                    switch press
                        case 1  % this is essentially reaction time
                            ToBePlanned = prep_perc(press:press+(id-1));
                            MT(B , Horizon).prs_numb(press ,1 :id) = press : press+id-1; % the presses you're planning
                            planNum = id;
                            PlanPerUnit = Tp;
                            tempid =  [];
%                             midl_plan_perc = 100;
                            for horz = 1:xx
                                Digit_start = horz+1;
                                midl_plan_perc(horz) = eval(Plan_Decay_start); % the planning top-up that you can add to the previous planned ones
                            end
                        otherwise
                            end_id = 2:size(MT(B , Horizon).com_plan , 2);
                            MT(B , Horizon).prs_numb(press ,1 :id) = press : press+id-1;
                            if ~isempty(end_id)
                                ToBePlanned = [MT(B , Horizon).com_plan(press-1 , end_id) prep_perc(press + id-1)];
                            else
                                ToBePlanned = [prep_perc(press)];
                            end
                            
                            planNum = id;
                            PlanPerUnit = Tp_plan_while_exe;
                            tempid =  ToBePlanned == 100;
%                             midl_plan_perc = 50;
                            for horz = 1:xx
                                Digit_mid = horz+1;
                                midl_plan_perc(horz) = eval(Plan_Decay_mid); % the planning top-up that you can add to the previous planned ones
                            end
                            
                    end
                    MT(B , Horizon).plan_Done(press,1:id)  = midl_plan_perc(1:id);
                    [B , Horizon]
                    MT(B , Horizon).plan_Done(press,:) = [MT(B , Horizon).plan_Done(press,1:id) zeros(1,xx-id)] + ToBePlanned; % MT(B , Horizon).plan_Done(press,:) is the planned movements after the planning limits while pressing have been added
                    % if after adding the execc plannig, a certian press exceeds 100, divide the axxecc planning between the other
                    temp3 = (MT(B , Horizon).plan_Done(press,:) > 100); % 
                    if sum(temp3)
                       temp4 = sum(MT(B , Horizon).plan_Done(press , temp3))  - 100*sum(temp3);
                       MT(B , Horizon).plan_Done(press,temp3) = 100;
                       temp5 = (MT(B , Horizon).plan_Done(press,:) < 100 & MT(B , Horizon).plan_Done(press,:)~=0);
                    end
                    temp3 = (MT(B , Horizon).plan_Done(press,:) > 100); % 
                    if sum(temp3)
                       temp4 = sum(MT(B , Horizon).plan_Done(press,temp3))  - 100*sum(temp3);
                       MT(B , Horizon).plan_Done(press,temp3) = 100;
                       temp5 = (MT(B , Horizon).plan_Done(press,:) < 100 & MT(B , Horizon).plan_Done(press,:)~=0);
                    end
                    % if the immediate press that you have to d o right after the press, is not fully prepared, you have prepare it fully then
                    
                    MT(B , Horizon).com_plan(press , :)        = MT(B , Horizon).plan_Done(press,:);
                    if press>1
                        MT(B , Horizon).dif_plan(press , :)           = MT(B , Horizon).plan_Done(press,:) - ToBePlanned;
                    else
                        MT(B , Horizon).dif_plan(press , :) = MT(B , Horizon).plan_Done(press,:); %right before the first press you are planning everything from scratch
                    end
                    PlanPerc                      = MT(B , Horizon).dif_plan(press , :);
%                   plantimePort                  = eval(plantimePortion);
                    plantimePort                = PlanPerc/100;
                    plantimePort(MT(B , Horizon).dif_plan(press , :) == 0) = 0;
                    MT(B , Horizon).PlanTimePerPress(press , :) =  eval(Tp_moreThanOne); % this is as if your planning each of them this way
                    
                    if MT(B , Horizon).plan_Done(press,1) < 100
                        PlanPerUnit = Tp; % because the execution has to wait for the  prep, so you
                        PlanPerc = 100*ones(size(MT(B , Horizon).plan_Done(press,:))) - MT(B , Horizon).plan_Done(press,:);
                        plantimePort     = eval(plantimePortion);
                        MT(B , Horizon).PlanTimePerPress(press , :) = MT(B , Horizon).PlanTimePerPress(press , :) + eval(Tp_moreThanOne); % this is as if your planning each of them this way
                    end
                    
                end
                for press = 1:N-1
                    if sum(MT(B , Horizon).PlanTimePerPress(press+1,:) , 2)
                        DIFPLAN  = sum(MT(B , Horizon).dif_plan(press+1 , :));
                        ExePerUnit = eval(Te_exe_while_plan);
                    else
                        ExePerUnit = Te;
                        DIFPLAN  = sum(MT(B , Horizon).dif_plan(press+1 , :));
                    end
                    MT(B , Horizon).Exec_Time(press) = eval(Te_X1_planned_ahead);
                end
                press = N;
                ExePerUnit = Te;
                DIFPLAN  = 0;
                MT(B , Horizon).Exec_Time(press) = eval(Te_X1_planned_ahead);
                
                MT(B , Horizon).Reac_Time = sum(MT(B , Horizon).PlanTimePerPress(1,:));
                MT(B , Horizon).Plan_Time = sum(MT(B , Horizon).PlanTimePerPress(2:end,:) , 2);
                MT(B , Horizon).Exec_Time(MT(B , Horizon).Plan_Time == 0) = Te;
                a = [0 ;cumsum(MT(B , Horizon).Plan_Time)]';
                b = cumsum(MT(B , Horizon).Exec_Time);
                for press = 1:N
                    temp(press) = b(press);%max(a(press) , b(press));
                end
                MT(B , Horizon).IPI = diff(temp);
                MT(B , Horizon).Move_Time = sum(MT(B , Horizon).IPI);
            end
        end
end


% ********************  Reaction time
figure('color' , 'white');
for h = 1:14
    for x = 1:13
        RT(x,h) = MT(x , h).Reac_Time; 
    end
end
for x = 1:13
    plot(RT(x ,:) ,'color',colors(x,:), 'LineWidth' , 3)
    hold on
end
title(['Initial reaction time for different buffer sizes vs. Horizons'])
ylabel('msec'  ,'FontSize' , 20)
xlabel('Horizon')
ax.XTick = [1:14];
legend({'B = 1' , 'B = 2' , 'B = 3', 'B = 4', 'B = 5', 'B = 6', 'B = 7', 'B = 8', 'B = 9', 'B = 10', 'B = 11', 'B = 12', 'B = 13', 'B = 14'})
set(gca ,'XTick'  ,[1:14], 'YTick',[100:100:2500],'Box' , 'off','FontSize' , 20 , 'YLim' , [0 2500])
grid on
% ******************** Planning time
figure('color' , 'white');
for h = 14
    for x = 1:13
        PT(x, :) = MT(x , h).Plan_Time; 
    end
end
for x = 1:13
    plot(PT(x ,:) ,'color',colors(x,:), 'LineWidth' , 3)
    hold on
end
title(['Planning time per press for different buffer sizes in full horizon'])
ylabel('msec'  ,'FontSize' , 20)
xlabel('Presses')
ax.XTick = [1:14];
legend({'B = 1' , 'B = 2' , 'B = 3', 'B = 4', 'B = 5', 'B = 6', 'B = 7', 'B = 8', 'B = 9', 'B = 10', 'B = 11', 'B = 12', 'B = 13', 'B = 14'})
set(gca ,'XTick'  , [1:14],'Box' , 'off','FontSize' , 20)
grid on
% ******************* Movement time
figure('color' , 'white');
for h = 1:14
    for x = 1:13
        mt(x,h) = MT(x , h).Move_Time; 
    end
end
for x = 1:13
    plot(mt(x ,:) ,'color',colors(x,:), 'LineWidth' , 3)
    hold on
end
title(['Movement time for different buffer sizes (B) vs. Horizons'])
ylabel('msec' )
xlabel('Horizon')
ax.XTick = [1:14];
set(gca ,'XTick'  , [1:14],'Box' , 'off','FontSize' , 20)
legend({'B = 1' , 'B = 2' , 'B = 3', 'B = 4', 'B = 5', 'B = 6', 'B = 7', 'B = 8', 'B = 9', 'B = 10', 'B = 11', 'B = 12', 'B = 13', 'B = 14'})
grid on

% ******************* IPI
for x = 3:8
    figure('color' , 'white');
    % buffer of interest
    
    for h = 1:14
        IPI(h , :) = MT(x , h).IPI    ;
    end
    subplot(1,2,1)
    for h = 1:14
        plot(IPI(h ,:) ,'color',colors(h,:), 'LineWidth' , 3)
        hold on
    end
    title(['IPIs for buffer size ' , num2str(x) , ' in different Horizons'])
    ylabel('msec'  ,'FontSize' , 20)
    xlabel('Presses')
    ax.XTick = [1:14];
    legend({'H = 1' , 'H = 2' , 'H = 3', 'H = 4', 'H = 5', 'H = 6', 'H = 7', 'H = 8', 'H = 9', 'H = 10', 'H = 11', 'H = 12', 'H = 13', 'H = 14'})
    grid on
    
    subplot(1,2,2)
    imagesc(IPI)
    title(['IPIs for buffer size ' , num2str(x) , ' in different Horizons'])
    ylabel('Horizon size')
    xlabel('IPI')
    set(gca , 'XTick'  , [1:13] , 'YTick' , [1:14] ,'Box' , 'off','FontSize' , 20)
end

%**************** Planning decay
figure('color' , 'white');
for B = 1:14
    x = 1:14;
    plot(x , 100*exp((x/B).*(1-x)), 'LineWidth' , 3);
    hold on
end
title('Planning decay function for different buffer sizes')
ylabel('Percent')
xlabel('Presses')
legend({'B = 1' , 'B = 2' , 'B = 3', 'B = 4', 'B = 5', 'B = 6', 'B = 7', 'B = 8', 'B = 9', 'B = 10', 'B = 11', 'B = 12', 'B = 13', 'B = 14'})
set(gca , 'XTick'  , [1:14] ,'Box' , 'off','FontSize' , 20)
grid on

%**************** Plannings slows down execuition 
figure('color' , 'white');
DP = 0:100;
plot(DP , 1+DP/100, 'LineWidth' , 3);
title('Increase of execution time as a function of accumulative planning percentage ')
xlabel('Percent')
ylabel('Increase fold')
set(gca ,'Box' , 'off','FontSize' , 20)
grid on

%**************** Planning more takes more time (function of accumulative % and number of digits)
PP = 1:100;
pport = (2.7-exp(.01*(100-PP)))/1.7;
plot(PP , 1+pport , 'LineWidth' , 3);



figure('color' , 'white');
for pN = 1:14
    T_inc_fold(:,pN) = (pport)*(beta*(pN-1));
    plot(PP , T_inc_fold(:,pN), 'LineWidth' , 3  ,'color',colors(pN,:));
    hold on
end
title('Increase of planning time as a function of accumulative planning percentage and buffer size')
xlabel('Percent')
ylabel('Increase fold')
set(gca ,'Box' , 'off','FontSize' , 20)
legend({'B = 1' , 'B = 2' , 'B = 3', 'B = 4', 'B = 5', 'B = 6', 'B = 7', 'B = 8', 'B = 9', 'B = 10', 'B = 11', 'B = 12', 'B = 13', 'B = 14'})
grid on

%**************** Having more digits planned ahead speeds up execution 
DP = 100;
figure('color' , 'white');
x = 1:14;
plot(x , 1 + .01*((100./DP)*exp(2*(2-x))), 'LineWidth' , 3);
title('Decrease of Execution time as a function of buffer size')
xlabel('Buffer size')
ylabel('Decrease fold')
grid on
set(gca ,'XTick'  , [1:14] , 'Box' , 'off','FontSize' , 20)

figure('color' , 'white');
x = 7;
imagesc(a(1:14,1:14))
colorbar
xlabel('Digits')
ylabel('Before press')
title(['Planning evolution in buffer size ' , num2str(x)])
set(gca ,'XTick'  , [1:14] , 'YTick'  , [1:14],'Box' , 'off','FontSize' , 20)


end
