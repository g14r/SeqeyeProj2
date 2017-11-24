function se2_Hypo(what)


Tp = 150; % Planning time per X ot less digits
Te = 200; % Execution time per each press
N = 14;   % number of presses in a sequence
colors = rand(N , 3);

switch what
    case  'Scenario1'
        % assumption 0: No parallel planning and execution. i.e when executing DO NOT plan and vice versa
        % assumption 1: Planning capacity :  maximum of X digits ahead
        % assumption 2: The Planning time for every X digits or less is fixed ---->   Tp =200ms
        % assumption 3: The Execution time per digit is fixed ----> Te  = 100ms
        
        % Implimentation of the model
        
        timeSteps  =1;
        for Horizon = 1:8
            for X = 1:8
                xx = min(X , Horizon);
                MT_plan(X , Horizon) = ceil(N/xx)*Tp;
                MT_exec(X , Horizon) = N*Te;
                MT(X , Horizon) = MT_plan(X , Horizon) + MT_exec(X , Horizon);
            end
        end
        
        
        for Horizon = 1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                L = N - xx*floor(N/xx);
                pln_TS =  floor(N/xx) + length(find(L));
                exe_TS = N;
                tot_TS = exe_TS + pln_TS;          % Total number of time steps
                
                TimeCourse_pln = zeros(tot_TS+1 , 4);
                TimeCourse_exe = zeros(tot_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:xx+1:tot_TS];
                plnTimes = zeros(size(log_T_pln));
                
                
                Tp_h     = Tp*exp(.15*(1-1)); % each planning time step is equal so just for consistency across code, set x1 to 1
                TimeCourse_pln(log_T_pln+1 , 2) = [1:xx:N]-1;
                TimeCourse_pln(log_T_pln+1 , 3) = TimeCourse_pln(log_T_pln+1 , 2) + xx ;
                if L
                    Tp_h_end = Tp*exp(.15*(1-1));  % planning time for the last L digits
                    plnTimes(1:end-1) = Tp_h;
                    plnTimes(end) = Tp_h_end;
                    TimeCourse_pln(log_T_pln(end)+1 , 3) = N;
                else
                    plnTimes(1:end) = Tp_h;
                end
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 2:tot_TS;
                log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                exeTimes = Te*ones(size(log_T_exe));
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(find(TimeCourse_pln(:,1)) , 4) = 1;
                TimeCourse_exe(find(TimeCourse_exe(:,1)) , 4) = 2;
                
                Movement_Time = TimeCourse_pln + TimeCourse_exe;
                Movement_Time(:,1) = cumsum(Movement_Time(:,1));
                
                A = cumsum(TimeCourse_pln(:,1) + TimeCourse_exe(:,1));
                TimeCourse_pln(:,1) = A;
                TimeCourse_exe(:,1) = A;
                
                for i = 1:xx
                    idx = find(TimeCourse_pln(:,3) == 0);
                    TimeCourse_pln(idx(2:end),3) = TimeCourse_pln(idx(2:end)-1,3);
                    
                    idx = find(TimeCourse_exe(:,3) == 0);
                    TimeCourse_exe(idx(2:end),3) = TimeCourse_exe(idx(2:end)-1,3);
                end
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
                legend({'planning time course' , 'execution time course'})
                yticks([0:14])
                ylabel('Presses')
                % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
                % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
                grid on
                
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        
        figure('color' , 'white');
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)) , [100,400]);
            yticks([1:8])
            ylabel('Horizon')
            xticks(1:14)
            xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Buffer = ' , num2str(X)])
        end
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
        
        
        
        
    case 'Scenario2_expo'
        % assumption 0: No parallel planning and execution. i.e when executing DO NOT plan and vice versa
        % assumption 1: Planning capacity :  maximum of X digits ahead
        % assumption 2: The Planning time for every single digit is fixed ---->   Tp =200ms and grows exponentially for more
        % assumption 3: After the first x1 digits, planning is done one ahead since only one at a time in revealed
        % assumption 4: The Execution time per digit is fixed ----> Te  = 100ms
        
        % Implimentation of the model
        for Horizon = 1:N
            for X = 1:N
                xx = min(X , Horizon);
                
                
                Tp_h = Tp*exp(.15*(xx-1)); % planning time needed for planning X digits   ----> Exponential growth
%                 Tp_h = Tp*(.7*x1);           % planning time needed for planning X digits   ----> Linear growth
                % if Horizon = X then Tp_h = Tp. this is the time needed to plan the firat h digits that appear
                % form the first h digits onward, the planning time in fixed because it's always one ahead
                L = N - xx*floor(N/xx);
                MT_plan(X , Horizon) = Tp_h*floor(N/xx) + Tp*exp(.15*(L-1)); % exponential growth
%                 MT_plan(X , Horizon) = Tp_h*floor(N/x1) + Tp*(.7*L);  % Linear growth
                MT_exec(X , Horizon) = N*Te;
                MT(X , Horizon) = MT_plan(X , Horizon) + MT_exec(X , Horizon);
            end
        end
        
        for Horizon = 1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                L = N - xx*floor(N/xx);
                pln_TS =  floor(N/xx) + length(find(L));
                exe_TS = N;
                tot_TS = exe_TS + pln_TS;          % Total number of time steps
                
                TimeCourse_pln = zeros(tot_TS+1 , 4);
                TimeCourse_exe = zeros(tot_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:xx+1:tot_TS];
                plnTimes = zeros(size(log_T_pln));
                
                
                 Tp_h     = Tp*exp(.15*(xx-1)); % each planning time step   exponential growth
%                 Tp_h  =  Tp*(.7*x1);           % planning time needed for planning X digits   ----> Linear growth
                TimeCourse_pln(log_T_pln+1 , 2) = [1:xx:N]-1;
                TimeCourse_pln(log_T_pln+1 , 3) = TimeCourse_pln(log_T_pln+1 , 2) + xx ;
                if L
                    Tp_h_end = Tp*exp(.15*(L-1));  % planning time for the last L digits  %exponential growth
%                     Tp_h_end = Tp*(.7*L);  % Linear growth
                    plnTimes(1:end-1) = Tp_h;
                    plnTimes(end) = Tp_h_end;
                    TimeCourse_pln(log_T_pln(end)+1 , 3) = N;
                else
                    plnTimes(1:end) = Tp_h;
                end
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 2:tot_TS;
                log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                exeTimes = Te*ones(size(log_T_exe));
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(find(TimeCourse_pln(:,1)) , 4) = 1;
                TimeCourse_exe(find(TimeCourse_exe(:,1)) , 4) = 2;
                
                Movement_Time = TimeCourse_pln + TimeCourse_exe;
                Movement_Time(:,1) = cumsum(Movement_Time(:,1));
                
                A = cumsum(TimeCourse_pln(:,1) + TimeCourse_exe(:,1));
                TimeCourse_pln(:,1) = A;
                TimeCourse_exe(:,1) = A;
                
                for i = 1:xx
                    idx = find(TimeCourse_pln(:,3) == 0);
                    TimeCourse_pln(idx(2:end),3) = TimeCourse_pln(idx(2:end)-1,3);
                    
                    idx = find(TimeCourse_exe(:,3) == 0);
                    TimeCourse_exe(idx(2:end),3) = TimeCourse_exe(idx(2:end)-1,3);
                end
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
                legend({'planning time course' , 'execution time course'})
                yticks([0:14])
                ylabel('Presses')
                % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
                % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
                grid on
                
%                 subplot(1,2,1)
%                 plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 title('Planning Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 yticks([0:14])
%                 ylabel('Presses')
%                 % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%                 % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
%                 grid on
%                 
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)));
            yticks([1:8])
            ylabel('Horizon')
            xticks(1:14)
            xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Buffer = ' , num2str(X)])
        end
        
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
        
     case 'Scenario2_Lin'
        % assumption 0: No parallel planning and execution. i.e when executing DO NOT plan and vice versa
        % assumption 1: Planning capacity :  maximum of X digits ahead
        % assumption 2: The Planning time for every single digit is fixed ---->   Tp =200ms and grows exponentially for more
        % assumption 3: After the first x1 digits, planning is done one ahead since only one at a time in revealed
        % assumption 4: The Execution time per digit is fixed ----> Te  = 100ms
        
        % Implimentation of the model
        for Horizon = 1:N
            for X = 1:N
                xx = min(X , Horizon);
                
                
                %Tp_h = Tp*exp(.15*(x1-1)); % planning time needed for planning X digits   ----> Exponential growth
                Tp_h = Tp+Tp*(.7*(xx-1));           % planning time needed for planning X digits   ----> Linear growth
                % if Horizon = X then Tp_h = Tp. this is the time needed to plan the firat h digits that appear
                % form the first h digits onward, the planning time in fixed because it's always one ahead
                L = N - xx*floor(N/xx);
                % MT_plan(X , Horizon) = Tp_h*floor(N/x1) + Tp*exp(.15*(L-1)); % exponential growth
                Tp_h_end = Tp+Tp*(.7*(L-1));  % Linear growth
                MT_plan(X , Horizon) = Tp_h*floor(N/xx) + Tp_h_end;  % Linear growth
                MT_exec(X , Horizon) = N*Te;
                MT(X , Horizon) = MT_plan(X , Horizon) + MT_exec(X , Horizon);
            end
        end
        
        for Horizon = 1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                L = N - xx*floor(N/xx);
                pln_TS =  floor(N/xx) + length(find(L));
                exe_TS = N;
                tot_TS = exe_TS + pln_TS;          % Total number of time steps
                
                TimeCourse_pln = zeros(tot_TS+1 , 4);
                TimeCourse_exe = zeros(tot_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:xx+1:tot_TS];
                plnTimes = zeros(size(log_T_pln));
                
                
                %  Tp_h     = Tp*exp(.15*(x1-1)); % each planning time step   exponential growth
                Tp_h  =  Tp+Tp*(.7*(xx-1));           % planning time needed for planning X digits   ----> Linear growth
                TimeCourse_pln(log_T_pln+1 , 2) = [1:xx:N]-1;
                TimeCourse_pln(log_T_pln+1 , 3) = TimeCourse_pln(log_T_pln+1 , 2) + xx ;
                if L
                    %  Tp_h_end = Tp*exp(.15*(L-1));  % planning time for the last L digits  %exponential growth
                    Tp_h_end = Tp+Tp*(.7*(L-1));  % Linear growth
                    plnTimes(1:end-1) = Tp_h;
                    plnTimes(end) = Tp_h_end;
                    TimeCourse_pln(log_T_pln(end)+1 , 3) = N;
                else
                    plnTimes(1:end) = Tp_h;
                end
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 2:tot_TS;
                log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                exeTimes = Te*ones(size(log_T_exe));
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(find(TimeCourse_pln(:,1)) , 4) = 1;
                TimeCourse_exe(find(TimeCourse_exe(:,1)) , 4) = 2;
                
                Movement_Time = TimeCourse_pln + TimeCourse_exe;
                Movement_Time(:,1) = cumsum(Movement_Time(:,1));
                
                A = cumsum(TimeCourse_pln(:,1) + TimeCourse_exe(:,1));
                TimeCourse_pln(:,1) = A;
                TimeCourse_exe(:,1) = A;
                
                for i = 1:xx
                    idx = find(TimeCourse_pln(:,3) == 0);
                    TimeCourse_pln(idx(2:end),3) = TimeCourse_pln(idx(2:end)-1,3);
                    
                    idx = find(TimeCourse_exe(:,3) == 0);
                    TimeCourse_exe(idx(2:end),3) = TimeCourse_exe(idx(2:end)-1,3);
                end
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
                legend({'planning time course' , 'execution time course'})
                yticks([0:14])
                ylabel('Presses')
                % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
                % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
                grid on
                
%                 subplot(1,2,1)
%                 plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 title('Planning Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 yticks([0:14])
%                 ylabel('Presses')
%                 % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%                 % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
%                 grid on
%                 
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)));
            yticks([1:8])
            ylabel('Horizon')
            xticks(1:14)
            xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Buffer = ' , num2str(X)])
        end
        
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
        
        
    case  'Scenario3'
        % assumption 0: Parallel planning and execution.
        % assumption 1: Planning and execution do not share resources. i.e. planning does not slow execution and vice versa
        % assumption 2: Planning capacity :  maximum of X digits ahead
        % assumption 3: The Planning time for every single digit is fixed ---->   Tp =200ms and grows exponentially for more
        % assumption 4: After the first x1 digits, planning is done one ahead since only one at a time in revealed
        % assumption 5: The Execution time per digit is fixed ----> Te  = 100ms
        % assumption 6: Execution can never preceed planning 
        
        
        % Implimentation of the model
        
        
        for Horizon = 1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                L = N - xx;
                pln_TS =  1 + L;     % the first x1 ones are planned together and then the rest are one at a time
                exe_TS =  N;
                
                
                TimeCourse_pln = zeros(pln_TS+1 , 4);
                TimeCourse_exe = zeros(exe_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:pln_TS];  % you dont have to wait for execution
                plnTimes = zeros(size(log_T_pln));
                
                
%                 Tp_h     = Tp*exp(.15*(x1-1)); % each planning time step
                Tp_h = Tp+Tp*(.7*(xx-1));           % planning time needed for planning X digits   ----> Linear growth
                TimeCourse_pln(log_T_pln+1 , 2) = [1,xx , xx+1:N-1]';
                TimeCourse_pln(log_T_pln+1 , 3) = xx:N ;
                
                plnTimes(1) = Tp_h;
                plnTimes(2:end) = Tp;
                
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 1:exe_TS;
                %                 log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                exeTimes = Te*ones(size(log_T_exe));
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(: , 4) = 1;
                TimeCourse_exe(: , 4) = 2;
                
                TimeCourse_pln(: , 1) = cumsum(TimeCourse_pln(: , 1));
                TimeCourse_exe(: , 1) = cumsum(TimeCourse_exe(: , 1));
                
                MT_plan(X , Horizon) = TimeCourse_pln(end,1) ;
                MT_exec(X , Horizon) = TimeCourse_exe(end,1);
                
                for press = 1:N
                    ide = find(TimeCourse_exe(:,3) == press , 1 , 'first');
                    idp = find(TimeCourse_pln(:,3) >= press , 1 , 'first');
                    if TimeCourse_exe(ide,1) < TimeCourse_pln(idp , 1)
                        difference = abs(TimeCourse_exe(ide,1) - TimeCourse_pln(idp , 1));
                        TimeCourse_exe(ide:end,1) = TimeCourse_exe(ide:end,1)+ difference;
                    end
                end
                
                MT_plan(X , Horizon) = Tp_h + L*(Tp);
                MT_exec(X , Horizon) = N*Te ;
                MT(X , Horizon) = TimeCourse_exe(end,1);
                
                
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
                legend({'planning time course' , 'execution time course'})
                yticks([0:14])
                ylabel('Presses')
                grid on
                
                
%                 subplot(1,2,1)
%                 plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 title('Planning Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 yticks([0:14])
%                 ylabel('Presses')
%                 % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%                 % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
%                 grid on
%                 
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
        
        figure('color' , 'white');
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)));
            yticks([1:8])
            ylabel('Horizon')
            xticks(1:14)
            xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Limit = ' , num2str(X)])
        end
        
        
        
        
    case 'Scenario4'
        % assumption 0: Parallel planning and execution.
        % assumption 1: Planning and execution do share resources. i.e. planning does slow execution and vice versa by the same factor ?
        % assumption 2: Planning capacity :  maximum of X digits ahead
        % assumption 3: The Planning time for every single digit is fixed ---->   Tp =200ms and grows exponentially for more
        % assumption 4: After the first x1 digits, planning is done one ahead since only one at a time in revealed
        % assumption 5: The Execution time per digit is fixed ----> Te  = 100ms
        % assumption 6: Execution can never preceed planning 
        
        
        % Implimentation of the model
     
        
        
        for Horizon = 2%1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                beta = 1.5;
                
                Tp_plan_while_exe = Tp*beta;
                Te_exe_while_plan = Te*3;
                
                L = N - xx;
                pln_TS =  1 + L;     % the first x1 ones are planned together and then the rest are one at a time
                exe_TS =  N;
                
                
                TimeCourse_pln = zeros(pln_TS+1 , 4);
                TimeCourse_exe = zeros(exe_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:pln_TS];  % you dont have to wait for execution
                plnTimes = zeros(size(log_T_pln));
                
                
                %Tp_h     = Tp * exp(.15*(x1-1)); % first planning time step   --> exponential growth
                Tp_h     = Tp + Tp*(.7*(xx-1));  % planning time needed for planning X digits   ----> Linear growth
                
                TimeCourse_pln(log_T_pln+1 , 2) = [1,xx , xx+1:N-1]';
                TimeCourse_pln(log_T_pln+1 , 3) = xx:N ;
                
                plnTimes(1) = Tp_h;
                plnTimes(2:end) = Tp_plan_while_exe;
                
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 1:exe_TS;
                %                 log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                
                exeTimes = [Te*ones(1,xx-1) Te_exe_while_plan*ones(1, length(log_T_exe)-xx) Te];
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(: , 4) = 1;
                TimeCourse_exe(: , 4) = 2;
                
                TimeCourse_pln(: , 1) = cumsum(TimeCourse_pln(: , 1));
                TimeCourse_exe(: , 1) = cumsum(TimeCourse_exe(: , 1));
                
                MT_plan(X , Horizon) = TimeCourse_pln(end,1) ;
                MT_exec(X , Horizon) = TimeCourse_exe(end,1);
                
                for press = 1:N
                    ide = find(TimeCourse_exe(:,3) == press , 1 , 'first');
                    idp = find(TimeCourse_pln(:,3) >= press , 1 , 'first');
                    if TimeCourse_exe(ide,1) < TimeCourse_pln(idp , 1)
                        difference = abs(TimeCourse_exe(ide,1) - TimeCourse_pln(idp , 1));
                        TimeCourse_exe(ide:end,1) = TimeCourse_exe(ide:end,1)+ difference;
                    end
                end
                

                MT(X , Horizon) = TimeCourse_exe(end,1);
                
                
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
                legend({'planning time course' , 'execution time course'})
                yticks([0:14])
                ylabel('Presses')
                grid on
                
                
%                 subplot(1,2,1)
%                 plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 title('Planning Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 yticks([0:14])
%                 ylabel('Presses')
%                 % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%                 % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
%                 grid on
%                 
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        figure('color' , 'white');
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)));
            yticks([1:8])
            ylabel('Horizon')
            xticks(1:14)
            xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Limit = ' , num2str(X)])
        end
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
        
        
    case 'Scenario5'
        % assumption 0: Parallel planning and execution.
        % assumption 1: Planning and execution do share resources. i.e. planning does slow execution and vice versa by the same factor ?
        % assumption 2: Planning capacity :  maximum of X digits ahead
        % assumption 3: The Planning time for every single digit is fixed ---->   Tp =200ms and grows exponentially for more
        % assumption 4: After the first x1 digits, planning is done one ahead since only one at a time in revealed
        % assumption 5: The Execution time per digit is fixed ----> Te  = 100ms
        % assumption 6: Execution can never preceed planning 
        
        % 1 - execution reduces buffer size by 1 so the planning limit will be X-1
        % 2 - if you have more transitions ahead planned, you are faster at executing
        % 3 - you are fastes at the end because you have everything ahead planned
        

        % Implimentation of the model
     
        
        
        for Horizon = 1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                gamma = 1.1;
                
                Tp_plan_while_exe = gamma*Tp;
                Te_exe_while_plan = gamma*Te;
                
                beta =  0.6667/13 ;
%                 Te_X1_planned_ahead = Te - Te*(beta*(x1-2));
                Te_X1_planned_ahead = Te + Te*exp(.3*(2-xx));
                
%                 Te_X1_oneless_planned_ahead = Te - Te*(beta*(max(0,x1-1)));
                Te_X1_oneless_planned_ahead = Te + Te*exp(.3*(2-(xx-1)));
                
%                 Te_X1_oneLess_planned_ahead_whilePlan = Te_exe_while_plan - Te_exe_while_plan*(beta*(max(0,x1-1)));
                Te_X1_oneLess_planned_ahead_whilePlan = Te_exe_while_plan + Te_exe_while_plan*exp(.3*(2-(xx-1)));
                
                pln_TS =  (N - xx)+1;     % the first x1 ones are planned together and then the rest are one at a time
                exe_TS =  N;
                
                
                TimeCourse_pln = zeros(pln_TS+1 , 4);
                TimeCourse_exe = zeros(exe_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:pln_TS];  % you dont have to wait for execution
                plnTimes = zeros(size(log_T_pln));
                
                
                %Tp_h     = Tp * exp(.15*(x1-1)); % first planning time step   --> exponential growth
                Tp_h        = Tp + Tp*(2*beta*(xx-1));  % planning time needed for planning X digits   ----> Linear growth
                theRestTp_h = Tp_plan_while_exe + Tp_plan_while_exe*(beta*(max(0,xx-1))); % this has alower coefficient cause ssome of them are alreadyplanned
                
                if xx == 1
                    TimeCourse_pln(log_T_pln+1 , 2) = [0:N-1]'; %starts of planning
                    TimeCourse_pln(log_T_pln+1 , 3) = [1:N]' ;
                    exeTimes = Te_X1_planned_ahead*ones(1,N);
                    plnTimes = Tp*ones(1,N);
                elseif xx==2
                    TimeCourse_pln(log_T_pln+1 , 2) = [0,xx:N-1]'; %starts of planning
                    TimeCourse_pln(log_T_pln+1 , 3) = [2:N]' ;
                    exeTimes = [Te_X1_planned_ahead*ones(1,min(1,xx-1)) Te_X1_oneLess_planned_ahead_whilePlan*ones(1, N-2-(xx-2)) Te_X1_oneless_planned_ahead*ones(1,max(1,xx-1))];
                    plnTimes = [Tp_h ones(1,pln_TS-1)*theRestTp_h];
                else
                    TimeCourse_pln(log_T_pln+1 , 2) = [0,3:N-(xx-2)]'; %starts of planning
                    TimeCourse_pln(log_T_pln+1 , 3) = [xx :N]' ;
                    exeTimes = [Te_X1_planned_ahead*ones(1,min(1,xx-1)) Te_X1_oneLess_planned_ahead_whilePlan*ones(1, N-2-(xx-2)) Te_X1_oneless_planned_ahead*ones(1,max(1,xx-1))];
                    plnTimes = [Tp_h ones(1,pln_TS-1)*theRestTp_h];
                end
                
                
                [X Horizon]
                
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 1:exe_TS;
                %                 log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                

                
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(: , 4) = 1;
                TimeCourse_exe(: , 4) = 2;
                TimeCourse_pln(:,1) 
                TimeCourse_exe(:,1)
                
                TimeCourse_pln(: , 1) = cumsum(TimeCourse_pln(: , 1));
                TimeCourse_exe(: , 1) = cumsum(TimeCourse_exe(: , 1));
                
                MT_plan(X , Horizon) = TimeCourse_pln(end,1) ;
                MT_exec(X , Horizon) = TimeCourse_exe(end,1);
                
                if xx ==1
                    TimeCourse_exe(:,1) = TimeCourse_exe(: , 1) + TimeCourse_pln(: , 1);
                else
                    for press = 1:N
                        ide = find(TimeCourse_exe(:,3) == press , 1 , 'first');
                        idp = find(TimeCourse_pln(:,3) >= press , 1 , 'first');
                        if TimeCourse_exe(ide,1) < TimeCourse_pln(idp , 1)
                            difference = abs(TimeCourse_exe(ide,1) - TimeCourse_pln(idp , 1));
                            TimeCourse_exe(ide:end,1) = TimeCourse_exe(ide:end,1)+ difference;
                        end
                    end
                end
                

                MT(X , Horizon) = TimeCourse_exe(end,1);
                MT(X , Horizon)
                
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                squeeze(IPI(Horizon , X ,:))'
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
                legend({'planning time course' , 'execution time course'})
%                 yticks([0:14])
                ylabel('Presses')
                grid on
                
                
%                 subplot(1,2,1)
%                 plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 title('Planning Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 yticks([0:14])
%                 ylabel('Presses')
%                 % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%                 % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
%                 grid on
%                 
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        figure('color' , 'white');
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)));
%             yticks([1:8])
            ylabel('Horizon')
%             xticks(1:14)
%             xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Limit = ' , num2str(X)])
        end
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
        
        
    case 'Scenario6'
        % assumption 0: Parallel planning and execution.
        % assumption 1: Planning and execution do share resources. i.e. planning does slow execution and vice versa by the same factor ?
        % assumption 2: Planning capacity :  maximum of X digits ahead
        % assumption 3: The Planning time for every single digit is fixed ---->   Tp =200ms and grows exponentially for more
        % assumption 4: After the first x1 digits, planning is done one ahead since only one at a time in revealed
        % assumption 5: The Execution time per digit is fixed ----> Te  = 100ms
        % assumption 6: Execution can never preceed planning 
        
        % 1 - execution reduces buffer size by 1 so the planning limit will be X-1
        % 2 - if you have more transitions ahead planned, you are faster at executing
        % 3 - you are fastes at the end because you have everything ahead planned
        

        % Implimentation of the model
        
        
        
        for Horizon = 1:8
            figure('color' , 'white');
            title(['Horizon' , num2str(Horizon)])
            for X = 1:8
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                xx = min(X , Horizon);
                
                gamma = 1.1; % the factor by which the planning time increases if sonr in synchrony with execution
                
                
                Tp_plan_while_exe = gamma*Tp;
                Te_exe_while_plan = 5*gamma*Te;
                
                beta =  0.6667/13 ;     %  the coefficient by which planning time increases when planning more than 1
                Te_plan_ahead_Coef = 2; % the factor used for the exponential growth of exection time the less you have planned
%                 Te_X1_planned_ahead = Te - Te*(beta*(x1-2));
                Te_X1_planned_ahead = Te + Te*exp(Te_plan_ahead_Coef*(2-xx));
                
%                 Te_X1_oneless_planned_ahead = Te - Te*(beta*(max(0,x1-1)));
                Te_X1_oneless_planned_ahead = Te + Te*exp(Te_plan_ahead_Coef*(2-(xx-1)));
                
%                 Te_X1_oneLess_planned_ahead_whilePlan = Te_exe_while_plan - Te_exe_while_plan*(beta*(max(0,x1-1)));
                Te_X1_oneLess_planned_ahead_whilePlan = Te_exe_while_plan + Te_exe_while_plan*exp(Te_plan_ahead_Coef*(2-(xx-1)));
                
                if xx==1
                    pln_TS = N;
                elseif xx==2
                    pln_TS = N-1;
                else
                    pln_TS =  max(N-2*xx+4 , 2);     % the first x1 ones are planned together and then the rest are one at a time
                end
                exe_TS =  N;
                
                TimeCourse_pln = zeros(pln_TS+1 , 4);
                TimeCourse_exe = zeros(exe_TS+1 , 4);
                
                % setting up planning time course
                log_T_pln = [1:pln_TS];  % you dont have to wait for execution
                plnTimes = zeros(size(log_T_pln));
                
                
                %Tp_h     = Tp * exp(.15*(x1-1)); % first planning time step   --> exponential growth
                Tp_h        = Tp + Tp*(beta*(xx-1));  % planning time needed for planning X digits   ----> Linear growth
                theRestTp_h = Tp_plan_while_exe + Tp_plan_while_exe*(beta*(max(0,xx-1))); % this has alower coefficient cause ssome of them are alreadyplanned
                
                if xx == 1
                    TimeCourse_pln(log_T_pln+1 , 2) = [0:N-1]'; %starts of planning
                    TimeCourse_pln(log_T_pln+1 , 3) = [1:N]' ;
                    exeTimes = Te_X1_oneless_planned_ahead*ones(1,N);
                    plnTimes = Tp*ones(1,N);
                elseif xx==2
                    TimeCourse_pln(log_T_pln+1 , 2) = [0,xx:N-1]'; %starts of planning
                    TimeCourse_pln(log_T_pln+1 , 3) = [2:N]' ;
                    exeTimes = [Te_X1_oneless_planned_ahead*ones(1,min(1,xx-1)) Te_X1_oneLess_planned_ahead_whilePlan*ones(1, N-2) Te_X1_oneless_planned_ahead];
                    plnTimes = [Tp_h ones(1,pln_TS-1)*theRestTp_h];
                else
                    TimeCourse_pln(log_T_pln+1 , 2) = [0,xx:N-(xx-2)]'; %starts of planning
                    TimeCourse_pln(log_T_pln+1 , 3) = [xx,2*xx-2:N]' ;
                    exeTimes = [];
                    % the execution time grows exponentially the less that you have planned
                    for xx = xx-1:-1:2
                        exeTimes = [exeTimes Te + Te*exp(Te_plan_ahead_Coef*(2-xx))];
                    end
                    xx = xx-1;
     
                    exeTimes = [exeTimes (Te_exe_while_plan + Te_exe_while_plan*exp(Te_plan_ahead_Coef*(2-xx)))*ones(1,pln_TS-1)];
                    
                    for xx = 1:xx-1
                        exeTimes = [exeTimes Te + Te*exp(Te_plan_ahead_Coef*(2-xx))];
                    end
                   
                    plnTimes = [Tp_h ones(1,pln_TS-1)*theRestTp_h];
                end
                
                
                [X Horizon]
                
                
                TimeCourse_pln(log_T_pln+1 , 1) = plnTimes;
                
                % setting up exe time course
                log_T_exe = 1:exe_TS;
                %                 log_T_exe = log_T_exe(~ismember(log_T_exe,log_T_pln));
                

                
                
                TimeCourse_exe(log_T_exe+1 , 2) = [1:N] -1;
                TimeCourse_exe(log_T_exe+1 , 3) = [1:N];
                
                
                TimeCourse_exe(log_T_exe+1 , 1) = exeTimes;
                
                % Mixing up planning and exe time
                TimeCourse_pln(: , 4) = 1;
                TimeCourse_exe(: , 4) = 2;
                Plan = TimeCourse_pln(:,1) 
                Exe  = TimeCourse_exe(:,1)
                
                TimeCourse_pln(: , 1) = cumsum(TimeCourse_pln(: , 1));
                TimeCourse_exe(: , 1) = cumsum(TimeCourse_exe(: , 1));
                
                MT_plan(X , Horizon) = TimeCourse_pln(end,1) ;
                MT_exec(X , Horizon) = TimeCourse_exe(end,1);
                
                if xx ==1
                    TimeCourse_exe(:,1) = TimeCourse_exe(: , 1) + TimeCourse_pln(: , 1);
                else
                    for press = 1:N
                        ide = find(TimeCourse_exe(:,3) == press , 1 , 'first');
                        idp = find(TimeCourse_pln(:,3) >= press , 1 , 'first');
                        if TimeCourse_exe(ide,1) < TimeCourse_pln(idp , 1)
                            difference = abs(TimeCourse_exe(ide,1) - TimeCourse_pln(idp , 1));
                            TimeCourse_exe(ide:end,1) = TimeCourse_exe(ide:end,1)+ difference;
                        end
                    end
                end
                

                MT(X , Horizon) = TimeCourse_exe(end,1);
                MT(X , Horizon)
                
                
                IPI(Horizon , X , 1) = TimeCourse_exe(2,1);
                for press = 1:N-1
                    t1 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press , 1 , 'first'),1);
                    t2 = TimeCourse_exe(find(TimeCourse_exe(:,3) == press+1 , 1 , 'first'),1);
                    IPI(Horizon , X , press+1) = t2 - t1;
                end
                squeeze(IPI(Horizon , X ,:))'
                
                subplot(2,4,X)
                plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                hold on
                plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 );%, 'color' , colors(X,:))
                title(['Planning Buffer Size = ' , num2str(X)]);
                legend({'planning time course' , 'execution time course'})
%                 yticks([0:14])
                ylabel('Presses')
                grid on
                
                
%                 subplot(1,2,1)
%                 plot(TimeCourse_pln(:,1) , TimeCourse_pln(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 title('Planning Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 yticks([0:14])
%                 ylabel('Presses')
%                 % xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
%                 % xticks([-1 -0.8 -0.2 0 0.2 0.8 1])
%                 grid on
%                 
%                 subplot(1,2,2)
%                 plot(TimeCourse_exe(:,1) , TimeCourse_exe(:,3) , 'LineWidth' , 3 , 'color' , colors(X,:))
%                 hold on
%                 yticks([0:14])
%                 ylabel('Presses')
%                 title('Execution Time Course')
%                 legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8'})
%                 grid on
                
            end
            
        end
        figure('color' , 'white');
        for X = 1:8
            subplot(4,2,X)
            imagesc(squeeze(IPI(:,X,:)));
%             yticks([1:8])
            ylabel('Horizon')
%             xticks(1:14)
%             xticklabels({'RT','IPI1' , 'IPI2', 'IPI3', 'IPI4', 'IPI5', 'IPI6', 'IPI7', 'IPI8', 'IPI9', 'IPI10', 'IPI11', 'IPI12', 'IPI13'})
            xlabel('IPIs')
            ylabel('Horizons')
            title(['Planning Limit = ' , num2str(X)])
        end
        
        figure('color' , 'white');
        title('IPI traces')
        for h = 1:8
            plot(squeeze(IPI(h,5,2:14)) , 'LineWidth' , 3 , 'color' , colors(h,:))
            hold on
        end
        grid on
        legend({'Horizon = 1','Horizon = 2','Horizon = 3','Horizon = 4','Horizon = 5','Horizon = 6','Horizon = 7','Horizon = 8'})
        title('planning buffer  = 5')
    case 'Scenario7'
        % assumption 1: you can only execute a digit that has been a 100% planned
        % assumption 2: execution pushes your planning percentages one down the exponential ladder
        % assumtion  3: you only worry about the xx digits ahead that you can plen/see. Even if you see more, or can plan more does not matter 
        Tp = 200; % time you need to plan the 100% of 1 press
        Te = 500; % time of executing onepress when you have one press ahead planned
        beta =  0.1 ;
        gamma = 1.1; % the factor by which the planning time increases if sonr in synchrony with execution
        Tp_plan_while_exe = gamma*Tp;
%         Te_exe_while_plan = 10*gamma*Te;
        
        
        for Horizon = 1:14
%             figure('color' , 'white');
%             title(['Horizon' , num2str(Horizon)])
            for X = 1:13
                
                xx = min(X , Horizon);
                syms Digit_start Digit_mid PlanPerc plantime planNum percPlanning PlanPerUnit plantimePort ExePerUnit DIFPLAN
                Plan_Decay_start = 28*exp(2-Digit_start);  % happen at the beggining, so no pressing
                Plan_Decay_mid   = 28*exp(2-Digit_mid);     % happen in the middle while pressing
                Te_X1_planned_ahead = (ExePerUnit + exp(4*DIFPLAN/100)*ExePerUnit)/2;
                Te_X1_planned_ahead = ExePerUnit + (DIFPLAN/100)*ExePerUnit*exp(2*(2-xx));
                Te_exe_while_plan   = Te + Te*(DIFPLAN/100);
                % Careful that plantime can be relating to different percentages of a planning a press
                plantimePortion         =  (2.7-exp(.01*(100-PlanPerc)))/1.7 ;   % the symbolic PlanPerc is the portion of 100 that you are planning a digit
                % plantime*PlanPerUnit is the time you need to plan any given percentage of a number
                % PlanPerUnit is either Tp (if not executing) or Tp_plan_while_exe (if executing)
                % plantimePort is the output of plantimePortion symbolic function
                
                Tp_moreThanOne   = (plantimePort * PlanPerUnit) + plantimePort*PlanPerUnit*(beta*(planNum-1));   
                
                clear log_T_pln log_T_exe TimeCourse_pln TimeCourse_exe time
                
                alpha = .5; % the percentage of presses in the horizon that you fully plan before execution
                full_plan_lim = ceil(alpha*xx);
                prep_perc = zeros(1,N);
                prep_perc(1:full_plan_lim) =  100;
                rem_plan = xx - full_plan_lim; 
                Digit_start = 1:rem_plan;
                prep_perc(full_plan_lim+1:full_plan_lim+rem_plan) = eval(Plan_Decay_start);
                
                MT(X , Horizon).com_plan = zeros(N , xx);
                MT(X , Horizon).prs_numb = zeros(N , xx);
                MT(X , Horizon).dif_plan = zeros(N , xx);
                MT(X , Horizon).PlanTimePerPress = zeros(N , xx);
                
                for press = 1:N
                    press
                    id = min(xx , N-press+1);
                    switch press
                        case 1  % this is essentially reaction time
                            ToBePlanned = prep_perc(press:press+(id-1));
                            MT(X , Horizon).prs_numb(press ,1 :id) = press : press+id-1; % the presses you're planning
                            planNum = id;
                            PlanPerUnit = Tp;
                            tempid =  [];
                            midl_plan_perc = 100;
                            for horz = 2:xx
                                Digit_start = horz;
                                midl_plan_perc(horz) = eval(Plan_Decay_start); % the planning top-up that you can add to the previous planned ones
                            end
                        otherwise
                            end_id = 2:size(MT(X , Horizon).com_plan , 2);
                            MT(X , Horizon).prs_numb(press ,1 :id) = press : press+id-1;
                            if ~isempty(end_id)
                                ToBePlanned = [MT(X , Horizon).com_plan(press-1 , end_id) prep_perc(press + id-1)];
                            else
                                ToBePlanned = [prep_perc(press)];
                            end
                            
                            planNum = id;
                            PlanPerUnit = Tp_plan_while_exe;
                            tempid =  ToBePlanned == 100;
                            midl_plan_perc = 50;
                            for horz = 2:xx
                                Digit_mid = horz;
                                midl_plan_perc(horz) = eval(Plan_Decay_mid); % the planning top-up that you can add to the previous planned ones
                            end
                            
                    end
                    MT(X , Horizon).plan_Done(press,1:id)  = midl_plan_perc(1:id);
                    [X , Horizon]
                    MT(X , Horizon).plan_Done(press,:) = [MT(X , Horizon).plan_Done(press,1:id) zeros(1,xx-id)] + ToBePlanned; % MT(X , Horizon).plan_Done(press,:) is the planned movements after the planning limits while pressing have been added
                    % if after adding the execc plannig, a certian press exceeds 100, divide the axxecc planning between the other
                    temp3 = (MT(X , Horizon).plan_Done(press,:) > 100); % 
                    if sum(temp3)
                       temp4 = sum(MT(X , Horizon).plan_Done(press , temp3))  - 100*sum(temp3);
                       MT(X , Horizon).plan_Done(press,temp3) = 100;
                       temp5 = (MT(X , Horizon).plan_Done(press,:) < 100 & MT(X , Horizon).plan_Done(press,:)~=0);
                    end
                    temp3 = (MT(X , Horizon).plan_Done(press,:) > 100); % 
                    if sum(temp3)
                       temp4 = sum(MT(X , Horizon).plan_Done(press,temp3))  - 100*sum(temp3);
                       MT(X , Horizon).plan_Done(press,temp3) = 100;
                       temp5 = (MT(X , Horizon).plan_Done(press,:) < 100 & MT(X , Horizon).plan_Done(press,:)~=0);
                    end
                    % if the immediate press that you have to d o right after the press, is not fully prepared, you have prepare it fully then
                    
                    MT(X , Horizon).com_plan(press , :)        = MT(X , Horizon).plan_Done(press,:);
                    if press>1
                        MT(X , Horizon).dif_plan(press , :)           = MT(X , Horizon).plan_Done(press,:) - ToBePlanned;
                    else
                        MT(X , Horizon).dif_plan(press , :) = MT(X , Horizon).plan_Done(press,:); %right before the first press you are planning everything from scratch
                    end
                    PlanPerc                      = MT(X , Horizon).dif_plan(press , :);
%                   plantimePort                  = eval(plantimePortion);
                    plantimePort                = PlanPerc/100;
                    plantimePort(MT(X , Horizon).dif_plan(press , :) == 0) = 0;
                    MT(X , Horizon).PlanTimePerPress(press , :) =  eval(Tp_moreThanOne); % this is as if your planning each of them this way
%                     MT(X , Horizon).PlanTimePerPress(press , :) =  MT(X , Horizon).PlanTimePerPress(press , :)/id;
                    
                    if MT(X , Horizon).plan_Done(press,1) < 100
                        PlanPerUnit = Tp; % because the execution has to wait for the  prep, so you
                        PlanPerc = 100*ones(size(MT(X , Horizon).plan_Done(press,:))) - MT(X , Horizon).plan_Done(press,:);
                        plantimePort     = eval(plantimePortion);
                        MT(X , Horizon).PlanTimePerPress(press , :) = MT(X , Horizon).PlanTimePerPress(press , :) + eval(Tp_moreThanOne); % this is as if your planning each of them this way
                    end
                    
                end
                
                for press = 1:N-1
                    if sum(MT(X , Horizon).PlanTimePerPress(press+1,:) , 2)
                        DIFPLAN  = sum(MT(X , Horizon).dif_plan(press+1 , :));
                        ExePerUnit = eval(Te_exe_while_plan);
                    else
                        ExePerUnit = Te;
                        DIFPLAN  = sum(MT(X , Horizon).dif_plan(press+1 , :));
                    end
                    MT(X , Horizon).Exec_Time(press) = eval(Te_X1_planned_ahead);
                end
                press = N;
                ExePerUnit = Te;
                DIFPLAN  = 0;
                MT(X , Horizon).Exec_Time(press) = eval(Te_X1_planned_ahead);
                
                MT(X , Horizon).Reac_Time = sum(MT(X , Horizon).PlanTimePerPress(1,:));
                MT(X , Horizon).Plan_Time = sum(MT(X , Horizon).PlanTimePerPress(2:end,:) , 2);
                MT(X , Horizon).Exec_Time(MT(X , Horizon).Plan_Time == 0) = Te;
                a = [0 ;cumsum(MT(X , Horizon).Plan_Time)]';
                b = cumsum(MT(X , Horizon).Exec_Time);
                for press = 1:N
                    temp(press) = b(press);%max(a(press) , b(press));
                end
                MT(X , Horizon).Move_Time = temp(end);
                MT(X , Horizon).IPI = diff(temp);
            end
        end
end


% ********************  Reaction time
figure('color' , 'white');
for h = 1:14
    for x = 1:8
        RT(x,h) = MT(x , h).Reac_Time; 
    end
end
for x = 1:8
    plot(RT(x ,:) ,'color',colors(x,:), 'LineWidth' , 3)
    hold on
end
title(['Initial reaction time for different buffer sizes vs. Horizons'])
ylabel('msec'  ,'FontSize' , 20)
xlabel('Horizon')
ax.XTick = [1:14];
legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8', 'X = 9', 'X = 10', 'X = 11', 'X = 12', 'X = 13', 'X = 14'})
set(gca ,'Box' , 'off','FontSize' , 20)
grid on
% ******************** Planning time
figure('color' , 'white');
for h = 14
    for x = 1:8
        PT(x, :) = MT(x , h).Plan_Time; 
    end
end
for x = 1:8
    plot(PT(x ,:) ,'color',colors(x,:), 'LineWidth' , 3)
    hold on
end
title(['Planning time per press for different buffer sizes vs. Horizons'])
ylabel('msec'  ,'FontSize' , 20)
xlabel('Presses')
ax.XTick = [1:14];
legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8', 'X = 9', 'X = 10', 'X = 11', 'X = 12', 'X = 13', 'X = 14'})
set(gca ,'Box' , 'off','FontSize' , 20)
grid on
% ******************* Movement time
figure('color' , 'white');
for h = 1:14
    for x = 1:8
        mt(x,h) = MT(x , h).Move_Time; 
    end
end
for x = 1:8
    plot(mt(x ,:) ,'color',colors(x,:), 'LineWidth' , 3)
    hold on
end
title(['Movement time time per press for different buffer sizes vs. Horizons'])
ylabel('msec' )
xlabel('Presses')
ax.XTick = [1:14];
set(gca ,'Box' , 'off','FontSize' , 20)
legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8', 'X = 9', 'X = 10', 'X = 11', 'X = 12', 'X = 13', 'X = 14'})
grid on

% ******************* IPI
figure('color' , 'white');
% buffer of interest
x = 5;
for h = 1:14
    IPI(h , :) = MT(x , h).IPI    ;
end
subplot(1,2,1)
for h = 1:14
    plot(IPI(h ,:) ,'color',colors(h,:), 'LineWidth' , 3)
    hold on
end
title(['IPIs for biffer size 5 in different Horizons'])
ylabel('msec'  ,'FontSize' , 20)
xlabel('Presses')
ax.XTick = [1:14];
legend({'X = 1' , 'X = 2' , 'X = 3', 'X = 4', 'X = 5', 'X = 6', 'X = 7', 'X = 8', 'X = 9', 'X = 10', 'X = 11', 'X = 12', 'X = 13', 'X = 14'})
grid on

subplot(1,2,2)
imagesc(IPI)
title('IPIs for biffer size 5 in different Horizons')
ylabel('Horizon size')
xlabel('IPI')
set(gca , 'XTick'  , [1:13] , 'YTick' , [1:14] , 'YTickLabels' ,'Box' , 'off','FontSize' , 20)
end



