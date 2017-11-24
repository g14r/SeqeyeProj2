function N = se2_timeNorm(Dall , what)

baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';
%baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/se1_data/analyze';
subj_name = {'AT1' , 'CG1' , 'HB1' , 'JT1' , 'CB1' , 'YM1' , 'NL1' , 'SR1' , 'IB1' , 'MZ1' , 'DW1'};
load([baseDir , '/CMB_34_1.mat'])
CMB = CMB_34_1;
Days  = {1 ,2 ,3 ,4 ,5,[1:5] ,[2:5] [2:3] [4:5] };
window  = 10;
N.BN = [];
N.TN = [];
N.Horizon = [];
N.seqNumb =[];


for sub = 1:length(subj_name)
    for d  = 1:length(Days)
        for h = [1:8 , 13]
            ANA = getrow(Dall , Dall.Horizon == h & Dall.SN == sub & ismember(Dall.Day , Days{d}) & ~Dall.isError & Dall.isgood & ismember(Dall.seqNumb , [0:6]));
            %%              first calculate then time normalize the velocities
            if ~isempty(ANA.Horizon)
                for s = 0:2
                    N_temp.Horizon = [];
                    [s h d sub]
                    id = find(ANA.seqNumb == s);
                    counter = 1;
                    for i = 1:length (id)
                        if sum(isnan(ANA.xEyePosDigit{id(i)})) < .4*length(ANA.xEyePosDigit{id(i)})
                            N_temp.Horizon(counter,1) = ANA.Horizon(id(i));
                            N_temp.seqNumb(counter,1) = ANA.seqNumb(id(i));
                            N_temp.BN(counter,1) = ANA.BN(id(i));
                            N_temp.TN(counter,1) = ANA.TN(id(i));
                            idd   = linspace(1 , ANA.AllPressIdx(id(i),14)-ANA.AllPressIdx(id(i),1)+1 , 1000);
                            idde = floor(linspace(ANA.AllPressIdx(id(i),1) , ANA.AllPressIdx(id(i),14) , 1001));
                            for j = 1:length(idde)-1
                                temp1(j) = nanmedian(ANA.xEyePosDigit{id(i)}(idde(j) : idde(j+1)));
                                temp2(j) = nanmedian(ANA.PressTimeSeries{id(i)}(idde(j) : idde(j+1)));
                                temp3(j) = nanmedian(ANA.xEyeAngVelocity{id(i)}(idde(j) : idde(j+1)));
                                temp4(j) = nanmedian(ANA.pressVelocity{id(i)}(idde(j) : idde(j+1)));
                            end
                            N_temp.eye(counter , :)      = inpaint_nans(temp1,3);
                            N_temp.press(counter , :)    = inpaint_nans(temp2);
                            N_temp.eyeveloc(counter , :) = inpaint_nans(temp3);
                            N_temp.prsveloc(counter , :) = inpaint_nans(temp4);
                            counter = counter + 1;
                        end
                    end
                    if length(N_temp.Horizon) > 1
                        N = addstruct(N , N_temp);
                        clear N_temp
                    else
                        clear N_temp
                    end
                end
            end
        end
    end
end
