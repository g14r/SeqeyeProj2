function ANA=se2_subj(subjname,fig,block,trial)



% ------------------------- Directories -----------------------------------
baseDir         ='/Volumes/MotorControl/data/SeqEye2/SEp/data';
% -------------------------------------------------------------------------
cd (baseDir)
mkdir('analyze');

if nargin<2
    fig=0;
end;
datafilename = ['SEp_' subjname '.dat'];

ANA=[];

D=dload(datafilename);
if (nargin<3)
    % ------ if not specified read all blocks and trials
    block = unique(D.BN);
    outfilename  = ['analyze/se1_' subjname '.mat'];
    trial  = [];
elseif (nargin<4)
    trial = [];
    outfilename  = ['analyze/se1_' subjname '_B',num2str(block),'.mat'];
else
    outfilename  = ['analyze/se1_' subjname '_B',num2str(block),'_T',num2str(trial),'.mat'];
    
end;
trcount = 1;
%define  number of trials
for b  = 1: length(block)
    disp(['Reading Block ' , num2str(block(b))])
    MOV   = movload(['SEp_' subjname '_' num2str(block(b),'%02d') '.mov']); % all trials of a block  
    if isempty(trial)
        trial = D.TN(find(D.BN==block(b)));
    end
    
    for i=1:length(trial) % loop over all trials
        [C]=se2_trial(MOV{D.TN(trial(i))},getrow(D,(D.BN == block(b) & D.TN == trial(i))),fig);
        if size(MOV , 2) == 14 % the first 3 subjects where I was not reading the PPDs
            C.state{1,:} = MOV{D.TN(trial(i))}(:,1);
            C.xEye{1,:}  = MOV{D.TN(trial(i))}(:,12);
            C.yEye{1,:}  = MOV{D.TN(trial(i))}(:,13);
            C.Pupil{1,:}  = MOV{D.TN(trial(i))}(:,14);
        else
            C.state{1,:}  = MOV{D.TN(trial(i))}(:,1);
            C.xEye{1,:}   = MOV{D.TN(trial(i))}(:,10);
            C.yEye{1,:}   = MOV{D.TN(trial(i))}(:,11);
            C.PPDx{1,:}   = MOV{D.TN(trial(i))}(:,12);
            C.PPDy{1,:}   = MOV{D.TN(trial(i))}(:,13);
            C.Pupil{1,:}  = MOV{D.TN(trial(i))}(:,14);
        end
        
        fprintf('%d %d\n',block(b),D.TN(trial(i)));
        ANA=addstruct(ANA,C);
    end

    trial  = [];
end

