function se2_targetfile_mixed(SubCode, GroupCode,genchunks, CMB)


% genchunks = 0;
% WantTheStars = 0;
% SubCode = 'SZ';
% GroupCode= 1;
% provode CMB from the subject's target file folder only if you already have the chunks and need to  produce more sequences
% baseDir = '/Users/nkordjazi/Documents/SeqEye/se1/SeqEye1/TargetFiles';
baseDir = '/Users/nkordjazi/Documents/SeqEye/SeqEye2/TargetFiles';
cd(baseDir)
Fname  = [SubCode , num2str(GroupCode)  , '_tgtFiles'];
mkdir(Fname)


if genchunks
    CMB = se2_getChunks(4,0,3,3, 0);
end
if GroupCode == 1
    CMB = rmfield(CMB, 'Group2');
    CMB.Group = CMB.Group1;   % just to have consistent filed name across groups
    CMB = rmfield(CMB, 'Group1');
    for i  = 1:length(CMB.Group)
        temp = CMB.DesiredCnhkargmnt(CMB.Group(i),:);
        temp1 = [];
        for j = 1:length(find(temp))
            temp1 = [temp1 1:temp(j)];
        end
        CMB.Seq2Chunk(i , :) = temp1;
    end
elseif GroupCode == 2
    CMB = rmfield(CMB, 'Group1');
    CMB.Group = CMB.Group2;   % just to have consistent filed name across groups
    CMB = rmfield(CMB, 'Group2');
    for i  = 1:length(CMB.Group)
        temp = CMB.DesiredCnhkargmnt(CMB.Group(i),:);
        temp1 = [];
        for j = 1:length(find(temp))
            temp1 = [temp1 1:temp(j)];
        end
        CMB.Seq2Chunk(i , :) = temp1;
    end
end
save([Fname , '/' , SubCode , num2str(GroupCode) , '_CMB'] , 'CMB')
MaxPress = size(CMB.Seq2Chunk ,2);

% make target file
% switch what
%     case 'targetfile'

rng('shuffle');


NumofChunks = 6;
SequenceLength = size(CMB.Seq2Chunk , 2);
Chunks = cell2mat(CMB.Chunks);
ChunkNumb = [103 203 303 104 204 304];


Cstar = {'' ,'***' , '****'};
repeating  = 0; % determins that every sequence in the CLAT and Intermixed blocks happens twice

FT = 1:6;

OrderFields = {'seqNumb','FT','press1','press2','press3','press4','press5','press6','press7','press8','press9','press10','press11','press12','press13','press14','hand','cueS','cueC','cueP','iti','sounds' , 'Horizon' , 'StimTimeLim'};



%% Chunk Training Block
%WantTheStars = 0;
Trials = 1:60;

%             NumStarTrials = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 9 9 9];
%             OrderFields = {'seqNumb','FT','press1','press2','press3','press4','press5','press6','press7','press8','press9','press10','press11','press12','press13','press14','hand','cueS','cueC','cueP','iti','sounds'};

clear ChunkTrain
StimTimeLim = zeros(length(Trials) , 1);
if ~ repeating
    for e = 1:15 % number of Blocks
        ChunkTrain.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkTrain.FT = 2 * ones(length(Trials),1);
        ChunkTrain.iti(1:length(Trials),:) = 500;
        ChunkTrain.hand(1:length(Trials),:) = 2;
        ChunkTrain.sounds(1:length(Trials),:) = 1;
        ChunkTrain.StimTimeLim = StimTimeLim;
        
        % making sure that all the chunks are ocurring equal times in each chunk trianing block
        X = [];
        for rep = 1: length(Trials)/NumofChunks
            X = [X  randperm(NumofChunks)];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible chunks in one run is not detected
        %                 X = kron(X, ones(NumStarTrials(e) + 1 , 1));
        %                 X = X(1:length(Trials));
        % Indeices in where the stras should apear as 0
        %                 starInd = repmat([1 ; zeros(NumStarTrials(e) ,1)] , length(Trials)/(NumStarTrials(e)+1) ,1);
        for i = 1:length(Trials)
            Numdigs = sum((Chunks(X(i),:) ~= 0 ));
            ChunkTrain.Horizon(i,:) = Numdigs - 1;
            ChunkTrain.cueP{i,1} = char(regexprep(cellstr(num2str(Chunks(X(i),1:Numdigs))),'\s',''));
            
            for press= 1:Numdigs
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = Chunks(X(i),press);'];
                eval(comnd);
            end
            for press= Numdigs + 1 : SequenceLength
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = 0;'];
                eval(comnd);
            end
            % Chunk length arrangement number
            ChunkTrain.seqNumb(i,1) = ChunkNumb(X(i));
            ChunkTrain.cueC(i,:) = {'£'} ;
        end
        name = [Fname ,'/'  , SubCode ,num2str(GroupCode), '_CT_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(ChunkTrain,OrderFields));
        clear x name
        
        clear ChunkTrain
        
    end
else
    Trials = 1: 60/(repeating+1);
    for e = 1:15 % number of Blocks
        ChunkTrain.cueS = cellstr(repmat('£',length(Trials),1));
        ChunkTrain.FT = 2 * ones(length(Trials),1);
        ChunkTrain.iti(1:length(Trials),:) = 500;
        ChunkTrain.hand(1:length(Trials),:) = 2;
        ChunkTrain.sounds(1:length(Trials),:) = 1;
        ChunkTrain.StimTimeLim = StimTimeLim;
        
        % making sure that all the chunks are ocurring equal times in each chunk trianing block
        X = [];
        for rep = 1: length(Trials)/NumofChunks
            X = [X  randperm(NumofChunks)];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible chunks in one run is not detected
        %                 X = kron(X, ones(NumStarTrials(e) + 1 , 1));
        %                 X = X(1:length(Trials));
        % Indeices in where the stras should apear as 0
        %                 starInd = repmat([1 ; zeros(NumStarTrials(e) ,1)] , length(Trials)/(NumStarTrials(e)+1) ,1);
        for i = 1:length(Trials)
            Numdigs = sum((Chunks(X(i),:) ~= 0 ));
            ChunkTrain.Horizon(i,:) = Numdigs - 1;
            ChunkTrain.cueP{i,1} = char(regexprep(cellstr(num2str(Chunks(X(i),1:Numdigs))),'\s',''));
            
            for press= 1:Numdigs
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = Chunks(X(i),press);'];
                eval(comnd);
            end
            for press= Numdigs + 1 : SequenceLength
                comnd  = [' ChunkTrain.press' , num2str(press) , '(i,1) = 0;'];
                eval(comnd);
            end
            % Chunk length arrangement number
            ChunkTrain.seqNumb(i,1) = ChunkNumb(X(i));
            ChunkTrain.cueC(i,:) = {'£'} ;
        end
        repChunkTrain = getrow(ChunkTrain ,  ChunkTrain.seqNumb == 0);
        for i = Trials
            tempChunkTrain = getrow(ChunkTrain , i);
            for repp = 1:repeating+1
                repChunkTrain  = addstruct(repChunkTrain , tempChunkTrain);
            end
        end
        
        name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CT_B' , num2str(e) , '.tgt'];
        dsave(name,orderfields(repChunkTrain,OrderFields));
        clear x name
        
        clear ChunkTrain repChunkTrain
        
    end
end

%% Chunk Arrangemnet Training Block
tempCL = rem(ChunkNumb , 100);

CLA = [1 ;2];
CMB.DesiredCnhkargmnt = [4 4 3 3 ; 3 4 3 4];
Trials = 1:2;
StimTimeLim = zeros(length(Trials) , 1);
%Generate 20 Full Horizons, and then 4 of each for 3 : 12
Horizon = [1:8 , 13];
clear ChunkArrangeLearn
for BN = 1:20
    for e = 1:length(Horizon)
        A(e).ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        for rep = 1: length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        A(e).ChunkArrangeLearn.seqNumb = X;
        A(e).ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = cell(1,length(find(ChunkArrange)) );
            UtempCL = unique(tempCL);
            for j = 1: length(UtempCL)
                numCL = sum(ChunkArrange == UtempCL(j));
                ch_id = find(ChunkArrange == UtempCL(j));
                chnkind{j} = [1 1];
                while diff(chnkind{j}) == 0
                    chnkind{j} = randi(sum(tempCL == UtempCL(j)),numCL,1);
                end
                for n = 1:numCL
                    seq{ch_id(n)} = CMB.Chunks{UtempCL(j)}(chnkind{j}(n) , 1:UtempCL(j));
                end
            end
            seq = cell2mat(seq);
            A(e).ChunkArrangeLearn.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            A(e).ChunkArrangeLearn.cueC(i,:) ={'£'};
            A(e).ChunkArrangeLearn.cueS(i,:) ={'£'};
            A(e).ChunkArrangeLearn.iti(1:i,:) = 500;
            A(e).ChunkArrangeLearn.hand(i,:) = 2;
            A(e).ChunkArrangeLearn.sounds(i,:) = 1;
            A(e).ChunkArrangeLearn.Horizon(i,:) = Horizon(e);
            for press= 1:MaxPress
                comnd  = [' A(e).ChunkArrangeLearn.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        clear x
    end
    ChunkArrangeLearn = addstruct(A(1).ChunkArrangeLearn , A(2).ChunkArrangeLearn);
    for k = 3:length(A)
        ChunkArrangeLearn = addstruct(ChunkArrangeLearn , A(k).ChunkArrangeLearn);
    end
    clear A
    idxAll = randperm(length(ChunkArrangeLearn.FT));
    for kk = 1:length(OrderFields)
        statement = ['ChunkArrangeLearn.' , OrderFields{kk},' = ','ChunkArrangeLearn.' , OrderFields{kk},'(idxAll);'];
        eval(statement);
    end
    
    name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_CLAT', '_B' , num2str(BN) , '.tgt'];
    dsave(name,orderfields(ChunkArrangeLearn,OrderFields));
end







%% Test blocks - Random Sequences
clear RandomSeq
Horizon = [1:8 , 13];

ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
Trials = 1:3;
StimTimeLim = zeros(length(Trials) , 1);

for BN = 1:20
    for e = 1:length(Horizon)
        A(e).RandomSeq.StimTimeLim = StimTimeLim;
        A(e).RandomSeq.seqNumb(1:length(Trials),1) = 0;
        A(e).RandomSeq.FT(1:length(Trials),1) = 2;
        A(e).RandomSeq.Horizon(1:length(Trials),:) = Horizon(e);
        for i = 1:length(Trials)
            if ElimChunkStart
                seq = Chunks(1,1:2); % makes sure to enter the first while loop
                seq = [seq  1 1 1 1];
                tempSeq  = diff(seq);
                tempSeq2 = diff(tempSeq);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) | sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            else
                seq = sample_wor([1:5],1,MaxPress);
                tempSeq = diff(seq);
                tempSeq2 = diff(tempSeq);
                while sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            end
            A(e).RandomSeq.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            A(e).RandomSeq.cueC(i,:) ={'£'};
            A(e).RandomSeq.cueS(i,:) ={'£'};
            A(e).RandomSeq.iti(1:i,:) = 500;
            A(e).RandomSeq.hand(i,:) = 2;
            A(e).RandomSeq.sounds(i,:) = 1;
            
            for press= 1:MaxPress
                comnd  = [' A(e).RandomSeq.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        
    end
    RandomSeq = addstruct(A(1).RandomSeq , A(2).RandomSeq);
    for k = 3:length(A)
        RandomSeq = addstruct(RandomSeq , A(k).RandomSeq);
    end
    clear A
    
    idxAll = randperm(length(RandomSeq.FT));
    for kk = 1:length(OrderFields)
        statement = ['RandomSeq.' , OrderFields{kk},' = ','RandomSeq.' , OrderFields{kk},'(idxAll);'];
        eval(statement);
    end
    
    
    name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_RAND' , '_B' , num2str(BN) , '.tgt'];
    dsave(name,orderfields(RandomSeq,OrderFields));
    
    clear x
    clear RandomSeq
end




%% The inrtermixed CLAT and Random blocks


tempCL = rem(ChunkNumb , 100);

CLA = [1 ;2];
CMB.DesiredCnhkargmnt = [4 4 3 3 ; 3 4 3 4];


StimTimeLim = zeros(length(Trials) , 1);
%Generate 20 Full Horizons, and then 4 of each for 3 : 12
Horizon = [1:8 , 13];
clear ChunkArrangeLearn
for BN = 1:20
    Trials = 1:2;
    StimTimeLim = zeros(length(Trials) , 1);
    for e = 1:length(Horizon)
        A(e).ChunkArrangeLearn.StimTimeLim = StimTimeLim;
        X = [];
        for rep = 1: length(Trials) / length(CLA)
            X = [X ; sample_wor([1:length(CLA)],length(CLA))];
        end
        X = X(randperm(length(X))); % so that the the cycle of going through all the possible CLAs in one run is not detected
        A(e).ChunkArrangeLearn.seqNumb = X;
        A(e).ChunkArrangeLearn.FT(1:length(Trials),1) = 2;
        for i = 1:length(X)
            ChunkArrange = CMB.DesiredCnhkargmnt(CLA(X(i)) , :);
            seq = cell(1,length(find(ChunkArrange)) );
            UtempCL = unique(tempCL);
            for j = 1: length(UtempCL)
                numCL = sum(ChunkArrange == UtempCL(j));
                ch_id = find(ChunkArrange == UtempCL(j));
                chnkind{j} = [1 1];
                while diff(chnkind{j}) == 0
                    chnkind{j} = randi(sum(tempCL == UtempCL(j)),numCL,1);
                end
                for n = 1:numCL
                    seq{ch_id(n)} = CMB.Chunks{UtempCL(j)}(chnkind{j}(n) , 1:UtempCL(j));
                end
            end
            seq = cell2mat(seq);
            A(e).ChunkArrangeLearn.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            A(e).ChunkArrangeLearn.cueC(i,:) ={'£'};
            A(e).ChunkArrangeLearn.cueS(i,:) ={'£'};
            A(e).ChunkArrangeLearn.iti(1:i,:) = 500;
            A(e).ChunkArrangeLearn.hand(i,:) = 2;
            A(e).ChunkArrangeLearn.sounds(i,:) = 1;
            A(e).ChunkArrangeLearn.Horizon(i,:) = Horizon(e);
            for press= 1:MaxPress
                comnd  = [' A(e).ChunkArrangeLearn.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
    end
        
        clear x
    
    
        
        
        
        
    for e = 1:length(Horizon)    
        
        ElimChunkStart  = 1; % 1 means that random seqs that start with a known chunk will be elimminated
        
        Trials = 1;
        StimTimeLim = zeros(length(Trials) , 1);
        
        A(e).RandomSeq.StimTimeLim = StimTimeLim;
        A(e).RandomSeq.seqNumb(1:length(Trials),1) = 0;
        A(e).RandomSeq.FT(1:length(Trials),1) = 2;
        A(e).RandomSeq.Horizon(1:length(Trials),:) = Horizon(e);
        for i = 1:length(Trials)
            if ElimChunkStart
                seq = Chunks(1,1:2); % makes sure to enter the first while loop
                seq = [seq  1 1 1 1];
                tempSeq  = diff(seq);
                tempSeq2 = diff(tempSeq);
                while length(unique([seq(1:2) ; Chunks(:,1:2)] , 'rows')) <= length(Chunks) | sum(tempSeq == 0) > 1 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            else
                seq = sample_wor([1:5],1,MaxPress);
                tempSeq = diff(seq);
                tempSeq2 = diff(tempSeq);
                while sum(tempSeq == 0) > 0 | sum(tempSeq == 1) > 1 | sum(tempSeq == -1) > 1
                    seq = sample_wor([1:5],1,MaxPress);
                    tempSeq = diff(seq);
                    tempSeq2 = diff(tempSeq);
                end
            end
            A(e).RandomSeq.cueP{i,:} = char(regexprep(cellstr(num2str(seq)),'\s',''));
            A(e).RandomSeq.cueC(i,:) ={'£'};
            A(e).RandomSeq.cueS(i,:) ={'£'};
            A(e).RandomSeq.iti(1:i,:) = 500;
            A(e).RandomSeq.hand(i,:) = 2;
            A(e).RandomSeq.sounds(i,:) = 1;
            
            for press= 1:MaxPress
                comnd  = [' A(e).RandomSeq.press' , num2str(press) , '(i,1) = seq(press);'];
                eval(comnd);
            end
        end
        
        
    end
    
    ChunkArrangeLearn = addstruct(A(1).ChunkArrangeLearn , A(2).ChunkArrangeLearn);
    for k = 3:length(A)
        ChunkArrangeLearn = addstruct(ChunkArrangeLearn , A(k).ChunkArrangeLearn);
    end
    
    
    RandomSeq = addstruct(A(1).RandomSeq , A(2).RandomSeq);
    for k = 3:length(A)
        RandomSeq = addstruct(RandomSeq , A(k).RandomSeq);
    end
    
    
    Intermixed = addstruct(RandomSeq , ChunkArrangeLearn);

    clear A
    
    idxAll = randperm(length(Intermixed.FT));
    for kk = 1:length(OrderFields)
        statement = ['Intermixed.' , OrderFields{kk},' = ','Intermixed.' , OrderFields{kk},'(idxAll);'];
        eval(statement);
    end
    
    
    name = [Fname ,'/'  , SubCode , num2str(GroupCode) , '_IM' , '_B' , num2str(BN) , '.tgt'];
    dsave(name,orderfields(Intermixed,OrderFields));
    
    clear x
    clear Intermixed
end








