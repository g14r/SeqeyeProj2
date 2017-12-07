function varargout  = se2_linearIPIModels(what , varargin)
% Function to fit and evaluate different linear models for the IPIs
% baseDir = '/Users/nedakordjazi/Documents/SeqEye/SeqEye2/analyze';     %macbook
baseDir = '/Volumes/MotorControl/data/SeqEye2/analyze';  % server
horzSize = {1,2,3,4,5,6,7,8,13,[5:13]};
modelTerms =  {[6]    [1 6] ,   [2 6] ,   [2:3 6]        [2:4 6]          ,[1:2  , 6] ,      [1:3 6]  ,    [1:4 , 6]};
modelNames = {'R'  'C+R',     '1st+R' ,'1st+2nd+R'    '1st+2nd+3rd+R' , 'C+1st+R'  ,  'C+1st+2nd+R'  ,  'Full'};

switch (what)
    case 'loadData'
        type =varargin{1};
        switch(type)
            case 'c'
                load([baseDir , '/se2_TranProb.mat'] , 'C')
                M = C;
                titleSuffix = 'Chunked';
            case 'r'
                load([baseDir , '/se2_TranProb.mat'] , 'R')
                M = R;
                titleSuffix = 'Random';
        end;
        varargout={M,titleSuffix};
    case 'makeDesignMatrix'
        M=varargin{1};  % Input cell array
        L = length(M.IPI);
        
        % I remove intercept here
        X1 = M.IPIarrangement; % within/between chunk
        X1(X1==0) = 1; % Random --> set to between
        X1(X1==2) = -1;  % within
        
        % Map 1st, 2nd 3rd transition probabilities to 5 bins and set them as predictors
        X2 = ceil(1.5*(1 + mapminmax(M.t2Prob_n')))'; % Unclear why this is a good way of splitting (make optional)
        X3 = ceil(1.5*(1 + mapminmax(M.t3Prob_n(:,1)')))';
        X4 = ceil(1.5*(1 + mapminmax(M.t4Prob_n(:,1)')))';
        
        % the learinng regressor - using normalized IPIs, we dont need this
        X5 = M.BN;
        
        % the repetition regressor
        X6 = ismember(M.t2 , [21:25])  +1;
        
        % the design matrix
        M.X = [X1 X2 X3 X4 X5 X6];
        M.Y = M.IPI_norm;
        varargout={M};
    case 'fit'
        [M,titleSuffix] = se2_linearIPIModels('loadData',varargin{1});
        M=se2_linearIPIModels('makeDesignMatrix',M);
        subj = unique(M.SN);
        RR=[]; % Structure to collect the answers
        for sn = subj'   % loop over subjects
            for h = 1:length(horzSize) % loop over horizons
                for dd = 1:5   % loop over days
                    T = getrow(M , ismember(M.SN , sn) & ismember(M.Horizon , horzSize{h}) & ismember(M.Day , dd));
                    R.numObs = length(T.Y);
                    R.SN = sn;
                    R.Horizon = h;
                    R.Day = dd;
                    X = bsxfun(@minus,T.X,mean(T.X)); % Faster and more compact than repmat
                    Y = T.Y-mean(T.Y);                  % Faster and more compact than repmat
                    for m=1:length(modelTerms)
                        b=X(:,modelTerms{m})\Y;                        % Linear regression is small an compact
                        Ypred= X(:,modelTerms{m})*b;                   % So if I don't need any bells and whistels I am coding from scatch
                        res = Y-Ypred;
                        R.numReg(1,m) = length(modelTerms{m});          % I am collecting them here in rows to make comparisions (t-test, normalization) 
                        R.R2(1,m) = 1-sum(res.^2)/sum(Y.^2);            % easier 
                        R.R(1,m)  = corr(Y,Ypred);
                        R.AIC(1,m) = 2*R.numReg(m) + R.numObs*log(sum(res.^2));
                    end;
                    RR=addstruct(RR,R);
                end;
            end;
        end;
        varargout={RR};
    case 'crossval'
        [M,titleSuffix] = se2_linearIPIModels('loadData',varargin{1});
        lambda = varargin{2}; % 0: OLS >0: Ridge
        M=se2_linearIPIModels('makeDesignMatrix',M);
        subj = unique(M.SN);
        RR=[]; % Structure to collect the answers
        for sn = subj'   % loop over subjects
            for h = 1:length(horzSize) % loop over horizons
                for dd = 1:5   % loop over days
                    T = getrow(M , ismember(M.SN , sn) & ismember(M.Horizon , horzSize{h}) & ismember(M.Day , dd));
                    R.numObs = length(T.Y);
                    R.SN = sn;
                    R.Horizon = h;
                    R.Day = dd;
                    CVI = crossvalind('Kfold', R.numObs, CVfol);
                    
                    % Subtract mean from each Fold in X and Y
                    for m=1:length(modelTerms)
                        % Copy here your core crossvalidation code.
                        % Using b=(X'*X+eye(N)*lamba)\Y you can have both Ridge
                        % and OLS regression
                        % I usually collect the predictions over all crossvalidation folds
                        % and then
                        res = Y-Ypred;
                        R.numReg(1,m) = length(modelTerms{m});
                        R.R2(1,m) = 1-sum(res.^2)/sum(Y.^2);
                        R.R(1,m)  = corr(Y,Ypred);
                    end;
                    RR=addstruct(RR,R);
                end
            end
        end
        varargout={RR};
        
        
end;