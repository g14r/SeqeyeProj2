function stats = se2_SigTest(Dall , what , varargin)

c=  1;
PoolSequences = 0;
while(c<=length(varargin))
    switch(varargin{c})
        case {'seqNumb'}
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'PoolSequences'}
            % whether to pool together all the sequences
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Day'}
            % defines the length of the sequence
            % default is 10
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'Horizon'}
            % defines the number of sequences to be simulated Default = 200
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        case {'whatIPI'}
            % required when what = 'IPI'
            % in what way to test the IPIs
            % 'testBegVsSS' test the start IPIs vs steady state
            % 'testEndVsSS' test the end IPIs vs steady state
            % 'steadyState' test Steady state
            % 'testBegHorz' test if there's an interaction between how many are faster at the beginning and horizon size
            % 'testEndHorz' test if there's an interaction between how many are faster at the end and horizon size
            eval([varargin{c} '= varargin{c+1};']);
            c=c+2;
        otherwise
            error(sprintf('Unknown option: %s',varargin{c}));
    end
end

ANA = getrow(Dall , Dall.isgood & ~Dall.isError & ...
    ismember(Dall.Horizon , Horizon) & ...
    ismember(Dall.Day , Day) & ismember(Dall.seqNumb , seqNumb));
ANA.RT = ANA.AllPressTimes(:,1);
ANA.seqNumb(ANA.seqNumb>1) = 1;
if PoolSequences
    ANA.seqNumb = zeros(size(ANA.seqNumb));
end
ANA.Day(ANA.Day == 3) = 2;
ANA.Day(ismember(ANA.Day , [4,5])) = 3;

Day(Day==3) = 2;
Day(ismember(Day , [4 5])) = 3;
Day = unique(Day);

factors = {'Horizon' , 'Day' , 'seqNumb'};
facInclude = [length(Horizon)>1 , length(Day)>1  , length(unique(ANA.seqNumb))>1];
FCTR =  factors(facInclude);

switch what
    case 'MT'
        var = [];
        for f = 1:length(FCTR)
            eval(['var = [var ANA.',FCTR{f},'];']);
        end
        stats = anovaMixed(ANA.MT  , ANA.SN ,'within',var ,FCTR,'intercept',1) ;
%         anovan(ANA.MT,var,'model','interaction','varnames',FCTR)  % between subject
        figure('color' , 'white')
        lineplot(var, ANA.MT , 'style_thickline');
        title(['Effect of ' , FCTR , ' on ' , what]);
    case 'IPI'
        switch whatIPI
            case 'steadyState'
                ipi = ANA.IPI(:,4:10);
                %         ipi = ANA.IPI(:,1:3);
                A.IPI = reshape(ipi , numel(ipi) , 1);
                A.Horizon = repmat(ANA.Horizon , size(ipi , 2) , 1);
                A.Day = repmat(ANA.Day , size(ipi , 2) , 1);
                A.seqNumb = repmat(ANA.seqNumb , size(ipi , 2) , 1);
                A.SN = repmat(ANA.SN , size(ipi , 2) , 1);
                var = [];
                for f = 1:length(FCTR)
                    eval(['var = [var A.',FCTR{f},'];']);
                end
                stats = anovaMixed(A.IPI  , A.SN ,'within',var ,FCTR,'intercept',1) ;
                anovan(A.IPI,var,'model','interaction','varnames',FCTR)  % between subject
                figure('color' , 'white')
                lineplot(var, A.IPI, 'style_thickline');
                title(['Effect of ' , FCTR , ' on ' , what]);
            case 'testBegVsSS'
                ipiss = ANA.IPI(:,4:10);
                nn= input('How many IPIs from the beggining to consider as start IPIs??');
                ipibeg = ANA.IPI(:,nn);
                
                label = [ones(numel(ipibeg) , 1) ;zeros(numel(ipiss) ,1)];
                A.IPI = [reshape(ipibeg , numel(ipibeg) , 1);reshape(ipiss , numel(ipiss) , 1)];
                A.Horizon = repmat(ANA.Horizon , size(ipibeg , 2)+size(ipiss , 2) , 1);
                A.Day = repmat(ANA.Day , size(ipibeg , 2)+size(ipiss , 2) , 1);
                A.seqNumb = repmat(ANA.seqNumb , size(ipibeg , 2)+size(ipiss , 2) , 1);
                A.SN = repmat(ANA.SN , size(ipibeg , 2)+size(ipiss , 2) , 1);
                var = [];
                for f = 1:length(FCTR)
                    eval(['var = [var A.',FCTR{f},'];']);
                end
                var = [var label];
                FCTR = [FCTR , 'Beggining/SS'];
                stats = anovaMixed(A.IPI  , A.SN ,'within',var ,FCTR,'intercept',1) ;
                anovan(A.IPI,var,'model','interaction','varnames',FCTR)  % between subject
                figure('color' , 'white')
                lineplot(var, A.IPI, 'style_thickline');
                title(['Effect of ' , FCTR , ' on ' , what]);
            case 'testEndVsSS'
                ipiss = ANA.IPI(:,4:10);
                nn= input('How many IPIs from the eng to consider as end IPIs??');
                ipiend = ANA.IPI(:,nn);
                label = [ones(numel(ipiend) , 1) ;zeros(numel(ipiss) ,1)];
                A.IPI = [reshape(ipiend , numel(ipiend) , 1);reshape(ipiss , numel(ipiss) , 1)];
                A.Horizon = repmat(ANA.Horizon , size(ipiend , 2)+size(ipiss , 2) , 1);
                A.Day = repmat(ANA.Day , size(ipiend , 2)+size(ipiss , 2) , 1);
                A.seqNumb = repmat(ANA.seqNumb , size(ipiend , 2)+size(ipiss , 2) , 1);
                A.SN = repmat(ANA.SN , size(ipiend , 2)+size(ipiss , 2) , 1);
                var = [];
                for f = 1:length(FCTR)
                    eval(['var = [var A.',FCTR{f},'];']);
                end
                var = [var label];
                FCTR = [FCTR , 'Ending/SS'];
                stats = anovaMixed(A.IPI  , A.SN ,'within',var ,FCTR,'intercept',1) ;
                anovan(A.IPI,var,'model','interaction','varnames',FCTR)  % between subject
                figure('color' , 'white')
                lineplot(var, A.IPI, 'style_thickline');
                title(['Effect of ' , FCTR , ' on ' , what]);
            case 'WithinBetweenChunk'
                ipiss = ANA.IPI(:,4:10);
                FCTR = FCTR(~strcmp(FCTR,'seqNumb'));
                nn= input('What IPIs to include?(0 random, 1 between, 2 within)');
                ipiLab = {'Random' , 'Between','Within'};
                ipiLab = ipiLab(nn+1);
                L = ipiLab{1};
                for l = 2:length(ipiLab)
                    L = [L,'/',ipiLab{l}];
                end
                A.IPIArr= reshape(ANA.IPIarrangement(:,4:10) , numel(ipiss) , 1);
                A.IPI = reshape(ipiss , numel(ipiss) , 1);
                A.Horizon = repmat(ANA.Horizon , size(ipiss , 2) , 1);
                A.Day = repmat(ANA.Day , size(ipiss , 2) , 1);
                A.seqNumb = repmat(ANA.seqNumb , size(ipiss , 2) , 1);
                A.SN = repmat(ANA.SN , size(ipiss , 2) , 1);
                A = getrow(A , ismember(A.IPIArr , nn));
                var = [];
                for f = 1:length(FCTR)
                    eval(['var = [var A.',FCTR{f},'];']);
                end
                var = [var A.IPIArr];
                FCTR = [FCTR , L];
                stats = anovaMixed(A.IPI  , A.SN ,'within',var ,FCTR,'intercept',1) ;
                anovan(A.IPI,var,'model','interaction','varnames',FCTR)  % between subject
                figure('color' , 'white')
                lineplot(var, A.IPI, 'style_thickline');
                title(['Effect of ' , FCTR , ' on ' , what]);
        end
    case 'RT'
        var = [];
        for f = 1:length(FCTR)
            eval(['var = [var ANA.',FCTR{f},'];']);
        end
        stats = anovaMixed(ANA.RT  , ANA.SN ,'within',var ,FCTR,'intercept',1) ;
        figure('color' , 'white')
        lineplot(var, ANA.RT-1500, 'style_thickline');
        title(['Effect of ' , FCTR , ' on ' , what]);
end