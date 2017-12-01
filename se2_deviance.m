% Deviation  = -2*ln[(likelihood of fitted model)/(likelihood of saturated model)]    
Mdl = fitglm(Train.X(:,params) , Train.IPI,'Intercept',false);
[Ypred,Posterior] = predict(Mdl,Train.X(:,params));

% sum(log(binopdf(y,n,yfit./n))) - sum(log(binopdf(y,n,y./n)))



pdf_Sat    = pdf('normal',Train.IPI , mean(Train.IPI) , std(Train.IPI));
pdf_fitted = pdf('normal',Train.IPI , mean(Ypred) , std(Ypred));



-2*log(pdf_fitted/pdf_Sat)
dev =  -2*(sum(log(pdf_Sat)) - sum(log(pdf_fitted)));

M = mle(Train.IPI,'distribution','normal')
M = mle(Ypred,'distribution','normal')

-2*log()

-2*sum(log())