function [Dev , lh_comp] = se2_crossval_DEVandLH(Yte , Ypred , Ypred0)

% Yte    = testing set
% Ypred  = fitted model's predition of the Test set
% Ypred0 = null model's predition of the Test set

% Deviation  = -2*ln[(likelihood of fitted model)/(likelihood of saturated model)]    
%%
SSR0 = sum((Ypred0 - mean(Yte)).^2); % Explained Variance by Null the model
SSR1 = sum((Ypred - mean(Yte)).^2); % Explained Variance by the model

SSE0 = sum((Yte - Ypred0).^2); % Residual Variance of the Null Model
SSE1 = sum((Yte-Ypred).^2); % Residual Variance of the  Model

ll1  =  log(exp(-0.5*SSE1));
ll0  =  log(exp(-0.5*SSE0));


pdf_Sat    = pdf('normal',Yte , mean(Yte) , std(Yte));
pdf_fitted = pdf('normal',Yte  , mean(Ypred) , std(Ypred));
pdf_null   = pdf('normal',Yte , mean(Ypred0) , std(Ypred0));



lh_comp =  max(-2*(sum(log(pdf_fitted)) - sum(log(pdf_null))) , 0);

Dev =  max(-2*(sum(log(pdf_Sat)) - sum(log(pdf_fitted))) , 0);

