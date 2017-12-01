function [R2]  = se2_R2ModelComp(Y , Y_predNull , YpredNew)

SST = sum((Y - mean(Y)).^2);  % Total Variance


SSE0 = sum((Y - Y_predNull).^2); % Residual Variance of the Null Model
SSE1 = sum((Y-YpredNew).^2); % Residual Variance of the  Model


SSR0 = sum((Y_predNull - mean(Y)).^2); % Explained Variance by Null the model
SSR1 = sum((YpredNew - mean(Y)).^2); % Explained Variance by the model




R2   = 1 - (SSE1/SSE0);


