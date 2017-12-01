function R2_adjusted  = se2_R2Adjusted(Y , Y_pred , numP)


SST = sum((Y - mean(Y)).^2);  % Total Variance


SSE1 = sum((Y-Y_pred).^2); % Residual Variance of the  Model


SSR1 = sum((Y_pred - mean(Y)).^2); % Explained Variance by the model

N = size(Y , 1);
R2_adjusted = 1 - (SSE1/SST)*((N -1 )/(N - numP -1 ));