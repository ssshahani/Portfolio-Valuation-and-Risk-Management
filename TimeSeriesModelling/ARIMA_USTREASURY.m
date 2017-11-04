
data = xlsread("USTREASURY-YIELD.xlsx",'D:D');

% Plotting original dataset
figure
title("USTREASURY-YIELD")
subplot(1,2,1)
plot(data)
title("D=0")
subplot(1,2,2)
autocorr(data)

% Check for stationarity 
adftest(data) % If logical 1 returned - data is stationary, otherewise non-stationary.
pptest(data);
kpsstest(data);
vratiotest(data);

% To make data stationary, uncomment this code
% if(adftest(data) ~= 1)
%     diffDataOnce = zeros(246,1); 
%     for i =1:length(data) 
%        if (i == length(data))
%            break;
%        end
%        diffDataOnce(i) = (data(i+1) - data(i));
%     end
% end
% splitting data into training set and testing set
%

% Splitting Data Set into training and testing data set.
training_data = data(1:181);
testing_data = data(182:184);
 
% Use AIC and BIC criterion to determine ARMA lags
LOGL = zeros(3,3); %Initialize
PQ = zeros(3,3);
aicManual = zeros(3,3);
bicManual =zeros(3,3);
for p = 1:3
    for q = 1:3
        mod = arima(p-1,0 ,q-1);
        [fit,~,logL] = estimate(mod,training_data);
        LOGL(p,q) = logL;
        PQ(p,q) = p+q - 2;
        aicManual(p,q) =  -2 * LOGL(p,q) + 2 * PQ(p,q);
        bicManual(p,q) =  -2 * LOGL(p,q) + log(length(training_data)) * PQ(p,q);
     end
end
% 
% Manual Matrix Print
aicManual
bicManual

LOGL = reshape(LOGL,9,1);
PQ = reshape(PQ,9,1);

[aic,bic] = aicbic(LOGL,PQ+1,length(training_data));
title("AIC MATRIX")
reshape(aic,3,3)
title("BIC MATRIX")
reshape(bic,3,3)
% 
% Best Model after analyzing the result is ARIMA (2,0,2). 
mod = arima(1, 0 ,1);
[fit,~,logL,residualsModel] = estimate(mod,training_data);
LOGL = logL;

% Plot the residuals
E = infer(fit,training_data(3:end),'Y0',training_data(1:2));
figure;
plot(E);
title 'Inferred Residuals';
% 
% % For predicting the response to testing data, we need the last residual of
% % the model fitted to training data.
residual2 = E(178);
residual1 = E(179);
% 
% Extracting the coefficients for ARIMA (2,0,2)
const = fit.Constant;
phi = fit.AR;
theta =fit.MA;
% muhat = const/(1 - (phi{1})) 
yt = training_data(end);
%  yprev = training_data((end-1));
% 
% Our model ARIMA(1,0,1) - Predicting next three data points.
yt1 = (const ) + phi{1}*(yt - const)  - theta{1}*residual1  
yt2 = (const ) + phi{1}*(yt1 - const)
yt3 = (const ) + phi{1}*(yt2 - const)
%  yt3 = (muhat ) + phi{1}*(yt2 - muhat) + phi{2}*(yt1 -muhat)

% Forecasting next 50 data points
% ytArray = zeros(50,1);
% for i = 1:length(ytArray)
%     if(i == 1)
%         ytArray(i) = (muhat ) + phi{1}*(yt3 - muhat) + phi{2}*(yt2 - muhat);
%     elseif(i == 2)
%         ytArray(i) = (muhat ) + phi{1}*(ytArray(i-1) - muhat) + phi{2}*(yt3 - muhat);
%     else
%         ytArray(i) = (muhat ) + phi{1}*(ytArray(i-1) - muhat) + phi{2}*(ytArray(i-2) - muhat);
%     end
% end

% Forecast Method
% ytArray = forecast(fit,50,'Y0', training_data);
