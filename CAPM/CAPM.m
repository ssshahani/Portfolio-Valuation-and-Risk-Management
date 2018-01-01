filename = "SPY-ASSETS.xlsx";
[covar_Asset_Market_array, var_Market] = calcBeta(filename);
beta = ones(length(covar_Asset_Market_array),1);

sum_covar = 0;
sum_var = 0;
for i=1:length(covar_Asset_Market_array)  
    beta(i,1) = covar_Asset_Market_array(i)/var_Market;
end

figure 
hist(beta)

Stocks = string({'AAPL';'ABBV';'COST';'CSCO';'CVX';'GE';'IBM';'JNJ';'JPM';'KO';'MCD';'MSFT';'PEP';'PFE';'PG';'SBUX';'T';'VZ';'WFC';'XOM'});
betaEstimate(filename);
traditional_beta = beta;

LeastTrimmedSquares_Beta = leastTrimmedSquares(filename);
shrinkage_beta = betaShrinkageEstimation(beta);
[Beta_EWMA_lamda1] = calcBetaUsingEWMA(filename, 0.94);
[Beta_EWMA_lamda2] = calcBetaUsingEWMA(filename, 0.97);

sorted_beta = sort(traditional_beta);
Stocks_Sorted = string({'KO';'T';'VZ';'PEP';'PG';'COST';'JNJ';'MCD';'PFE';'XOM';'IBM';'GE';'SBUX';'ABBV';'CSCO';'AAPL';'WFC';'CVX';'MSFT';'JPM'});
riskiness = strings(0);
for i=1:length(sorted_beta)
    if (sorted_beta(i) > 1)
        riskiness = [riskiness 'Aggressive'];
    elseif(sorted_beta(i) == 1)
        riskiness = [riskiness 'AverageRisk'];
    elseif (sorted_beta(i) < 1)
        riskiness = [riskiness 'NonAggressive'];
    end
end

risk_reward = transpose(riskiness);
t_traditional = table(Stocks_Sorted,sorted_beta,risk_reward);
t_LTS = table(Stocks,traditional_beta,LeastTrimmedSquares_Beta,shrinkage_beta,Beta_EWMA_lamda1, Beta_EWMA_lamda2); 

function [beta_LTS] = leastTrimmedSquares(filename)
    alpha = 0.05
    data = xlsread(filename);
    [rows,col] = size(data);
    marketReturns = assetLogReturns(data(:,1));
    constant = ones(rows-1,1);
    X = [constant marketReturns];
    t_X =transpose(X);
    beta_LTS = ones(col-1,1);
    for j = 2:col
        assetReturns = assetLogReturns(data(:,j));
        y_Ret = assetReturns;
        sorted_yRet =  sort(y_Ret);
        
        quantile_alpha = round(alpha*length(sorted_yRet));
        quantile_1minusAlpha = round((1-alpha)*length(sorted_yRet));

        quantile_alpha_value = sorted_yRet(quantile_alpha);
        quantile_1minusAlpha_value = sorted_yRet(quantile_1minusAlpha);
        if(j == 2)
            figure
            hold on
            xlabel('Actual Asset Returns for AAPL');
            ylabel('Estimated Asset Returns for AAPL');
            title('Before Winsorization for AAPL');
            hist(y_Ret);
            hold off
        end
        for i =1:length(y_Ret)
            if (y_Ret(i) < quantile_alpha_value)
                y_Ret(i) = quantile_alpha_value;
            elseif(y_Ret(i) > quantile_1minusAlpha_value)
                y_Ret(i) = quantile_1minusAlpha_value;
            else
                y_Ret(i) = y_Ret(i);
            end
        end
        if (j==2)
            figure
            hold on
            xlabel('Actual Asset Returns for AAPL')
            ylabel('Estimated Asset Returns for AAPL');
            title('After Winsorization for AAPL');
            hist(y_Ret);
            hold off
        end
        betaVec = inv((t_X * X))*t_X*y_Ret;
        beta_LTS(j-1) =  betaVec(2,1);
        winsorized_yRet_hat = X * betaVec;
    end
    % mdl = fitlm(marketReturns, representativeAssetReturns);
    % representativeAssetReturns_hat = mdl.Fitted;

    % residuals = winsorized_yRet - winsorized_yRet_hat;

    % define loss function

    % row_eps = ones(length(residuals),1);
    % for i = 1:length(residuals)
    %     if(residuals(i) < quantile_alpha_value)
    %         row_eps(i,1) = quantile_alpha_value^2;
    %     elseif (and (residuals(i) > quantile_alpha_value, residuals(i) < quantile_1minusAlpha_value))
    %         row_eps(i,1) = residuals(i)^2; 
    %     else
    %         row_eps(i,1) = quantile_1minusAlpha_value^2;
    %     end
    % end
    % 
    % temp = row_eps*residuals;
    % 
    % 
    % summation = sum(temp);
end

function betaEstimate(filename)
    data = xlsread(filename);
    [rows,col] = size(data);
    spy_etf_data = data(:,1);
    [market_logReturns] = assetLogReturns(spy_etf_data);
    mean_market = mean(market_logReturns);
    sum_var_market = 0;
    beta_vec = zeros(length(col-1));
    for i =1:length(market_logReturns)
        sum_var_market = sum_var_market + (market_logReturns(i)- mean_market)^2;
    end
    for j=2:col
        adj_close_price = data(:,j);
        [asset_logReturns] = assetLogReturns(adj_close_price);
        mean_asset = mean(asset_logReturns);
        sum_asset_market = 0;
        for i =1:length(market_logReturns)
            sum_asset_market = sum_asset_market + (asset_logReturns(i) - mean_asset)*(market_logReturns(i) -mean_market);
        end
        beta_vec(j-1) = sum_asset_market/sum_var_market;
    end
    portfolio_beta = mean(beta_vec);
    display(portfolio_beta);
end


function [covar_Asset_Market_array, var_Market] = calcBeta(filename)
    data = xlsread(filename);
    [rows,col] = size(data);
    daily_risk_free_rate = ones(rows-1,1)*0.003;
    logReturnsMatrix = ones(rows-1,2);
    covar_Asset_Market_array = ones(col-1,1);
        for j=1:col
            adj_close_price = data(:,j);
            asset_logReturns = assetLogReturns(adj_close_price);
            asset_logReturns = asset_logReturns - daily_risk_free_rate;
            if j==1                  
                 logReturnsMatrix(:,1)=asset_logReturns;
                 var_Market = (std(asset_logReturns)^2);
            else
                logReturnsMatrix(:,2)=asset_logReturns;
                covar_Asset_Market = cov(logReturnsMatrix);
                covar_Asset_Market_array(j-1,1) = covar_Asset_Market(1,2);
            end
        end
portfoliobeta();        
end

function [asset_log_returns] = assetLogReturns(adj_close_price)
    logReturns = ones(length(adj_close_price)-1,1);
    for i = 1:(length(adj_close_price)-1)
        returns = (adj_close_price(i+1)/adj_close_price(i)) - 1;
        logReturns(i,1) = log(1+returns) * 100;
    end
    asset_log_returns = logReturns;
end

function [beta_JS] = betaShrinkageEstimation(ols_Beta)
    mean_olsBeta = mean(ols_Beta);
    const_alpha = 2/3;
    % James-Stein Estimator
    beta_JS = ones(length(ols_Beta),1);
    for i=1:length(ols_Beta)
        beta_JS(i) = mean_olsBeta + const_alpha*(ols_Beta(i) - mean_olsBeta);
    end
end

function [EWMA] = calcBetaUsingEWMA(filename, lamda)
    data = xlsread(filename);
    [rows, col] = size(data);
    marketReturns = assetLogReturns(data(:,1));
    variance_market_returns_t = zeros(size(marketReturns));
    variance_market_returns_t(1,1) = (1-lamda)*marketReturns(1,1)^2;
    EWMA = ones(col-1,1);
    for j=2:col
        assetReturns = assetLogReturns(data(:,j));
        covar_Asset_Market_time_t = zeros(size(assetReturns));
        covar_Asset_Market_time_t(1,1) = 0;
        for t =2:length(assetReturns)
            covar_Asset_Market_time_t(t,1) = (1-lamda) * assetReturns(t-1) * marketReturns(t-1) + lamda * covar_Asset_Market_time_t(t-1,1);                   
            variance_market_returns_t(t,1) = (1-lamda) * marketReturns(t,1)^2 + lamda*variance_market_returns_t(t-1,1);
        end
        mean_covar_daily =  mean(covar_Asset_Market_time_t);
        mean_var_daily = mean(variance_market_returns_t);
        EWMA(j-1) = mean_covar_daily/mean_var_daily;
    end
end