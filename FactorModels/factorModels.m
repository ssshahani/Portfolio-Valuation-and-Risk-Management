capm_filename = "CapmData.xlsx";
french_fama_filename = "new_data_FF.xlsx";
capm_F1_filename = "Capm_Factor1Data.xlsx"
capm_F2_filename = "Capm_Factor2Data.xlsx";
capm_F1_F2_filename = "Capm_Factor1_Factor2Data.xlsx";
french_fama_F1_filename = "FrenchFamaF1_Data.xlsx";
french_fama_F2_filename = "FrenchFama_F2Data.xlsx";
french_fama_F1_F2Data_filename = "FrenchFama_F1_F2Data.xlsx";

[CAPMVR, CAPMER, CAPMCov, rSquarecapmModel] = capmModel(capm_filename);
CAPMVR = CAPMVR .* 210;
[FFVR, FFER, FFCov, rSquarefrenchFamaModel]  = frenchFamaModel(french_fama_filename);
FFVR = FFVR .* 210;

data = xlsread(capm_filename);
[rows,col] = size(data);
asset_returns = data(:,2:end);
historical_returns = ones(col-1,1);
historical_variance_estimates = ones(col-1,1);
annulized_expected_returns = ones(col-1,1);

for i=1:col-1
    temp =  mean(asset_returns(:,i));
    historical_returns(i) = temp* 210;
    historical_variance_estimates(i) = var(asset_returns(:,i)) * 210;
end

historical_returns_covariance = cov(asset_returns);

[CAPMF1VR, CAPMF1ER, CAPM49IndustryPortfolios, rSquarecapmF1Model] = capmF1Model(capm_F1_filename);
[CAPMF2VR, CAPMF2ER, CAPMTermSpread, rSquarecapmF2Model] = capmF2Model(capm_F2_filename);
[CAPMF1F2VR, CAPMF1F2ER, CAPM49IndustryPortfoliosTermSpread, rSquarecapmF1F2Model] = capmF1F2Model(capm_F1_F2_filename);

CAPMF1VR = CAPMF1VR .* 210;
CAPMF2VR =  CAPMF2VR .* 210;
CAPMF1F2VR = CAPMF1F2VR .* 210;

[FFF1VR, FFF1ER, FF49IndustryPortfolios,rSquarefrenchFama_F2Model] = frenchFama_F2Model(french_fama_F2_filename);
[FFF2VR, FFF2ER, FFTermSpread, rSquarefrenchFama_F1_F2Model] = frenchFama_F1_F2Model(french_fama_F1_F2Data_filename);
[FFF1F2VR, FFF1F2ER, FF49IndustryPortfoliosTermSpread, rSquarefrenchFama_F1Model] = frenchFama_F1Model(french_fama_F1_filename);

FFF1VR = FFF1VR .* 210;
FFF2VR = FFF2VR .* 210;
FFF1F2VR = FFF1F2VR .* 210;

models = string({'CAPM';'CAPM_F1'; 'CAPM_F2'; 'CAPM_F1_F2'; 'FF';'FF_F1'; 'FF_F2'; 'FF_F1_F2' });
values = {rSquarecapmModel;rSquarecapmF1Model;rSquarecapmF2Model;rSquarecapmF1F2Model; rSquarefrenchFamaModel;rSquarefrenchFama_F1Model;rSquarefrenchFama_F2Model; rSquarefrenchFama_F1_F2Model};
t_RSquare = table(models, values);

Stocks = string({'AAPL';'ABBV';'COST';'CSCO';'CVX';'GE';'IBM';'JNJ';'JPM';'KO';'MCD';'MSFT';'PEP';'PFE';'PG';'SBUX';'T';'VZ';'WFC';'XOM'});
t_ExpectedReturnsModels = table(Stocks,CAPMER, CAPMF1ER,CAPMF2ER, CAPMF1F2ER,FFER, FFF1ER, FFF2ER, FFF1F2ER, historical_returns);
t_VarianceReturnsModels = table(Stocks, CAPMVR  , CAPMF1VR , CAPMF2VR , CAPMF1F2VR ,FFVR , FFF1VR , FFF2VR , FFF1F2VR , historical_variance_estimates );
cov_list= {CAPMCov, CAPM49IndustryPortfolios, CAPMTermSpread, CAPM49IndustryPortfoliosTermSpread,FFCov, FF49IndustryPortfolios, FFTermSpread, FF49IndustryPortfoliosTermSpread, historical_returns_covariance}; 
[norm_mat] = compare_cov(cov_list);
row_names = {'CAPM', 'CAPM_F1', 'CAPM_F2', 'CAPM_F1_F2','FF', 'FF_F1', 'FF_F2', 'FF_F1_F2', 'Historical_Returns'};
col_names = {'CAPM', 'CAPM_F1', 'CAPM_F2', 'CAPM_F1_F2','FF', 'FF_F1', 'FF_F2', 'FF_F1_F2', 'Historical_Returns'};
sTable = array2table(norm_mat,'RowNames',row_names,'VariableNames',col_names)


selected_model_1_FFF1 = tangencyPortfolioWeights(FF49IndustryPortfolios,FFF1ER);
selected_model_2_FFF1F2 = tangencyPortfolioWeights(FF49IndustryPortfoliosTermSpread,FFF1F2ER);
historical_weights = tangencyPortfolioWeights(historical_returns_covariance, historical_expected_returns);

row_names = {'AAPL';'ABBV';'COST';'CSCO';'CVX';'GE';'IBM';'JNJ';'JPM';'KO';'MCD';'MSFT';'PEP';'PFE';'PG';'SBUX';'T';'VZ';'WFC';'XOM'};
t_OptimalWeights = table(Stocks, selected_model_1_FFF1, selected_model_2_FFF1F2, historical_weights);


function [CAPMF1VR, CAPMF1ER, CAPM49IndustryPortfolios, rSquarecapmF1Model] = capmF1Model(filename)
 data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    market_returns = data(:,2) - daily_risk_free_rate_vector;
    factor_data = data(:,1);
    Sigma_F = cov(factor_data,market_returns);
    X_Factor = [constant factor_data market_returns];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(2 ,col-2);
    alpha_vector = ones(col-2,1);
    expected_returns_asset = ones(col-2,1);
    variance_returns_asset = ones(col-2,1);
    ten_months_expected_returns = ones(col-2,1);
    rsquared = ones(col-2,1);
    for i=3:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-2) = beta_vec(1);
        Beta_Matrix(:,i-2) = beta_vec(2:end);
        expected_returns_on_asset = (daily_risk_free_rate + alpha_vector(i-2,1)) + (Beta_Matrix(1,i-2) * mean(factor_data)) + (Beta_Matrix(2,i-2) * mean(market_returns));
        ten_months_expected_returns(i-2) = expected_returns_on_asset*210;  
        estimated_returns_of_asset = daily_risk_free_rate_vector + alpha_vector(i-2) + (Beta_Matrix(1,i-2) * factor_data) + (Beta_Matrix(2,i-2) * market_returns);
        rsquared(i-2) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-2) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-2)' * Sigma_F * Beta_Matrix(:,i-2)) + residual'*residual;
        variance_returns_asset(i-2) = var_asset_return;      
    end
    
    covariance = ones(col-2,col-2);
    for i=1:col-2
        for j = (i+1):col-2
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end

    x = triu(covariance);
    covariance = x + x';
    for i=1:col-2
        covariance(i,i) =  variance_returns_asset(i);
    end
    CAPM49IndustryPortfolios = covariance;
    CAPMF1ER = ten_months_expected_returns;
    CAPMF1VR =  variance_returns_asset;
    rSquarecapmF1Model = mean(rsquared);

end

function [CAPMF2VR, CAPMF2ER, CAPMTermSpread, rSquarecapmF2Model] = capmF2Model(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    market_returns = data(:,2) - daily_risk_free_rate_vector;
    factor_data = data(:,1);
    Sigma_F = cov(factor_data,market_returns);
    X_Factor = [constant factor_data market_returns];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(2 ,col-2);
    alpha_vector = ones(col-2,1);
    expected_returns_asset = ones(col-2,1);
    variance_returns_asset = ones(col-2,1);
    ten_months_expected_returns = ones(col-2,1);
    rsquared = ones(col-2,1);
    for i=3:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-2) = beta_vec(1);
        Beta_Matrix(:,i-2) = beta_vec(2:end);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-2,1) + Beta_Matrix(1,i-2) * mean(factor_data) + Beta_Matrix(2,i-2) * mean(market_returns);
        ten_months_expected_returns(i-2) = expected_returns_on_asset*210;  
        estimated_returns_of_asset = daily_risk_free_rate_vector + alpha_vector(i-2) + Beta_Matrix(1,i-2) * factor_data + Beta_Matrix(2,i-2) * market_returns;
        rsquared(i-2) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-2) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-2)' * Sigma_F * Beta_Matrix(:,i-2)) + residual'*residual;
        variance_returns_asset(i-2) = var_asset_return;      
    end
    
    covariance = ones(col-2,col-2);
    for i=1:col-2
        for j = (i+1):col-2
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-2
        covariance(i,i) =  variance_returns_asset(i);
    end
    CAPMTermSpread = covariance;
    CAPMF2ER = ten_months_expected_returns; 
    CAPMF2VR = variance_returns_asset;
    rSquarecapmF2Model = mean(rsquared);
end

function [CAPMF1F2VR, CAPMF1F2ER, CAPM49IndustryPortfoliosTermSpread, rSquarecapmF1F2Model] =capmF1F2Model(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    market_returns = data(:,3) - daily_risk_free_rate;
    factor1_data = data(:,1);
    factor2_data = data(:,2);
    X_Factor = [factor1_data factor2_data market_returns];
    Sigma_F = cov(X_Factor);
    X_Factor = [constant factor1_data factor2_data market_returns];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(3 ,col-3);
    alpha_vector = ones(col-3,1);
    expected_returns_asset = ones(col-3,1);
    variance_returns_asset = ones(col-3,1);
    ten_months_expected_returns = ones(col-3,1);
    rsquared = ones(col-3,1);
    for i=4:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-3) = beta_vec(1);
        Beta_Matrix(:,i-3) = beta_vec(2:end);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-3) + Beta_Matrix(1,i-3) * mean(factor1_data) + Beta_Matrix(2,i-3) * mean(factor2_data) + Beta_Matrix(3,i-3) * mean(market_returns);
        ten_months_expected_returns(i-3) = expected_returns_on_asset * 210;  
        estimated_returns_of_asset = daily_risk_free_rate_vector + alpha_vector(i-3) + Beta_Matrix(1,i-3) * factor1_data + Beta_Matrix(2,i-3) * factor2_data + Beta_Matrix(3,i-3) * market_returns;
        rsquared(i-3) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-3) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-3)' * Sigma_F * Beta_Matrix(:,i-3)) + residual'*residual;
        variance_returns_asset(i-3) = var_asset_return;      
    end
    
    covariance = ones(col-3,col-3);
    for i=1:col-3
        for j = (i+1):col-3
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-3
        covariance(i,i) =  variance_returns_asset(i);
    end
    CAPM49IndustryPortfoliosTermSpread = covariance;
    CAPMF1F2ER = ten_months_expected_returns;
    CAPMF1F2VR =  variance_returns_asset;
    rSquarecapmF1F2Model = mean(rsquared);
end

function [CAPMVR, CAPMER, CAPMCov, rSquarecapmModel] =capmModel(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    market_returns = data(:,1) - daily_risk_free_rate_vector;
    Sigma_F = cov(market_returns);
    X_Factor = [constant  market_returns];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(1 ,col-1);
    alpha_vector = ones(col-1,1);
    expected_returns_asset = ones(col-1,1);
    variance_returns_asset = ones(col-1,1);
    ten_months_expected_returns = ones(col-1,1);
    rsquared = ones(col-1,1);
    for i=2:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-1) = beta_vec(1);
        Beta_Matrix(:,i-1) = beta_vec(2);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-1,1) + Beta_Matrix(1,i-1) * mean(market_returns);
        ten_months_expected_returns(i-1) = expected_returns_on_asset*210;  
        estimated_returns_of_asset = daily_risk_free_rate_vector + alpha_vector(i-1) + Beta_Matrix(1,i-1) * market_returns;
        rsquared(i-1) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-1) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-1)' * Sigma_F * Beta_Matrix(:,i-1)) + residual'*residual;
        variance_returns_asset(i-1) = var_asset_return;      
    end
    
    covariance = ones(col-1,col-1);
    for i=1:col-1
        for j = (i+1):col-1
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-1
        covariance(i,i) =  variance_returns_asset(i);
    end
    CAPMCov = covariance;
    CAPMER = ten_months_expected_returns; 
    CAPMVR = variance_returns_asset;
    rSquarecapmModel = mean(rsquared);
end

function [Rsquared] = calculateRSquared(y, y_hat)
    [sum_RSS, total_variance] = anovaTotal(y,y_hat);
    Rsquared = sum_RSS/total_variance;
end

function [sum_RSS, total_variance]  = anovaTotal(y, y_hat)
    mean_y = mean(y);
    sum_RSS = 0;
    for i =1:size(y)
        sum_RSS = sum_RSS + (y_hat(i) - mean_y)^2;
    end
    
    sum_SSE = 0;
    for i =1:size(y)
        sum_SSE = sum_SSE + (y_hat(i) - y(i))^2;
    end
    total_variance = sum_RSS + sum_SSE;
end

function [FFVR, FFER, FFCov, rSquareFrenchFamaModel] =frenchFamaModel(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    excess_market_returns = data(:,1) - daily_risk_free_rate_vector;
    small_big_returns = data(:,2);
    high_low_returns = data(:,3);
    X_Factor = [excess_market_returns small_big_returns high_low_returns];
    Sigma_F = cov(X_Factor);
    X_Factor = [constant X_Factor];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(3 ,col-3);
    alpha_vector = ones(col-3,1);
    expected_returns_asset = ones(col-3,1);
    variance_returns_asset = ones(col-3,1);
    ten_months_expected_returns = ones(col-3,1);
    rsquared = ones(col-3,1);
    for i=4:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-3) = beta_vec(1);
        Beta_Matrix(:,i-3) = beta_vec(2:end);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-3,1) + Beta_Matrix(1,i-3) * mean(excess_market_returns) + Beta_Matrix(2,i-3) * mean(small_big_returns) + Beta_Matrix(3,i-3) * mean(high_low_returns) ;
        ten_months_expected_returns(i-3) = expected_returns_on_asset*210;  
        estimated_returns_of_asset = daily_risk_free_rate_vector + alpha_vector(i-3) + Beta_Matrix(1,i-3) * excess_market_returns + Beta_Matrix(2,i-3) * small_big_returns + Beta_Matrix(3,i-3) * high_low_returns;
        rsquared(i-3) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-3) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-3)' * Sigma_F * Beta_Matrix(:,i-3)) + residual'*residual;
        variance_returns_asset(i-3) = var_asset_return;      
    end
    
    covariance = ones(col-3,col-3);
    for i=1:col-3
        for j = (i+1):col-3
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-3
        covariance(i,i) =  variance_returns_asset(i);
    end
    FFCov = covariance;
    FFER = ten_months_expected_returns; 
    FFVR = variance_returns_asset;
    rSquareFrenchFamaModel = mean(rsquared);
end

function [FFF1VR, FFF1ER, FF49IndustryPortfolios, rSquarefrenchFama_F1Model] =frenchFama_F1Model(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    excess_market_returns = data(:,1) - daily_risk_free_rate_vector;
    small_big_returns = data(:,2);
    high_low_returns = data(:,3);
    factor1_data = data(:,4);
    X_Factor = [excess_market_returns small_big_returns high_low_returns factor1_data];
    Sigma_F = cov(X_Factor);
    X_Factor = [constant X_Factor];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(4 ,col-4);
    alpha_vector = ones(col-4,1);
    expected_returns_asset = ones(col-4,1);
    variance_returns_asset = ones(col-4,1);
    ten_months_expected_returns = ones(col-4,1);
    rsquared = ones(col-4,1);
    for i=5:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-4) = beta_vec(1);
        Beta_Matrix(:,i-4) = beta_vec(2:end);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-4) + Beta_Matrix(1,i-4) * mean(excess_market_returns) + Beta_Matrix(2,i-4) * mean(small_big_returns) + Beta_Matrix(3,i-4) * mean(high_low_returns) + Beta_Matrix(4,i-4) * mean(factor1_data);
        ten_months_expected_returns(i-4) = expected_returns_on_asset * 210;  
        estimated_returns_of_asset = daily_risk_free_rate + alpha_vector(i-4) + Beta_Matrix(1,i-4) * excess_market_returns + Beta_Matrix(2,i-4) * small_big_returns + Beta_Matrix(3,i-4) * high_low_returns + Beta_Matrix(4,i-4) * factor1_data;
        rsquared(i-4) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-4) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-4)' * Sigma_F * Beta_Matrix(:,i-4)) + residual'*residual;
        variance_returns_asset(i-4) = var_asset_return;      
    end
    covariance = ones(col-4,col-4);
    
    for i=1:col-4
        for j = (i+1):col-4
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-4
        covariance(i,i) =  variance_returns_asset(i);
    end
    FF49IndustryPortfolios = covariance;
    FFF1ER = ten_months_expected_returns;
    FFF1VR =  variance_returns_asset;
    rSquarefrenchFama_F1Model = mean(rsquared);
end

function [FFF2VR, FFF2ER, FFTermSpread, rSquarefrenchFama_F2Model] = frenchFama_F2Model(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    excess_market_returns = data(:,1) - daily_risk_free_rate_vector;
    small_big_returns = data(:,2);
    high_low_returns = data(:,3);
    factor2_data = data(:,4);
    X_Factor = [excess_market_returns small_big_returns high_low_returns factor2_data];
    Sigma_F = cov(X_Factor);
    X_Factor = [constant X_Factor];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(4 ,col-4);
    alpha_vector = ones(col-4,1);
    expected_returns_asset = ones(col-4,1);
    variance_returns_asset = ones(col-4,1);
    ten_months_expected_returns = ones(col-4,1);
    rsquared = ones(col-4,1);
    for i=5:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-4) = beta_vec(1);
        Beta_Matrix(:,i-4) = beta_vec(2:end);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-4) + Beta_Matrix(1,i-4) * mean(excess_market_returns) + Beta_Matrix(2,i-4) * mean(small_big_returns) + Beta_Matrix(3,i-4) * mean(high_low_returns) + Beta_Matrix(4,i-4) * mean(factor2_data);
        ten_months_expected_returns(i-4) = expected_returns_on_asset * 210;  
        estimated_returns_of_asset = daily_risk_free_rate + alpha_vector(i-4) + Beta_Matrix(1,i-4) * excess_market_returns + Beta_Matrix(2,i-4) * small_big_returns + Beta_Matrix(3,i-4) * high_low_returns + Beta_Matrix(4,i-4) * factor2_data;
        rsquared(i-4) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        expected_returns_asset(i-4) = expected_returns_on_asset;
        var_asset_return = (Beta_Matrix(:,i-4)' * Sigma_F * Beta_Matrix(:,i-4)) + residual'*residual;
        variance_returns_asset(i-4) = var_asset_return;      
    end
    covariance = ones(col-4,col-4);
    
    for i=1:col-4
        for j = (i+1):col-4
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-4
        covariance(i,i) =  variance_returns_asset(i);
    end
    FFTermSpread = covariance;
    FFF2ER = ten_months_expected_returns;
    FFF2VR = variance_returns_asset;
    rSquarefrenchFama_F2Model = mean(rsquared);

end

function [FFF1F2VR, FFF1F2ER, FF49IndustryPortfoliosTermSpread, rSquarefrenchFama_F1_F2Model] = frenchFama_F1_F2Model(filename)
    data = xlsread(filename);
    [rows col] = size(data);
    constant = ones(rows,1);
    daily_risk_free_rate = 0.003;
    daily_risk_free_rate_vector = ones(rows,1) * daily_risk_free_rate;
    excess_market_returns = data(:,1) - daily_risk_free_rate_vector;
    small_big_returns = data(:,2);
    high_low_returns = data(:,3);
    factor1_data = data(:,4);
    factor2_data = data(:,5);
    X_Factor = [excess_market_returns small_big_returns high_low_returns factor1_data factor2_data];
    Sigma_F = cov(X_Factor);
    X_Factor = [constant X_Factor];
    t_X_Factor = X_Factor';
    Beta_Matrix = ones(5 ,col-5);
    alpha_vector = ones(col-5,1);
    expected_returns_asset = ones(col-5,1);
    variance_returns_asset = ones(col-5,1);
    ten_months_expected_returns = ones(col-5,1);
    rsquared = ones(col-5,1);
    for i=6:col
        actual_asset_returns = data(:,i) - daily_risk_free_rate_vector;
        beta_vec = inv(t_X_Factor*X_Factor)*t_X_Factor*actual_asset_returns;
        alpha_vector(i-5) = beta_vec(1);
        Beta_Matrix(:,i-5) = beta_vec(2:end);
        expected_returns_on_asset = daily_risk_free_rate + alpha_vector(i-5) + Beta_Matrix(1,i-5) * mean(excess_market_returns) + Beta_Matrix(2,i-5) * mean(small_big_returns) + Beta_Matrix(3,i-5) * mean(high_low_returns) + Beta_Matrix(4,i-5) * mean(factor1_data)+ Beta_Matrix(5,i-5) * mean(factor2_data);
        expected_returns_asset(i-5) = expected_returns_on_asset;
        ten_months_expected_returns(i-5) = expected_returns_on_asset * 210;  
        estimated_returns_of_asset = daily_risk_free_rate + alpha_vector(i-5) + Beta_Matrix(1,i-5) * excess_market_returns + Beta_Matrix(2,i-5) * small_big_returns + Beta_Matrix(3,i-5) * high_low_returns + Beta_Matrix(4,i-5) * factor1_data + Beta_Matrix(5,i-5) * factor2_data;
        rsquared(i-5) = calculateRSquared(actual_asset_returns, estimated_returns_of_asset);
        residual = actual_asset_returns - estimated_returns_of_asset; 
        var_asset_return = (Beta_Matrix(:,i-5)' * Sigma_F * Beta_Matrix(:,i-5)) + residual'*residual;
        variance_returns_asset(i-5) = var_asset_return;      
    end
    covariance = ones(col-5,col-5);
    
    for i=1:col-5
        for j = (i+1):col-5
            beta_i = Beta_Matrix(:,i);
            beta_j = Beta_Matrix(:,j);
            temp_cov = beta_i' * Sigma_F * beta_j;
            covariance(i,j) = temp_cov; 
        end
    end
    x = triu(covariance);
    covariance = x + x';
    for i=1:col-5
        covariance(i,i) =  variance_returns_asset(i);
    end
    FF49IndustryPortfoliosTermSpread = covariance;
    rSquarefrenchFama_F1_F2Model = mean(rsquared);
    FFF1F2ER = ten_months_expected_returns;
    FFF1F2VR = variance_returns_asset.*sqrt(210);
end

function [norm_mat] = compare_cov(cov_list)
    len = length(cov_list);
    euclidean_norm = zeros(len,len);
    
    for i=1:len-1
        cov_1 = cov_list(1,i);
        for j = (i+1):len
            cov_2 = cov_list(1,j);
            covar_1 = cell2mat(cov_1(1,1));
            covar_2 = cell2mat(cov_2(1,1));
            temp = norm((covar_1 - covar_2), 2);
            euclidean_norm(i,j) = temp; 
        end
    end
    
    t_euclidean_norm =  euclidean_norm';
    complete_matrix = t_euclidean_norm + euclidean_norm;
    
    for i =1:length(cov_list)
        complete_matrix(i,i) = 0.0;
    end
    norm_mat = complete_matrix;
end

function [omega_T] = tangencyPortfolioWeights(covar,muP)
    i_covar = inv(covar);
    mu_f = 0.0275;
    %mu_f = 0.0101;
    one = ones(length(muP),1);
    t_one = transpose(one);
    omegaBar = i_covar * ( muP -  mu_f * one );
    omega_T = omegaBar / (t_one * omegaBar);
end
