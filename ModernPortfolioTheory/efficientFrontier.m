bmu =[];
covar = [];
bmu_p = [];

filename1 = "AllStocks1Year.xlsx";
filename2 = "AllStocks2Years.xlsx";
filename3 = "AllStocks3Years.xlsx";
filename4 = "AllStocks1YearWithVXX.xlsx";

[covar1, muP1] = annulizeData(filename1);
plotEfficientFrontier(covar1,muP1);
tangencyPortfolioWeights(covar1,muP1);

[covar2, muP2] = annulizeData(filename2);
plotEfficientFrontier(covar2,muP2);
[covar3, muP3] = annulizeData(filename3);
plotEfficientFrontier(covar3,muP3);
[covar4, muP4] = annulizeData(filename4);
plotEfficientFrontier(covar4,muP4);
% figure
% plot(sigmaP1, muP1);
% hold on
% plot(sigmaP2, muP2);
% plot(sigmaP3, muP3);
% hold off

function [covar, muP] = annulizeData(filename)
    data = xlsread(filename);
%     dataset = []
%     for i = 1:3
%       records =  size(data,1)*i/3
%       dataset = data(1:end,1:20);
        [rows,col] = size(data);
        Rvec = ones(col,1);
        annulizeVec = ones(col,1);
        temp = ones(rows,1);
        R = ones(rows-1,col);
        for j=1:col
            temp = data(:,j);
            [logReturns ,test] = expectationOfLogReturns(temp);
            R(:,j)=logReturns;
            annulizeVec(j,1) = ((1+test/100)^365 - 1);
        end
        muP = annulizeVec;
        covar = cov(R);

%         for i=1:col
%             for j=1:rows
%                 temp(j,1) = data(j,i);
%             end
%             [logReturns ,Rvec(i,1)] = expectationOfLogReturns(temp);
%             R(:,i)=logReturns;
%             annulizeVec(i,1) = ((1+Rvec(i,1))^365 - 1);
%         end
%         muP = annulizeVec;
%         covar = cov(R);
%         %[sigmaP, muP] = plotEfficientFrontier(covar,bmu);
     %end
    end

function [optW,g,h] = optimalWeights(mu,covariance_mat,mu_p)
    one = ones(length(mu),1);
    if det(covariance_mat) == 0
        disp('Omega is not reversible!!!');
        return
    end
    i_covariance_mat = inv(covariance_mat);
    t_one = transpose(one);
    t_mu = transpose(mu);
    A = t_one * i_covariance_mat * mu;
    B = t_mu *  i_covariance_mat * mu;
    C = t_one * i_covariance_mat * one;
    D = B*C - A^2;
    g = B / D * i_covariance_mat*one - A / D * i_covariance_mat * mu;
    h = C / D * i_covariance_mat * mu - A / D * i_covariance_mat * one; 
    optW = g+mu_p*h;
end

function [sigma_p, mu_p] = plotEfficientFrontier(covariance_mat, bmu)
    mu_p = linspace(min(bmu),max(bmu),100);
    sigmaP = [];
    for i =1:length(mu_p)
        [optW_mu_p,g,h] = optimalWeights(bmu,covariance_mat,mu_p(i));
        sigmaP(i) =  sqrt(transpose(optW_mu_p) * covariance_mat * optW_mu_p);
    end
    sigma_p = sigmaP;
    h_t = transpose(h);
    g_t = transpose(g);
    gh = g_t * covariance_mat * h;
    hh = h_t * covariance_mat * h;
    gg = g_t * covariance_mat * g;
    mu_min = - gh/hh;
    sd_min  = sqrt( gg - gh^2/hh);
    efficient_mu =  [];
    efficient_sigma = [];
    above_min_index = []; 
    unefficient_mu =  [];
    unefficient_sigma = [];
    below_min_index = []; 
    j = 1;
    k = 1;
    %index = find(mu_p == mu_min)
    for i=1:length(mu_p)
        if mu_p(i) >  mu_min
            above_min_index(j) = i;
            efficient_mu(j) = mu_p(i);
            efficient_sigma(j)  = sigmaP(i);
            j =j+1;
        else 
            below_min_index(k) = i;
            unefficient_mu(k)  = mu_p(i);
            unefficient_sigma(k)  = sigmaP(i);
            k = k+1;
        end
    end
    
    figure
    plot(sigmaP(above_min_index), mu_p(above_min_index), 'r', 'Linewidth',2);  
    hold on
    xlabel('Risk in portfolio (Standard Deviation)')
    ylabel('Expected Returns of the Portfoloio');
    title('Efficient Frontier');
    plot(sigmaP(below_min_index), mu_p(below_min_index),'b--','Linewidth',2);
    plot(sd_min,mu_min,'k-x', 'Markersize', 20); 
    hold off
    bargraphs(bmu,covariance_mat,mu_min)
end

function bargraphs(bmu,covariance_mat,mu_min)
    [optW_mu_p] = optimalWeights(bmu,covariance_mat,mu_min);
    su = sum(optW_mu_p); 
    figure
    if(length(optW_mu_p) == 20)
        x = [1:20];
       labels = {"AAPL","ABBV", "COST", "CSCO","CVX","GE","IBM","JNJ","JPM","KO","MCD","MSFT","PEP","PFE","PG","SBUX","T","VZ","WFC","XOM"};
    else
        x = [1:21];
        labels = {"AAPL","ABBV","COST","CSCO","CVX","GE","IBM","JNJ","JPM","KO","MCD","MSFT","PEP","PFE","PG","SBUX","T","VZ","WFC","XOM","VXX"};
    end
    bar(optW_mu_p)
    %xt = get(gca, 'XTick');
    text(x,optW_mu_p, labels, 'rotation',90, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
    % set(gca, 'XTick', xt, 'XTickLabel', {"AAPL","ABBV", "COST", "CSCO","CVX","GE","IBM","JNJ","JPM","KO","MCD","MSFT","PEP","PFZR","PG","SBUX","T","VZ","WFC","XOM"});
    title('Weigths of the minimum variance portfolio')
end

function [R, mu] = expectationOfLogReturns(prices)
    logReturns = ones(length(prices)-1,1);
    for i = 1:(length(prices)-1)
        returns = (prices(i+1)/prices(i)) - 1;
        logReturns(i,1) = log(1+returns) * 100;
    end
    mu = (mean(logReturns));
    R = logReturns;
end

function tangencyPortfolioWeights(covar1,muP1)
i_covar = inv(covar1);
mu_f = 0.0288
one = ones(length(muP1),1);
t_one = transpose(one);
%for i =1:length(muP1)
omegaBar = i_covar * ( muP1 -  mu_f * one );
omega_T = omegaBar / (t_one * omegaBar);
end