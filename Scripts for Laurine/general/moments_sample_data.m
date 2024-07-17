function [ mu, sigma, a ] = moments_sample_data( X )
% moments of a dataset

mu = mean(X);
sigma = std(X,0); %unbiased
a = skewness(X,0); %unbiased


end

