function [ mean, std, skewness ] = moments_TypeI( alpha,u )
% moments of a dataset

mean = u + 0.577/alpha;
std = pi/alpha*1/sqrt(6);
skewness = 1.14;


end

