function [ mu_X, sigma_X ] = parameters_LN( mu_Y, sigma_Y )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

mu_X = exp(mu_Y+1/2*sigma_Y^2);
sigma_X = mu_X*sqrt(exp(sigma_Y^2)-1);

end

