function [ mean, std, skewness ] = moments_TypeIII(parGEV)
% moments of a dataset
ksi_Y = parGEV(1,1); 
sigma_Y = parGEV(2,1); 
mu_Y = parGEV(3,1);
mean = mu_Y + sigma_Y*((gamma(1-ksi_Y)-1)/ksi_Y);
std = sqrt(sigma_Y^2*((gamma(1-2*ksi_Y)-gamma(1-ksi_Y)^2)/ksi_Y^2));
skewness = -(gamma(1-3*ksi_Y)-3*gamma(1-ksi_Y)*gamma(1-2*ksi_Y)+2*gamma(1-ksi_Y)^3)/((gamma(1-2*ksi_Y)-gamma(1-ksi_Y)^2)^(3/2));

end




