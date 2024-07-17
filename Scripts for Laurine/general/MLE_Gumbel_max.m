function [ alpha, u, nlogL ] = MLE_Gumbel_max( X )
% MLE estimation of model parameters of the Gumbel distribution (maxima)

parmhat = evfit(-X);   %This form of the probability density function is suitable for modeling the minimum value. To model the maximum value, use the negative of the original values.
alpha = 1/parmhat(2);
u = -parmhat(1);


nlogL = evlike(parmhat,-X);



end

