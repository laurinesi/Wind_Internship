function [ F_N_max ] = F_N_max( X )
% empirical distribution function (maxima)

F_N_max = [];
N = length(X);
for i = 1:N
    F_N_max(i,1) = i/(N+1);  % ik heb hier gedaan i/(n+1) -- Die man van Raphael gebruikt hier iets anders, namelijk (i+1.5)/(n+2.5) ofzo?
end 


end

