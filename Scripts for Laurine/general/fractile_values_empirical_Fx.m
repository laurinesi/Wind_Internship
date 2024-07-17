function [X_fractile ] = fractile_values_empirical_Fx( X, Fx, P_fractiles )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i = 1: length(P_fractiles)
    tmp = [];
    tmp = abs(Fx-P_fractiles(i)) ;
    [~, I] = min(tmp);
    X_fractile(i) = X(I);

end

