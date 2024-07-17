function  P  = cdf_Gumbel_max(alpha1,u1,C)
%Gumbel distribution. Kans op een waarde kleiner dan C.    
for i=1:length(C)
       
P(i) = exp(-exp(-alpha1*(C(i)-u1)));

end
P=transpose(P);

