function [P, Q] = QfromX (X)
% 
% Subroutine voor het berekenen van de over/onderschrijdingskans bij een
% gegeven betrouwbaarheidsindex voor een normale verdeling.
% Approximation Formula 26.2.17 for Normal Probability Function,
% Handbook of Mathematical Functions, Abramowitz & Stegun, page 932
%     
% Parameters
%   X       I  Betrouwbaarheidsindex (beta)
%   P       O  Kans
%   Q       O  Overschrijdingskans
%
% TNO Bouw Sept 2006 SNH
% ----------------------------------------------------------------------
%
% Constanten
   r  =  0.2316419D0;
   b1 =  0.319381530D0;
   b2 = -0.356563782D0;
   b3 =  1.781477937D0;
   b4 = -1.821255978D0;
   b5 =  1.330274429D0;
%
   ntal = length(X);
   if (ntal>1)
        P = zeros(size(X));
        Q = zeros(size(X));
        %
        % Kansdichtheidsfunctie
           z = exp(-0.50 * X.*X) / 2.50662827463100D0;
        %
        % Kans en Overschrijdingskans
           t = 1.0./(1.0 + r * abs(X));
           b = t.*(b1 + t.*(b2 + t.*(b3 + t.*(b4 + t.*b5))));
           indices = find(X >= 0.0);
           if size(indices,1)>0
              Q(indices(:,1)) = z(indices(:,1)).*b(indices(:,1));
              P(indices(:,1)) = 1.0 - Q(indices(:,1));
           end
           indices = find(X < 0.0);
           if size(indices,1)>0
              P(indices(:,1)) = z(indices(:,1)).*b(indices(:,1));
              Q(indices(:,1)) = 1.0 - P(indices(:,1));
           end
        %
   else
        %
        % Kansdichtheidsfunctie
           z = exp(-0.50 * X * X) / 2.50662827463100D0;
        %
        % Kans en Overschrijdingskans
           t = 1.0D0 / (1.0D0 + r * abs(X));
           b = t * (b1 + t * (b2 + t * (b3 + t * (b4 + t * b5))));
           if (X >= 0.0D0)
              Q = z * b;
              P = 1.0D0 - Q;
           else
              P = z * b;
              Q = 1.0D0 - P;
           end
        %
   end
