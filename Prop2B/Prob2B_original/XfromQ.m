function [X] = XfromQ (Q)
%
% Subroutine voor het berekenen van de betrouwbaarheidsindex bij een
% gegeven overschrijdingskans voor een normale verdeling.
% Iteratief gecorrigeerd met subroutine QfromX.
%      
% Parameters
%   Q       I  Overschrijdingskans
%   X       O  Betrouwbaarheidsindex (beta)
%
% TNO Bouw Sept 2006 SNH
% ----------------------------------------------------------------------
%
% Parameters
%
% Afbreekcriterium en maximum aantal iteraties
   eps    = 1.D-20;
   k2max  = 25;
%
% Overgaan op overschrijdingskansen tussen 0 en 0.5
   if (Q > 0.5D0) 
      Qsr = 1.0D0 - Q;
   else
      Qsr = Q;
   end
%
% Beginschatting voor Xpa
   if (Qsr < 0.15D0) 
      if (Qsr < 1.0D-300), Qsr = 1.0D-300; end;
      X = sqrt(-2.0 * log(5.0 * Qsr));
   else
      X = 1.464795D0 - 2.929590D0 * Qsr;
   end
%
% Beschouw ligging van schatting
%    CALL QfromX(X, P, Qs)
   [P Qs] = QfromX(X);
   if (Qs > Qsr) 
      X1  = X;
      Q1  = Qs;
      koQ = 1;
      X   = X + 0.1D0;
   else
      X2  = X;
      Q2  = Qs;
      koQ = 2;
      X   = X - 0.1D0;
   end
%
% Bepaal punt aan andere zijde van de correcte waarde
   success = false;
   while (~success)
       [P Qs] = QfromX(X);
       if (Qs > Qsr) 
          X1 = X;
          Q1 = Qs;
          if (koQ == 2) 
              success=true;
          else
              X = X + 0.1;
          end
       else
          X2 = X;
          Q2 = Qs;
          if (koQ == 1) 
              success = true;
          else
              X = X - 0.1D0;
          end
       end
   end
%
% Benader X d.m.v. lineaire interpolatie
   k2 = 0;
   success = false;
   while ~success
       k2 = k2 + 1;
       if (k2 > (2 * k2max)) 
           success = true;
       else
           if (k2 > k2max) 
              if (abs(Qs - Qsr) < eps), success=true; end;
              eps = 10.0D0 * eps;
              if (eps > 1.0D-8), eps = 1.0D-8; end;
           end
           if (Q1 == Q2), success=true; end;
           if (~success)
               X = X1 + (Q1 - Qsr) / (Q1 - Q2) * (X2-X1);
               [P Qs] = QfromX(X);
               if (Qs > Qsr) 
                  X1 = X;
                  Q1 = Qs;
               else
                  X2 = X;
                  Q2 = Qs;
               end
           end
       end
   end
%
% Terugtransformatie indien overschrijdingskans > 0.5
   if (Q > 0.5D0), X = -X; end;

