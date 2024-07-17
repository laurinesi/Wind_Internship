function [Pfvv ierr] = Hodepo(betav, Pfu, rho, epsZ, epsH)           
%
% Subroutine Hohenbichler
%
% Parameters
%   betav   I  Kleinste Betrouwbaarheidsindex
%   Pfu     I  Kleinste waarde van de kans
%   Pfvv    O  Kans vv
%   rho     I  Correlatiecoefficient
%   ierr    O  Switch voor foutmelding
%
% TNO Bouw Sept 2006 SNH
% ---------------------------------------------------------------------- 
%
% Initialisatie 
   ierr   = 0;
   iterm  = 50;
   relaxf = 0.5;
   rexalf = 1.0 - relaxf;
   Q_new  = 0.0;
%
% Bij onvoldoende correlatie geen Hohenbichler nodig
   if (abs(rho) < 1e-8)
%       CALL QfromX (betav, P, Pfvv)
%       Pfvv = normcdf(-betav,0,1);
        [P Pfvv] = QfromX(betav);
      return
   end
%
% Constanten uitrekenen
   if (rho^2 >= 1)
      wortel = 0;
   else
      wortel = sqrt(1 - rho^2);
   end
   dZdW   = - wortel;
% 
% Loop aantal maal
   for nmaal = 1:1:3
      kosigZ = 1;
%
% Startwaarden
      Uster  = 0.0;
      Wster  = 0.0;
      dUster = -0.10;
%
% Nulde stap
%       CALL QfromX (-Uster, P, PhiU)
%       PhiU = normcdf(Uster,0,1);
        [P PhiU] = QfromX(-Uster);
      
      Psi = Pfu * PhiU;
%       CALL XfromQ (Psi, Uacc)
%       Uacc = norminv(1-Psi,0,1);
      Uacc = XfromQ (Psi);
      ZUstrv = -rho * Uacc;
      Uster = Uster + dUster;
%
% Loop aantal iteraties      
      for iter = 1:1:iterm
         Q_old = Q_new;
%          CALL QfromX (-Uster, P, PhiU)
%          PhiU = normcdf(Uster,0,1);
         [P, PhiU] = QfromX(-Uster);
         Psi = Pfu * PhiU;
%          CALL XfromQ (Psi, Uacc)
%          Uacc = norminv(1-Psi,0,1);
         Uacc = XfromQ (Psi);
         ZUster = -rho * Uacc;
         Zster  = betav + ZUster - wortel * Wster;
         dZdV = (ZUster - ZUstrv) / dUster;
         gemZ = Zster - dZdV * Uster - dZdW * Wster;
         sigZ = sqrt(dZdV^2 + dZdW^2);
         if (sigZ == 0.0)
            kosigZ = 2;
            break;
         end
         beta = gemZ / sigZ;
%          CALL QfromX (beta, P, Q_new)
%          Q_new = normcdf(-beta,0,1);
         [P, Q_new] = QfromX(beta);
%
% Convergentie controleren
         if (abs(Zster) < epsZ)
            if (abs(Zster / sigZ) < 1e-2)
               if (abs(Q_old) > 1e-30)
                  if (abs(1 - Q_new / Q_old) < epsH)
                      Pfvv = Q_new;
                      return;
                  end
               else
                  if (abs(Q_new - Q_old) < epsH)
                      Pfvv = Q_new;
                      return;
                  end
               end
            end
         end
%
% Volgende iteratie voorbereiden
         ZUstrv = ZUster;
         Usterv = Uster;
         Ustern = - dZdV * beta / sigZ;
         Uster  = relaxf * Ustern + rexalf * Usterv;
         dUster = Uster - Usterv;
         Wster  = - dZdW * beta / sigZ;
%
% Toetsen voldoende aantal iteraties      
         if (abs(dUster) < 1e-30) 
            Pfvv = Q_new;
            return
         end
%
% Volgende iteratie      
      end
      relaxf = 0.5D0 * relaxf;
      rexalf = 1.0D0 - relaxf;
      iterm  = 2 * iterm;
   end
%
% Foutmeldingen
   if (kosigZ == 2)
      ierr = 1;
      disp (['De st. afw. van Z is 0 na ' num2str(iter) ' iter. in subroutine HODEPO']);
   else
      ierr = 0;
   end
%
% Kans vv bepalen 
  Pfvv = Q_new;
