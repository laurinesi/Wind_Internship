function [alfa, beta, Pf] = CombiCor (cfunc, nstoch, alfa1, beta1, Pf1, ...
                                                alfa2, beta2, Pf2, rho12, epsZ, epsH)
%
% Subroutine voor het combineren van de twee mechanismen
% Er wordt aangenomen dat stochasten met gelijk volgnummer een gelijke betekenis hebben
% en er correlatie bestaat. De waarde hiervan wordt berekend op basis van
% de alfa's en een op te geven correlatie
% Parameter
%   cfunc   I  'And' of 'Or' kans (string)
%   nstoch  I  Aantal stochasten
%   alfa1   I  Alfa waarden mechanisme 1
%   beta1   I  Betrouwbaarheidsindex mechanisme 1
%   Pf1     I  Kans mechanisme 1
%   alfa2   I  Alfa waarden mechanisme 2
%   beta2   I  Betrouwbaarheidsindex mechanisme 2
%   Pf2     I  Kans mechanisme 2
%   rho12   I  Correlatie tussen de stochasten van mecahnisme 1 en 2
%   alfa    O  Alfa waarden gecombineerd mechanisme
%   beta    O  Betrouwbaarheidsindex gecombineerd mechanisme
%   Pf      O  Kans gecombineerd mechanisme
%
% TNO Bouw Sept 2006 SNH
% ---------------------------------------------------------------------- 
%
    alfa = zeros(nstoch,1);
% Correlatie uitrekenen
% De formule voor de correlatie is rho = som i = 1 tot n [alfa i*alfa i]
   rho = sum(alfa1.*alfa2.*rho12);
   if (rho >= 0.999999), rho = 0.999999; end;
%
% Aanroepen procedure Hohenbichler voor bepaling betrouwbaarheidsindex
   if (beta1 >= beta2) 
      Pfu   = Pf1;
      betav = beta2;
   else
      Pfu   = Pf2;
      betav = beta1;
   end
   if (beta1 > 37.05) && (beta2 > 37.05)
      beta = beta1;
      Pf   = Pf1;
      alfa = alfa1;
      return;
   end
   [Pfvv ierr] = Hodepo(betav, Pfu, rho, epsZ, epsH);           
   if (ierr > 0), return; end;
%
% Faalkans berekenen
   if (strcmp(upper(strtrim(cfunc)),'OR'))
      Pf = Pf1 + Pf2 - Pfu * Pfvv;
   else
      Pf = Pfu * Pfvv;
   end
%
% Betrouwbaarheidsindex bepalen
%    CALL XfromQ (Pf, beta)
%     beta = norminv(1-Pf,0,1);
   beta = XfromQ(Pf);
%
% Aanroepen procedure Hohenbichler voor bepaling van alfa waarden
for i = 1:1:nstoch
    %%% gecorreleerde deel
      beta1k = beta1 - alfa1(i) * 0.01;
%       CALL QfromX (beta1k, P, Pf1k)
%       Pf1k = normcdf(-beta1k,0,1);
        [P, Pf1k] = QfromX(beta1k);
        
      beta2k = beta2 - alfa2(i) * 0.01*rho12(i);
%       CALL QfromX (beta2k, P, Pf2k)
%       Pf2k = normcdf(-beta2k,0,1);
        [P, Pf2k] = QfromX(beta2k);
%
% Grootste waarde beta1k of beta2k
      if (beta1k >= beta2k)
         Pfu   = Pf1k;
         betav = beta2k;
      else
         Pfu   = Pf2k;
         betav = beta1k;
      end
      [Pfvv, ierr] = Hodepo(betav, Pfu, rho, epsZ, epsH);           
      if (ierr > 0) , return; end;
%
% Faalkans Pfk berekenen
      if (strcmp(upper(strtrim(cfunc)),'OR'))
         Pfk = Pf1k + Pf2k - Pfu * Pfvv;
      else
         Pfk = Pfu * Pfvv;
      end
%
% Betrouwbaarheidsindex bepalen
%       CALL XfromQ (Pfk, betak)
%       betak = norminv(1-Pfk,0,1);
       betak = XfromQ(Pfk);
%
% Alfa1 bepalen
      alfax1 = (beta - betak) / 0.01;
        
    % % % Ongecorreleerde deel
      if (rho12(i)*rho12(i) <= 1.0)
	     beta1k = beta1;
%          CALL QfromX (beta1k, P, Pf1k);
%          Pf1k = normcdf(-beta1k,0,1);
         [P, Pf1k] = QfromX(beta1k);
         beta2k = beta2 - alfa2(i) * 0.01 * sqrt(1.0-rho12(i)*rho12(i));
%          CALL QfromX (beta2k, P, Pf2k);
%          Pf2k = normcdf(-beta2k,0,1);
         [P, Pf2k] = QfromX(beta2k);

%  Grootste waarde beta1k of beta2k
         if (beta1k >= beta2k)
            Pfu   = Pf1k;
            betav = beta2k;
         else
            Pfu   = Pf2k;
            betav = beta1k;
         end

%  Aanroepen procedure Hohenbichler voor bepaling invloedscoefficient
         [Pfvv, ierr] = Hodepo(betav, Pfu, rho, epsZ, epsH);           
         if (ierr > 0), return; end;
%          CALL Hodepo (betav, Pfu, Pfvv, rho, ferr, werr);
% 	     if (ferr .NE. 0) RETURN

%  Faalkans Pfk berekenen
         if (strcmp(upper(strtrim(cfunc)),'OR'))
            Pfk = Pf1k + Pf2k - Pfu * Pfvv;
         else
            Pfk = Pfu * Pfvv;
         end 
 
%   Betrouwbaarheidsindex bepalen
%          CALL XfromQ (Pfk, betak)
%          betak = norminv(1-Pfk,0,1);
         betak = XfromQ(Pfk);

%  Alfa2 bepalen
         alfax2 = (beta - betak) / 0.01;
      else
	     alfax2 = 0.00;
      end

%  Alfa bepalen
      alfa(i) = alfax1*alfax1+alfax2*alfax2;
      if (alfa(i) > 0.)
         if (abs(alfax1) > abs(alfax2))
            alfa(i) = sqrt(alfa(i))*sign(alfax1);
         else
            alfa(i) = sqrt(alfa(i))*sign(alfax2);
         end
      end
     
end
%
% Alfa normeren
   som = alfa'*alfa;
   if (som > 0)
      som = sqrt(som);
      alfa = alfa./som;
   else
      alfa = ones(nstoch,1)./sqrt(nstoch);
   end
%
