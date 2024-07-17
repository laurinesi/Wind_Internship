function [results] = OmTijd (iopt, settings, nrgt, alfa, Pfin, beta, inrhot)
% 
% Subroutine voor het in rekening brengen van de tyd
%
% Parameters
%   nrgt           I  Aantal getijen
%   alfa           I  Alfa waarden
%   beta           I  Betrouwbaarheidsindex
%   Pf             I  Faalkans
%   inrhot         I  Correlatie in de tijd
%   settings       I  Structure: specifieke settings voor methode
%   results        O  Structure: Pf, beta, cnvg
%      
% TNO Bouw Sept 2006 SNH
% ----------------------------------------------------------------------
%
epsH = settings.epsH;
epsZ = settings.epsZ;
% Correlatie bepalen
rhot = sum(alfa.*alfa.*inrhot);
%
% Kans bepalen
if(iopt == 1 && rhot > 0.001)
  if (rhot > 0.9999999999),	rhot = 0.9999999999; end
  uonder = (beta - 8.5D0 * sqrt(1.-rhot)) / sqrt(rhot);
  fac = 1.e-15;
  udelta = Pfin * fac;
  uboven = XfromQ(udelta);
% %   uboven = norminv(1-udelta,0,1);
%   while uboven==Inf
%       fac = 10*fac;
%       udelta = Pfin * fac;
% %       uboven = norminv(1-udelta,0,1);
%         uboven = XfromQ(udelta);
%   end
%   uboven = max(uboven,beta);
%   if (uonder>=uboven)
%       if (uboven>beta)
%           uonder = beta - (uboven-beta);
%       else
%           uonder = beta - 1;
%       end
%   end
      

%   uonder = -8;
%   uboven = 8;
  udelta = (uboven - uonder) / 1000;
  Pf = 0.0D0;
  for i = 1:1:1001
     u = uonder + udelta * (i-1);
     beta_st = (beta -  sqrt(rhot) * u) / sqrt(1.0D0-rhot);
%          CALL QfromX (x, p, q)
     [p q] = QfromX(beta_st);
     Pf = Pf + (1 - p^nrgt) * (exp(-u * u / 2.0D0) / sqrt(2 * pi)) * udelta;
  end
else
%        CALL Hodepo (beta, Pf, Pfvv, rhot, ierr)
   [Pfvv, ierr] = Hodepo(beta, Pfin, rhot, epsZ, epsH);
%    [Pfvv ierr] = Hodepo2(beta, Pfin, rhot, epsZ, epsH);           
   if (ierr > 0) , return; end;
   Pf = Pfin + (nrgt - 1.0D0) * (Pfin - Pfvv * Pfin);
   if (Pf > 1.0), Pf = 1.; end;
   
%     m.beta = beta;
%     m.rho = rhot;
%     m.p = Pfin;
%     m.ngetij = 60;
%     
%     par = struct( ...
%     'names', {'u'       ;'w'     }, ...
%     'gemX',  {0         ;0}, ...
%     'sigX',  {1         ;1}, ...
%     'data',  {[]        ;[]}, ...
%     'type_v',{'NOR     ';'NOR     '}, ...
%     'x',     {0         ;0}, ...
%     'u',     {0         ;0}, ...
%     'rhot',  {0         ;0}, ...
%     'alfa',  {0         ;0});
% 
%     form_settings = struct('itmax',50,'relaxf',0.25,'epsB',0.0001,'epsZ',0.0001);
%     
%     [FormR FormP] = FormAPT1(@ZHohen, par, [], m, form_settings);
%     Pf = Pfin + (nrgt - 1.0D0) * (Pfin - FormR.Pf * Pfin);
%     if (Pf > 1.0), Pf = 1.; end;
   
   
end

if (Pf<Pfin)
    %    CALL XfromQ(Pf, beta)
    results.Pf = Pfin;
    results.beta = beta;
    results.cnvg = '';
    results.Z = NaN;
else
    %    CALL XfromQ(Pf, beta)
    %    beta = norminv(1-Pf,0,1);
    beta = XfromQ(Pf);
    results.Pf = Pf;
    results.beta = beta;
    results.cnvg = '';
    results.Z = NaN;
end
