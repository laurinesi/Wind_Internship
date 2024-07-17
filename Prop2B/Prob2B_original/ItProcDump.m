function [itsom, inorde, alfan, beta, gemZ, Pf, sigZ, x, u, z] = ...
                ItProcDump(func, iterm, itsom, nstoch, parameters, trmH, ...
                       relaxf, u, epsB, epsZ, arg)
%                                                   
% Subroutine voor het uitvoeren van het iteratieproces
%
% Parameters
%   func    I  handle naar z-functie
%   iterm   I  Maximum aantal iteraties
%   itsom  I/O Som iteraties
%   inorde  O  Switch voor convergentie
%   nstoch  I  Aantal stochasten
%   alfan   O  Alfa waarden
%   beta    O  Betrouwbaarheidsindex
%   Pf      O  Kans
%   relaxf  I  Relaxatie factor
%   sigX    I  Standaard afwijking X
%   sigZ    O  Standaard afwijking Z
%   x       O  Waarden X
%   u      I/O Vector van U-waarden
%   z       O  Waarde Z-functie
%   parameters  I  Verdelingstype en data
%   trmH    I   transformation matrix for correlated parameters
%   arg     I  free form: user parameters voor z-functie en BepXfU
%
% TNO Bouw Sept 2006 SNH
% ----------------------------------------------------------------------
%                   
% Initialisatie 
itsom = 0;
inorde = 0;
dU   = 0.3D0;
Umat = zeros(iterm,nstoch);
dZdUmat = zeros(iterm,nstoch);
Zmat = zeros(iterm,1);
gemZmat = zeros(iterm,1);
sigZmat = zeros(iterm,1);
%      
% Iteraties uitvoeren. Maximaal iterm maal
dZdU = zeros(nstoch,1);
for iter = 1:1:iterm
% Totaal aantal iteraties      
    itsom = itsom + 1;
disp([' iteration ' num2str(itsom)]);
    
% Afgeleiden (helling) z-functie bepalen voor iedere stochast
    for k = 1:1:nstoch
         u(k) = u(k) - dU / 2.D0;
         [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
         for l=1:nstoch
             disp(['istoch ' num2str(l) ' u ' num2str(u(l)) ' x ' num2str(x(l))]);
         end
         Zm = func(x, arg);
         disp([' Zmin ' num2str(Zm)]);
         u(k) = u(k) + dU;
         [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
         for l=1:nstoch
             disp(['istoch ' num2str(l) ' u ' num2str(u(l)) ' x ' num2str(x(l))]);
         end
         Zp = func(x, arg);
         disp([' Zplus ' num2str(Zp)]);
         dZdU(k,1) = (Zp - Zm) / dU;
         u(k) = u(k) - dU / 2.D0;
    end
     for l=1:nstoch
         disp(['istoch ' num2str(l) ' dZdU ' num2str(dZdU(l))]);
     end
% Waarde z-functie bepalen
    [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
     for l=1:nstoch
         disp(['istoch ' num2str(l) ' u ' num2str(u(l)) ' x ' num2str(x(l))]);
     end
    z = func(x, arg);
    disp([' z ' num2str(z)]);
% Gemiddelde waarde Z bepalen
    gemZ = z - dZdU'*u;
    disp([' gemZ ' num2str(gemZ)]);
% Standaardafwijking Z bepalen
    sigZ = dZdU'*dZdU;
    disp([' sigZ ' num2str(sigZ)]);
    if (sigZ<=0.)
        sigZ = 0.D0;
    else
        sigZ = sqrt(sigZ);
    end
    if (abs(sigZ)<1.e-8)
        inorde = 0;
        break
    end;
% Beta bepalen
    beta = gemZ / sigZ;
    disp([' beta ' num2str(beta)]);
    
% Alfan bepalen
    alfan = dZdU/sigZ;
     for l=1:nstoch
         disp(['istoch ' num2str(l) ' alfan ' num2str(alfan(l))]);
     end
    
% Kans bepalen 
%     Pf = normcdf(-beta,0,1);
    [P Pf] = QfromX(beta);
% Toets op convergentie
    conv1 = (abs(z/sigZ)<epsZ);
    if (conv1)
        beta2 = beta^2;
        som = u'*u;
        conv2 = ((beta2>=(1.-epsB)*som) && (beta2<=(1.+epsB)*som));
        if (conv2)
            conv3 = (abs(beta)<19.999);
            if (conv1&&conv2&&conv3)
                inorde = 1;
                break;
            end
        end
    end
    Un = -beta.*alfan;
    u = relaxf.*Un + (1-relaxf).*u;

Umat(iterm,:) = u';
dZdUmat(iterm,:) = dZdU';
Zmat(iterm) = z;
gemZmat(iterm) = gemZ;
sigZmat(iterm) = sigZ;
    
% Einde loop
end
%
% Beta consistent met u-waarden maken bij niet convergentie
if inorde==0
    beta2 = u'*u;
    if (gemZ>0.)
      beta = sqrt(beta2);
    else
      beta = -sqrt(beta2);
    end
%     Pf = normcdf(-beta,0,1);
    [P Pf] = QfromX(beta);
end
% X, U en alfa waarden 100% consistent maken (ook bij convergentie)
u = -beta*alfan;
[u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
z = func(x, arg);
      