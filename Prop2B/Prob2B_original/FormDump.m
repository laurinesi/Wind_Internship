function [results, parametersOut] = FormDump(func, parameters, correlations, settings, arg)
%
% Subroutine voor het uitvoeren van een FORM berekening
%
% Parameters
%   func           I  function handle naar z-functie
%   parameters     I  structure: verdelingstype en data
%   correlations   I  correlation matrix for parameters
%   settings       I  structure: specifieke settings voor methode
%   arg            I  free form: user parameters voor z-functie en BepXfU
%   parametersOut  O  structure: verdelingstype en data, u , x, alfa
%   results        O  structure: Pf, beta, cnvg
%
% TNO Bouw Sept 2006 SNH
% ----------------------------------------------------------------------
%
itmax = settings.itmax;
relaxf = settings.relaxf;
epsB = settings.epsB;
epsZ = settings.epsZ;
nloop = settings.maxloop;
nstoch = size(parameters,1);
parametersOut = parameters;
%
% Startwaarden bepalen
[Uopt] = FrmIni(func, parameters, arg);
%
if ~isempty(correlations)
    trmH = getCorrTrans(correlations);
else
    trmH = eye(nstoch);
end
% Initialiseren string voor convergentie en loops
cnvg   = ' ';
itsom  = 0;
iterm  = itmax;
%
% Loop voor 3 pogingen met verschillend max. aantal iteraties en relaxatiefactor
for nmaal = 1:1:nloop
%
    disp([' iloop ' num2str(nmaal) ' itmax ' num2str(iterm) ' relaxf ' num2str(relaxf)]);
% Startwaarden copieren     
    u = Uopt;
    [itsom, inorde, alfan, beta, gemZ, Pf, sigZ, x, u, z] = ...
    ItProcDump(func, iterm, itsom, nstoch, parameters, trmH, relaxf, ...
    u, epsB, epsZ, arg);
%
% Copieren data voor iteratiemethode 3     
    UMeth3 = u;
    iMeth3 = inorde;
%
% Controleren resultaten iteratieproces        
    if (inorde==1), break; end;
%
% Nog eens proberen met nieuwe waarden
    relaxf = 0.5D0 * relaxf;
    iterm  = 2 * iterm;
end
%
% Foutmelding indien er geen convergentie is
if (inorde~=1)
    disp(['Form  GEEN CONVERGENTIE']);
    cnvg = '*';
end
%
% Alfawaarden normeren
somalf = alfan'*alfan;
if (somalf>0.)
    somalf = sqrt(somalf);
    alfan = alfan/somalf;
    u = u/somalf;
end

% Resultaten
results.Pf = Pf;
results.beta = beta;
results.cnvg = cnvg;
results.Z = z;
results.sigZ = sigZ;
for i1=1:1:nstoch
   parametersOut(i1,1).u = u(i1);
   parametersOut(i1,1).x = x(i1);
   parametersOut(i1,1).alfa = alfan(i1);
end
   



