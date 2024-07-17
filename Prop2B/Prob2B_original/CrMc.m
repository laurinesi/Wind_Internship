function [results parametersOut] = CrMc(func, parameters, correlations, settings, arg)
%
% Subroutine voor het uitvoeren van een Monte Carlo berekening
%
% Parameters
%   func           I  function handle naar z-functie
%   parameters     I  structure: verdelingstype en data
%   settings       I  structure: specifieke settings voor methode
%   arg            I  free form: user parameters voor z-functie en BepXfU
%   parametersOut  O  structure: verdelingstype en data, u , x, alfa
%   results        O  structure: Pf, beta, cnvg
%
% TNO Bouw Sept 2006 SNH
% ----------------------------------------------------------------------
%
iseed = settings.iseed;
DSmxsp = settings.nmax;
MCVarPf = settings.varpfail;
MCVarNP = settings.varpnonfail;
if (isfield(settings, 'logger'))
    logger = settings.logger;
else
    logger = 0;
end

nstoch = size(parameters,1);
parametersOut = parameters;
%
% Initialiseren string voor convergentie en loops
cnvg = ' ';
rmin = 200.0;
ntot = 0;
%

if ~isempty(correlations)
    trmH = getCorrTrans(correlations);
else
    trmH = eye(nstoch);
end


% Initialisatie random generator
if (iseed == 0) 
    rand('twister', sum(100*clock));
else
    rand('twister', iseed);
end
%
if logger == 1
    filename = [pwd '\logger.txt'];
    fileID = fopen(filename,'w');
    formatstr = '%12.4e';
    for ii = 2:2*nstoch+1
        formatstr = [formatstr ', %12.4e'];
    end
    formatstr = [formatstr '\n'];
else
    fileID = [];
end

% Loop voor aantal samples
for nmaal = 1:1:DSmxsp
%
% U-waarden trekken
    p = rand(nstoch,1);
    u = norminv(p,0,1);
%
% Z waarde bepalen
    [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
    z = func(x, arg);
    if logger == 1
        fprintf(fileID, formatstr, [u' x' z]);
    end
%
% Totale bezwijkkans bepalen
    if (z < 0.) 
         ntot = ntot + 1;
%
% Min waarde r en alfa vastleggen
	     rbeta = u'*u ;
		 rbeta = sqrt(rbeta);
         if (rbeta < rmin) 
            rmin = rbeta;
            alfan = -u./rbeta;
         end
    end
    
     Pf = 1.0 * ntot / nmaal;
     if (Pf>0.) && (Pf < 1.0)
        varPf = sqrt(1.0 / nmaal * (1.0 / Pf - 1.0));
        varNP = sqrt(1.0 / nmaal * (1.0 / (1.0 - Pf) - 1.0));
        if (varPf<MCVarPf) && (varNP<MCVarNP)
            break;
        end
     end
%
% Einde loop voor samples
end
%
% % Kans en Beta bepalen
% Pf = 1.0 * ntot / nmaal;
% beta = norminv(1-Pf,0,1);
% if (ntot == 0) 
%      alfan=-sqrt(1/nstoch)*ones(nstoch,1);
% end
% % Consistente set van U en X-waarden bepalen
% u = -alfan * beta;
% [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
% z = func(x, arg);

if (ntot == 0) 
    % Kans en Beta bepalen
    Pf = 0;
    beta = 8;
    alfan=-sqrt(1/nstoch)*ones(nstoch,1);
else
    % Kans en Beta bepalen
    Pf = 1.0 * ntot / nmaal;
    beta = norminv(1-Pf,0,1);
end
% Consistente set van U en X-waarden bepalen
u = -alfan * beta;
[u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
z = func(x, arg);
if logger == 1
    fprintf(fileID, formatstr, [u' x' z]);
end


% Resultaten
results.Pf = Pf;
results.beta = beta;
results.cnvg = cnvg;
results.Z = z;
results.N = nmaal;

for i1=1:1:nstoch
   parametersOut(i1,1).u = u(i1);
   parametersOut(i1,1).x = x(i1);
   parametersOut(i1,1).alfa = alfan(i1);
end

if logger == 1
    fclose(fileID);
end
