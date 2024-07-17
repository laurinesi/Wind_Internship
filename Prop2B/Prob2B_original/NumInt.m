function [results, parametersOut] = NumInt(func, parameters, correlations,  settings, arg)
%
% Subroutine voor het uitvoeren van een numerieke integratie berekening
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
umin = settings.umin;
umax = settings.umax;
nstep = settings.nstep;
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
Pf   = 0;
u = zeros(nstoch,1);
nip = zeros(nstoch,1);
%
if ~isempty(correlations)
    trmH = getCorrTrans(correlations);
else
    trmH = eye(nstoch);
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

% Loops voor numerieke integratie
if (nstoch > 5), return; end;
if (nstoch >= 5) 
  n5 = nstep;
  if strcmp(parameters(5).type_v(1:3),'DET'), n5 = 1;end;
else 
  n5 = 1;
end
if (nstoch >= 4) 
  n4 = nstep;
  if strcmp(parameters(4).type_v(1:3),'DET'), n4 = 1;end;
else 
  n4=1;
end
if (nstoch >= 3) 
  n3 = nstep;
  if strcmp(parameters(3).type_v(1:3),'DET'), n3 = 1;end;
else 
  n3=1;
end
if (nstoch >= 2) 
  n2 = nstep;
  if strcmp(parameters(2).type_v(1:3),'DET'), n2 = 1;end;
else 
  n2 =1;
end
n1 = nstep;
if strcmp(parameters(1).type_v(1:3),'DET'), n1 = 1;end;
for i5 = 1:1:n5
  if (nstoch >= 5)  
     u(5) = umin + (umax - umin) * (2 * i5 - 1) / n5 / 2;
     p1 = 1-normcdf(-(umin + (umax - umin) * i5 / n5),0,1);         
     p2 = 1-normcdf(-(umin + (umax - umin) * (i5-1) / n5),0,1);         
     nip(5) = p1 - p2;
  end
  for i4=1:1:n4
      if (nstoch >= 4);
         u(4) = umin + (umax - umin) * (2 * i4 - 1) / n4 / 2;
         p1 = 1-normcdf(-(umin + (umax - umin) * i4 / n4),0,1);         
         p2 = 1-normcdf(-(umin + (umax - umin) * (i4-1) / n4),0,1);         
         nip(4) = p1 - p2;
      end
      for i3=1:1:n3
         if (nstoch >= 3);
             u(3) = umin + (umax - umin) * (2 * i3 - 1) / n3 / 2;
             p1 = 1-normcdf(-(umin + (umax - umin) * i3 / n3),0,1);         
             p2 = 1-normcdf(-(umin + (umax - umin) * (i3-1) / n3),0,1);         
             nip(3) = p1 - p2;
         end
         for i2=1:1:n2
             if (nstoch >= 2);
                 u(2) = umin + (umax - umin) * (2 * i2 - 1) / n2 / 2;
                 p1 = 1-normcdf(-(umin + (umax - umin) * i2 / n2),0,1);         
                 p2 = 1-normcdf(-(umin + (umax - umin) * (i2-1) / n2),0,1);         
                 nip(2) = p1 - p2;
             end
             for i1=1:1:n1
                 u(1) = umin + (umax - umin) * (2 * i1 - 1) / n1 / 2;
                 p1 = 1-normcdf(-(umin + (umax - umin) * i1 / n1),0,1);         
                 p2 = 1-normcdf(-(umin + (umax - umin) * (i1-1) / n1),0,1);         
                 nip(1) = p1 - p2;
%
% Z waarde bepalen
% Waarde z-functie bepalen
                 [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
                 Z = func(x, arg);
                 if logger == 1
                     fprintf(fileID, formatstr, [u' x' Z]);
                 end
%
% Totale bezwijkkans bepalen
                 if (Z < 0) 
                    Pf = Pf + prod(nip);
%
% Min waarde r en alfa vastleggen
                    rbeta = sqrt(sum(u.*u));
                    if (rbeta < rmin) 
                       rmin = rbeta;
                       alfan = -u / rbeta;
                    end
                 end
             end
         end
      end
  end
end
%
% Kans en Beta bepalen
results.Pf = Pf;
results.beta = norminv(1-Pf,0,1);
results.cnvg = cnvg;
if (Pf == 0) 
   alfan = -ones(nstoch,1)/sqrt(nstoch);
end
%
% Consistente set van U en X-waarden bepalen
u = -alfan * results.beta;
[u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
Z = func(x, arg);
if logger == 1
    fprintf(fileID, formatstr, [u' x' Z]);
end
results.Z = Z;
for i1=1:1:nstoch
   parametersOut(i1,1).u = u(i1);
   parametersOut(i1,1).x = x(i1);
   parametersOut(i1,1).alfa = alfan(i1);
end


if logger == 1
    fclose(fileID);
end

   
