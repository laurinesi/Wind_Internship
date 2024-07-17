function [results, parametersOut] = Sorm(func, parameters, correlations, settings, arg)
%
% Subroutine voor het uitvoeren van een SORM berekening
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
sormoption = 1;
itmax = settings.itmax;
relaxf = settings.relaxf;
epsB = settings.epsB;
epsZ = settings.epsZ;
nstoch = size(parameters,1);
parametersOut = parameters;
if ~isempty(correlations)
    trmH = getCorrTrans(correlations);
else
    trmH = eye(nstoch);
end

[FormRes, FormPar] = Form(func, parameters, correlations,settings, arg);
%
% Foutmelding indien er geen convergentie is
if (strcmp(FormRes.cnvg,'*'))
    disp(['Form  GEEN CONVERGENTIE']);
else
%
end

u = zeros(nstoch,1);
x = zeros(nstoch,1);
alfa = zeros(nstoch,1);
for i1=1:1:nstoch
   u(i1) = FormPar(i1,1).u;
   x(i1) = FormPar(i1,1).x;
   alfa(i1) = FormPar(i1,1).alfa;
end
beta = FormRes.beta;
Pf = FormRes.Pf;
sigZ = FormRes.sigZ;

% % if ((alfa(1)~=0)||(alfa(nstoch)~=0))
% %check exceedance prob of Q*v
%     par(1,1) = parameters(1);
%     par(2,1) = parameters(nstoch);
%     meta.Qv = x(1)*x(nstoch);
%     meta.ngetij = metadata.ngetij;
%     [QvRes QvPar] = FormAPT1(@Qv, par, [], meta, settings);
%     disp([ 'P(Qv>Qv_d) : ' num2str(QvRes.Pf)]);
% %end




% Rotatie-Transformatie matrix D bepalen
   IX1 = linspace(1,nstoch,nstoch)';
   D1 = zeros(nstoch,nstoch);
% Eerste kolom
   D1(:,1) = -alfa;
% Andere kolommen
   for k = 2:1:nstoch
      D1(k,k) = 1.0;
      for l = 1:1:k-1
         inpr = D1(:,k)'* D1(:,l);
         D1(:,k) = D1(:,k) - inpr * D1(:,l);
      end 
      somalf = sum(D1(:,k).* D1(:,k));
      if (somalf > 0.0) 
         somalf = somalf^0.5;
         D1(:,k) = D1(:,k) / somalf;
      else
         werr = 600;
      end
   end


   [B,IX2] = sort(abs(alfa),1,'descend');
   D2 = eye(nstoch);
   for k = 2:1:nstoch
       if (alfa(IX2(k))~=0)
           T = eye(nstoch);
           cs = sqrt(alfa(IX2(1:k-1))'* alfa(IX2((1:k-1))));
           if (cs==0)
               gamma = sign(-alfa(IX2(k)))*pi/2;
           else
               gamma = atan(-alfa(IX2(k))/cs);
           end
%             if (IX2(k)>IX2(1))
               T(IX2(1),IX2(1)) =  cos(gamma);
               T(IX2(1),IX2(k)) = -sin(gamma);
               T(IX2(k),IX2(1)) =  sin(gamma);
               T(IX2(k),IX2(k)) =  cos(gamma);
%            else
%                T(IX2(1),IX2(1)) =  cos(gamma);
%                T(IX2(1),IX2(k)) =  sin(gamma);
%                T(IX2(k),IX2(1)) = -sin(gamma);
%                T(IX2(k),IX2(k)) =  cos(gamma);
%            end
           D2 = D2*T;
       end
   end

   IX3 = linspace(1,nstoch,nstoch)';
   D3 = eye(nstoch);
   for k = 2:1:nstoch
       if (alfa(k)~=0)
           T = eye(nstoch);
           cs = sqrt(alfa(1:k-1)'* alfa((1:k-1)));
           if (cs==0)
               gamma = sign(-alfa(k))*pi/2;
           else
               gamma = atan(-alfa(k)/cs);
           end
           T(1,1) =  cos(gamma);
           T(1,k) = -sin(gamma);
           T(k,1) =  sin(gamma);
           T(k,k) =  cos(gamma);
           D3 = D3*T;
       end
   end

   if (nstoch<=3)
       IX4 = linspace(1,nstoch,nstoch)';
       D4 = zeros(nstoch);
       % % Eerste kolom
       D4(:,1) = -alfa;
       if (nstoch>1)
           % Tweede kolom
           ind = find(D4(:,1)==0);
           if ~isempty(ind)
               D4(ind(1),2) = 1;
           else
               ind = find(D4(:,1)~=0);
               if ~isempty(ind)
                   D4(ind(1),2) =  D4(ind(2),1);
                   D4(ind(2),2) = -D4(ind(1),1);
               else
                   werr = 600;
               end
           end
           for k = 3:1:nstoch
               D4(:,k) = cross(D4(:,k-2),D4(:,k-1));
           end
       end
   end


  if (sormoption==1)
       D = D1;
       IX = IX1;
  else
      if (sormoption==2)
           D = D2;
           IX = IX2;
      else
          if (sormoption==3)
               D = D3;
               IX = IX3;
          else
              if (sormoption==4)
                   D = D4;
                   IX = IX4;
              else
              end
          end
      end
  end

  nrot = 1;
  Dtmp = D;
  for irot = 1:1:nrot
       gamma = (irot-1)*2*pi/nrot;
       T = eye(nstoch);
       T(2,2) =  cos(gamma);
       T(2,3) = -sin(gamma);
       T(3,2) =  sin(gamma);
       T(3,3) =  cos(gamma);
       D = Dtmp*T;
%        disp(num2str(D(1,:)));
%        disp(num2str(D(2,:)));
%        disp(num2str(D(3,:)));
%        disp(num2str(D3(:,1)'*D3(:,2)));
%        disp(num2str(D3(:,1)'*D3(:,3)));
%        disp(num2str(D3(:,2)'*D3(:,3)));
%        disp(num2str(D3(:,1)'*D3(:,1)));
%        disp(num2str(D3(:,2)'*D3(:,2)));
%        disp(num2str(D3(:,3)'*D3(:,3)));
       dV = 0.3;%0.3;
       v = zeros(nstoch,1);
       Gk = zeros(nstoch-1,nstoch-1);
       %% default
       % Krommingen bepalen
       v(IX(1))=beta;
       for i1 = 2:1:nstoch
          i = IX(i1);
          for k1 = 2:1:nstoch
              k = IX(k1);
                    v(i) = v(i) - dV / 2.0D0;
                    v(k) = v(k) - dV / 2.0D0;
                    u = D*v;
                    [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
                Zm1 = func(x, arg);
                    v(k) = v(k) + dV;
                    u = D*v;
                    [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
                Zp1 = func(x, arg);
             Hm = (Zp1 - Zm1) / dV;
                    v(i) = v(i) + dV;
                    v(k) = v(k) - dV;
                    u = D*v;
                   [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
                Zm2 = func(x, arg);
                    v(k) = v(k) + dV;
                    u = D*v;
                    [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);
                Zp2 = func(x, arg);
             Hp = (Zp2 - Zm2) / dV;
             Gk(i1-1,k1-1) = -(Hp - Hm) / dV / sigZ;
%              Gk(k1-1,i1-1) = Gk(i1-1,k1-1);
                    v(k) = 0;
                    v(i) = 0;
          end
       end
       % Hoofdkrommingen bepalen
       %    CALL Jacoba (Gk, nstoch-1, NT%mstoch, Kr, GkV, werr)
      Kr = eig(Gk);
   disp(num2str(Gk(1,:)));
   disp(num2str(Gk(2,:)));


      disp([' gamma : ' num2str(gamma') ' krommingen : ' num2str(Kr') ...
          ' factor     : ' num2str(prod((1. - beta * Kr).^(-.5)))]);
      disp([' factoren   : ' num2str(((1. - beta * Kr).^(-.5))')]);
      disp([' factor     : ' num2str(prod((1. - beta * Kr).^(-.5)))]);
      disp([' beta Form  : ' num2str(beta)]);
      disp([' Pf   Form  : ' num2str(Pf)]);
  end


% Kans corrigeren volgens SORM
   for i = 1:1:nstoch - 1
      if (beta * Kr(i) > 1.)
         Pf = Pf * (1.013 + .1556 + 1.2896);
      else
         if (beta * Kr(i) > 0.6)
             Pf = Pf * (1.013 + .1556 * beta * Kr(i) + 1.2896 * (beta * Kr(i))^2); 
         else
             Pf = Pf * (1. - beta * Kr(i))^(-.5);
         end
      end
   end
%    beta = XfromQ(Pf);
    beta = norminv(1-Pf,0,1);

% Consistente set u, x en z waarden bepalen bij designpunt
   u = -alfa * beta;
   [u, x] = BepXfU(u, parameters, nstoch, trmH, arg);

results = FormRes;
parametersOut = FormPar;

% Resultaten
results.Pf = Pf;
results.beta = beta;
for i1=1:1:nstoch
   parametersOut(i1,1).u = u(i1);
   parametersOut(i1,1).x = x(i1);
end
