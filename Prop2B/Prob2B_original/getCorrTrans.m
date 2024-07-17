function [trmH] = getCorrTrans(correlations)    
%
% Subroutine voor het opgeven van de correlaties en het berekenen 
% van de inverse transformatie matrix H
%
% Parameters
%   correlations  I   Covariantie matrix
%   trmH          O   Transformation matrix
%
% TNO Bouw November 2004 SNH
% --------------------------------------------------------------------- 
%
% Eigenwaarden en vectoren uitrekenen
%    CALL Jacoba(NT%Vcov, maxnrc, NT%mstoch, D, V, werr)

maxnrc = size(correlations);
trmH = eye(maxnrc);
%
[V Dtmp ] = eig(correlations);
for i=1:1:maxnrc,D(i)=Dtmp(i,i);end;
%
% Eigenwaarden en vectoren sorteren
for i = 1:1:maxnrc
    kol = i;
    for j = i + 1:1:maxnrc
        if (abs(V(i,j)) > abs(V(i, kol))), kol = j;end;
    end
    if (kol > i)
        tmp    = D(i);
        D(i)   = D(kol);
        D(kol) = tmp;
        tmp      = V(:,i);
        V(:,i)   = V(:,kol);
        V(:,kol) = tmp;
    end  
    if (V(i,i) < 0)
      V(:,i) = -V(:,i);
    end
end
%
% Inverse transformatie matrix H berekenen
for i = 1:1:maxnrc
  if (D(i) < 0)
     disp('getCorrTrans; neagtive eigenvalue found');
     return;
  end
  if (D(i) > 0)
    tmp = sqrt(D(i));
  else
    tmp = 0.;
  end
  trmH(:,i) = V(:,i) * tmp;
end

% trmH(5,5) = 1;
% trmH(5,6) = 0;
% trmH(6,5) = 0.820;
% trmH(6,6) = 0.572;
end
%
