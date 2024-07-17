function [u0, x] = BepXfU (u0, parameters, nstoch, trmH, arg)
% $ Declare
%    SUBROUTINE BepXfU (u, x, gemX, sigX, parameters(i).type_v, ipoint, nstoch)
%
% Subroutine die afhankelijk van de verdeling X uit U berekent
%
% Parameters
%   u      I/O Vector van U-waarden
%   x       O  Waarde stochast
%   parameters  I  Verdelingstype en data
%   nstoch  I  Aantal stochasten
%   trmH    I   transformation matrix for correlated parameters
%   arg     I  free form: user parameters voor z-functie en BepXfU

%
% TNO Bouw Sept 2006 SNH
% ---------------------------------------------------------------------- 
%
if (nargin==3)
    trmH = eye(nstoch);
end
% Omrekenen van ongecorreleerd naar gecorreleerde variabelen
u = trmH * u0;
%
x = zeros(nstoch,1);
% Voor alle stochasten (bekende verdelingen)
for i = 1:1:nstoch
% Normale verdeling 
    if strcmp(parameters(i).type_v(1:3),'NOR')
      x(i) = parameters(i).gemX + u(i) * parameters(i).sigX;
    else
        
    if strcmp(parameters(i).type_v(1:3),'LOG')
        mean = parameters(i).gemX;
        stdev = parameters(i).sigX;
        shift = parameters(i).epsX;
% 		double logmean = 0.5* Math.log( Math.pow(mean-shift, 4)/ (Math.pow (mean-shift,2) +Math.pow(stdev, 2)));
% 		double logsd = Math.sqrt ( Math.log( 1 +  Math.pow(stdev,2)/Math.pow(mean-shift, 2)));
% 		if (Double.isNaN(logmean) || Double.isNaN(logsd)) return logmean + logsd;
% 		if(p < 0 || p > 1 || logsd <= 0) {
% 			throw new java.lang.ArithmeticException("Math Error: DOMAIN");
% 		//            return Double.NaN;
% 		}
% 		if (p == 0) return shift;
% 		if (p == 1) return Double.POSITIVE_INFINITY;
% 		return shift + java.lang.Math.exp(normal.quantile(p, logmean, logsd));
		logmean = 0.5 * log( ((mean-shift)^4) / ((mean-shift)^2 + stdev^2));
		logsd = sqrt( log( 1 +  (stdev^2)/((mean-shift)^2)));
        p = normcdf(u(i),0,1);
		x(i) = shift + exp(norminv(p, logmean, logsd));
        x(i) = shift + exp(logmean + u(i) * logsd);
    else
      
        if strcmp(parameters(i).type_v(4:8),'TRUNC')
% Trunc verdeling
            x(i) = parameters(i).gemX + u(i) * parameters(i).sigX;
            if (x(i)<0.)
                x(i)= 0;
                u(i) = (x(i) - parameters(i).gemX)/parameters(i).sigX;
            end
        else
            if strcmp(parameters(i).type_v(1:3),'EXP')
% Exponentionele verdeling
%                 [p q ] = QfromX (u(i))
%                 if (q < 1e-300); q = 1e-300; end;
%                 x(i) = parameters(i).sigX - log(q) * parameters(i).gemX;
                p = normcdf(u(i),0,1);
                x(i) = expinv(p,parameters(i).gemX) + parameters(i).sigX;
            else
                if strcmp(parameters(i).type_v(1:3),'UNI')
% Uniforme verdeling
                     [p q ] = QfromX (u(i));
%                      p = normcdf(u(i),0,1);
                     x(i) = parameters(i).gemX + (p - 0.5D0) * parameters(i).sigX;
                else
                    if strcmp(parameters(i).type_v(1:3),'TAB')
% tabel
                         p = normcdf(u(i),0,1);
                         x(i) = interp1(parameters(i).data(:,2),parameters(i).data(:,1),p);
                    else
                        if strcmp(parameters(i).type_v(1:3),'DET')
% Deterministisch
                            x(i) = parameters(i).gemX;
                        else
                            if strcmp(parameters(i).type_v(1:3),'GUM')
% Gumbel
                               [p q] = QfromX(u(i));
                                 x(i) = parameters(i).gemX - log(-log(p)) / parameters(i).sigX  ;
%                                if (u(i)>5.4)
%                                    x(i) = parameters(i).gemX - log(q) / parameters(i).sigX  ;
%                                else
%                                    x(i) = parameters(i).gemX - log(-log(p)) / parameters(i).sigX  ;
%                                end
                            else
                                    if strcmp(parameters(i).type_v(1:3),'UXT')
% UXT
                                         x(i) = interp1(parameters(i).data(:,1),parameters(i).data(:,2),u(i));
                                    else
                                    if strcmp(parameters(i).type_v(1:3),'MIX')
% MIX
%                                        [p q] = QfromX(u(i));
                                       p = normcdf(u(i),0,1);
                                       omega = parameters(i).epsX;
                                       w = parameters(i).gemX;
                                       k = parameters(i).sigX;
                                       x(i) = omega - (omega-w)*(-log(p))^(1/k);
                                    else
                                    if strcmp(parameters(i).type_v(1:3),'QUA')
% QUA
                                       [p q] = QfromX(u(i));
                                       x(i) = 1 - (1-p)^2; 
                                    else
                                        if strcmp(parameters(i).type_v(1:3),'LIN')
% LIN
                                           [p q] = QfromX(u(i));
                                           x(i) = p; 
                                        else
                                            if strcmp(parameters(i).type_v(1:3),'BET')
% BETA
                                               [p q] = QfromX(u(i));
                                               x(i) = betainv(p,parameters(i).gemX,parameters(i).sigX); 
                                            else
                                                if strcmp(parameters(i).type_v(1:3),'WEI')
% WEIBULL
                                                   [p q] = QfromX(u(i));
                                                   x(i) = wblinv(p,parameters(i).gemX,parameters(i).sigX); 
                                                else
                                                    disp(['BepXfU: distr type not found ' parameters(i).type_v(i,1:3)]);
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                            end
                        end
                        end
                    end
                end
            end
        end
    end
end
end
