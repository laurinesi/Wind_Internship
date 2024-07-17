%% Explanation

%% Input
extremes = 1;
rho_steel = 1.0/100; 
fck = 40*10^6;
fyk = 500*10^6;
gamma_c = 1.5;
gamma_s = 1.15;

Ed = Md_base;    
if extremes == 2
    EdQ = Qd_base;
    
end

%% Calculations
b_square = (Ed/(rho_steel*fyk/gamma_s)*1/(1-0.5*rho_steel*(fyk/gamma_s)/(fck/gamma_c)))^(1/3); % Required width for massive concrete beam

fun = @(b) rho_steel*(b^2-(b-2*0.5)^2)*(fyk/gamma_s)*(b-0.39*((rho_steel*(b^2-(b-2*0.5)^2)*fyk/gamma_s)/(0.75*fck/gamma_c)/b)) - Ed;
b0=10;
b = fzero(fun,b0);
xu = ((rho_steel*(b^2-(b-2*0.5)^2)*fyk/gamma_s)/(0.75*fck/gamma_c)/b);
if extremes == 2
Asw_s = EdQ/((b-0.39*xu)*2.5*fyk/gamma_s);
end

% %% Check ultimate design value after all calculations. 
% b_core = 8;
% x = design_values_G;
% Rd = x(1)*x(2)*(b_core^2-(b_core-2*x(5))^2)*x(3)*(b_core-0.39*((x(2)*(b_core^2-(b_core-2*x(5))^2)*x(3))/(0.75*x(4)*b_core)));
