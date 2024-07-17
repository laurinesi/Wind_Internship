

%% determination of design strength according to Eurocode
% INPUT PARAMETERS
    % input parameters differing for each tap, real parameters
    cp_ww= 0.8;
    cp_lw= -0.65; 
    cf = 2.4;                                                               % Force coefficient for aspect ratio case study building

    % general input parameters, the same for all taps
    rho = 1.25;
    rho_conc = 2400;                                                        % [kg/m3] Density of concrete
    cdir = 1;
    cseason= 1;
    cprob = 1;
    vb0 = 27.0;
    kI = 1;
    co = 1;
    z0 = 0.5;
    z0II = 0.05; 
   
    gamma_wind = 1.5;
    % building geometry, the same for all taps
    h = 120; %[m]
    b = 30; %[m]
    n_1x = 46/h; 
    % information, the same for all taps
    location = 'Schiphol';
    wind_area = 'II'; 
% CALCULATIONS IN EUROCODE
  run('mu_guess');
    % standard values for all directions
    vb=cdir*cseason*cprob*vb0;
    kr = 0.19*(z0/z0II)^0.07;
    cr_ze=kr*log((0.6*h)/z0);                                       %Terrain roughness factor at href

    % Calculations for background factor
    vm_ze=cr_ze*co*vb;                                                      % Mean wind velocity at href
                Lt=300;                                                     % Reference length scale (EN 1991-1-4 B.1)
                zt=200;                                                     % Reference height (EN 1991-1-4 B.1)
                alpha=0.67+0.05*log(z0);                                    % (EN 1991-1-4 B.1)
            L=Lt*(0.6*h/zt)^alpha;                                          % [-] Turbulence length scale for href=0.6h for vertical structure like buildings
    B2=1/(1+3/2*sqrt((b/L)^2+(h/L)^2+((b*h)/L^2)^2));                       %[-] Background factor
    Iv_ze=(kr*vb*kI)/vm_ze;                                                 % Turbulence intensity at href

    % Calculations for resonant response factor         
            S_l_ze = 6.8*(n_1x*L/vm_ze)/((1+10.2*(n_1x*L/vm_ze))^(5/3));    % Non-dimensional power spectral density function of wind at reference height (EN1991-1-4 B.2)
                Gy = 1/2;                                                   % Constant depending on mode shape in width direction
                Gz = 3/8;                                                   % Constant depending on mode shape in height direction: linear mode shape assumption
                %Gz = 5/18;                                                  % Parabolic mode shape assumption
                phi_y = 11.5*b*n_1x/vm_ze;                                  % Non-dimensional frequency parameter
                phi_z = 11.5*h*n_1x/vm_ze;                                  % Non-dimensional frequency parameter
            K_s_ze = 1/(1+sqrt((Gy*phi_y)^2+(Gz*phi_z)^2+((2/pi())*Gy*phi_y*Gz*phi_z)^2)); % Size reduction function (EN1991-1-4 C.3)
            delta_s = 0.10;                                                 % Logaritmic decrement for structural damping for reinforced concrete structures               
                mu = m_guess;                                    % Mass per unit area
                % Calculatios equivalent mass
                fun = @(z) mu.*((z./h).^1.5).^2;
                fun2 = @(z) ((z./h).^1.5).^2;
                q1 = integral(fun,0,h);
                q2 = integral(fun2,0,h);
                mu_e = q1/q2;                                               % Equivalent mass per unit length given the mode shape     
            delta_a = cf*rho*vm_ze/(2*n_1x*mu_e);                           % Logarithmic decrement for aerodynamic damping
            delta = delta_s + delta_a;                                      % Total logarithmic decrement for damping
    R2 = (pi()^2/(2*delta))*S_l_ze*K_s_ze;                                  %[-] Resonant response factor method Annex C
% General factors
        nu = n_1x*sqrt(R2/(B2+R2));                                         % Upcrossing frequency 
    kp = sqrt(2*log(nu*600))+0.6/sqrt(2*log(nu*600));                       % Peak factor (EN1991-1-4 (B.4)
        cs=(1+7*Iv_ze*sqrt(B2))/(1+7*Iv_ze);                                    % [-] Size factor
        cd_dyn = (1+2*kp*Iv_ze*sqrt(B2+R2))/(1+7*Iv_ze*sqrt(B2));               % Dynamic factor for dynamic analysis
    cscd=cs*cd_dyn;                                                             % [-] Structural factor
    if cscd<0.85; 
        cscd=0.85; % [-] Minimum value for structural factor
    end

    % depending on reference location
    h_strip=5;                                                            %[m] to get a Aref=10 [m2] with b=30 [m] NECESSARY?
    Aref = b*h_strip;
    ref=h/h_strip;                                                          %amount of strips over height
    z=[];
    for j=1:ref 
    z(j+1)=j*h_strip;                                                             % Height vector for case study building
    end
    z(1) = 0;
    z_e = [];
    for i = 2: length(z)
    if z(i) < b
            z_e(i) = b;
    
    elseif z(i) > (h-b)
            z_e(i) = h;
    else 
            z_e(i) = z(i);
    end
    cr(i)=kr*log(z_e(i)/z0);
    vm(i)=cr(i)*co*vb;
    Iv(i)=(kr*vb*kI)/vm(i);
    G(i)=1+7*Iv(i);
    qp(i)=G(i)*1/2*rho*vm(i)^2;
    Fchar_ww(i)=qp(i)*cp_ww*cscd*Aref;
    Fchar_lw(i)=qp(i)*cp_lw*cscd*Aref;
    Fd_max(i)=gamma_wind*0.85*(Fchar_ww(i)-Fchar_lw(i));
    Md_max(i) = Fd_max(i)*((z(i)+z(i-1))/2);
    end
Qd_base = sum(Fd_max);
Md_base = sum(Md_max);
summary_input_parameters = transpose([z; z_e; cr; vm; Iv; qp; Fd_max; Md_max]);

% figure
% stairs(Fd_max,z);
% figure
% stairs(Md_max,z);

    





