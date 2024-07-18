
b_core = 10;                                                                % [m] Width of concrete core true case study building similar to a reference building in Riga, Latvia
t_core = 0.5;                                                               % [m] Thickness of concrete core case study building
h_floor = 0.5;
perc_floor = 1.0;
N_floors = h/3;
m_guess = 1.2*((b_core^2-(b_core-t_core)^2)*h*rho_conc+(b^2-b_core^2)*h_floor*rho_conc*perc_floor*N_floors)/h;                       % [kg] Estimation mass of true case study building
