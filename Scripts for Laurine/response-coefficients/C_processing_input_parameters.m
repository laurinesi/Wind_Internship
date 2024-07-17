%% values resulting fom input parameters  
    
    % averaged reference wind speed in wind tunnel
        Pref=[];
        J = c_i(:,tapsRef);
        Pref = sum(sum(J)) / nnz(J);
        v_wt = sqrt(Pref/(0.5*rho_air)); % [m/s] resulting average wind speed at top building in wind tunnel 
    % averaged reference wind speed in full scale 
        v_fs = 0.19*(z0/z0_II)^0.07 * log(h_fs/z0)*v_pot; % [m/s] wind speed at height 120 [m] and terrain roughness 0.8 [m] (ASSUMPTION)
    % scaling values
        lambda_g = h_wt / h_fs; % [-] geometric scale 
        lambda_v = v_wt / v_fs; % [-] wind speed scale  
        lambda_t = 1/lambda_v * lambda_g; % [-] time scale
        lambda_f = 1/lambda_t; % [-] frequency scale
    % Sampling time
        t_sample=1/f_sampling;
        
%% minimum possible reference area

        f_max_wt = f_sampling/2; % [1/s] maximum possible frequency that is of importance of the load effect should be at maximum half of the sampling frequency [CUR, aanvulling]
        t_min_wt = 1/f_max_wt; % [s] the minimum gust length value in the wt scale
        t_min_fs = t_min_wt/lambda_t; % [s] the minimum gust length value in the fs scale
        L_min_ref_fs = t_min_fs*v_fs/1.5; % [m] the minimum reference length value
        A_min_ref_fs = (L_min_ref_fs/sqrt(2))^2; % [m^2] the minimum reference area      
 

%% required sampling frequency
    
        L_ref_fs = tap_area_dim(18,2); % [m] reference length value 
        t_gust_fs = 1.5*L_ref_fs/v_fs; % [s] duration of gust which is of importance for the load effect [Lawson]
        f_gust_fs = 1/t_gust_fs; % [Hz]
        f_wt_required = 2*lambda_f*f_gust_fs; % [Hz] minimum required sampling frequency corresponding to A_ref

        
%% required smoothing characteristics

        t_gust_wt = t_gust_fs*lambda_t;
        N_gust_wt = round(t_gust_wt*f_sampling); % number of samples that need to be averaged, round to bottom because of chosen A_ref
        
        
%% required duration of partial registrations

        T_wt = T_fs * lambda_t; % determination of the duration of partial registration in wt scale
        N_partial_wt = round(T_wt * f_sampling); % number of samples corresponding to one partial registration of time T_wt / T_fs
        
        