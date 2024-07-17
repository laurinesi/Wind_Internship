%% settings and preferences

% name of data
data_mom = 'data\c_m_24hr.mat';
data_for = 'data\c_f_24hr.mat';
data_stat = 'data\c_m_stat.mat';
% name of the building
    building = 'Ref. building';

% what times for partial registrations do you want to investigate?
    T_fs  = 20; %[10 20 3600];
    T_3600 = 3600;
% required fracile values
    P_cook = 0.78;
    P_code = 0.95;
    P_d = 0.999;
    P_fractiles = [P_cook P_code P_d];
% saving
    saving_directory = 'figures';
% atmospheric circumstances input parameters  
    rho_air = 1.225; % [kg/m^3] air density in the wind tunnel at time of measurement (check with temperature?)        
% geometry input parameters
    h_wt = 0.48; % [m] height wind tunnel
    h_fs = 120;  % [m] height full scale        
% wind speed input parameters
    v_pot = 27.0; % [m/s] potential wind speed at schiphol [windcursus] (corresponds to 50 year reference period, z0=0.03 [m] and z=10 [m])% hier eventueel de echte waarde pakken die ik zelf heb afgeleid?
    z0 = 0.8; % [m] terrain roughness 
    z0_II = 0.03; % [m] standard terrain roughness (terrain category II Eurocode ,moet nu 0.03 zijn toch?)   
% time related input parameters
    f_sampling = 400;   % [Hz] sampling frequency 
    data_skip = 250; % Datapoints that are skipped due to discontinuities at beginning of DIANA time model
% reference pressure 
    tapsRef = [87 88 176 177];
    
%% 
    

load('data\c_i.mat');
load(data_mom);
load(data_for);
load(data_stat);
c_m = c_m_24hr;
c_f = c_f_24hr;
c_m_stat = c_m_serie;
load('data\tap_area_dim.mat');
run('C_processing_input_parameters');

%% determine dependency of extremes (maxima and minima)
        
if f_wt_required < f_sampling % calculation only needs to be executed if f_wt_required < f_sampling
   X_t_mom = [];
   X_t_for = [];
   for T = T_fs  
        % simplify sample data
            X_t_mom = c_m(:);     
            X_t_for = c_f(:);   
            X_t_stat = c_m_stat(:);
        % required duration of partial registrations
            T_wt = T * lambda_t; % determination of the duration of partial registration in wt scale
            N_partial_wt = round(T_wt * f_sampling); % number of samples corresponding to one partial registration of time T_wt / T_fs        
        % normalise the data (divide by reference pressure)
            X_norm_mom = X_t_mom(data_skip+1:end) ;   
            X_norm_for = X_t_for(data_skip+1:end) ;   
            X_norm_stat = X_t_stat(data_skip+1:end) ;
        % smooth the data 
            X_norm_smooth_mom = X_norm_mom;
            X_norm_smooth_for = X_norm_for;
            X_norm_smooth_stat = X_norm_stat;
            N_total = size(X_norm_mom,1);
            N_total_stat = size(X_norm_stat,1);
        % number of partial registrations (blocks)
            N_blocks = floor(N_total / N_partial_wt); % maximum number of blocks that can be obtained from the dataset
            N_rows_max = N_blocks * N_partial_wt; % remove last rows
            N_blocks_stat = floor(N_total_stat / N_partial_wt); % maximum number of blocks that can be obtained from the dataset
            N_rows_max_stat = N_blocks_stat * N_partial_wt; % remove last rows
            X_norm_smooth_red_mom = X_norm_smooth_mom(1 : N_rows_max,:);
            X_norm_smooth_red_for = X_norm_smooth_for(1 : N_rows_max,:);
            X_norm_smooth_red_stat = X_norm_smooth_stat(1 : N_rows_max_stat,:);
            N_red = size(X_norm_smooth_red_mom,1);
        % break op data into partial registrations
            X_norm_smooth_red_blocks_mom = reshape(X_norm_smooth_red_mom, [N_partial_wt, N_blocks]);
            X_norm_smooth_red_blocks_for = reshape(X_norm_smooth_red_for, [N_partial_wt, N_blocks]);
            X_norm_smooth_red_blocks_stat = reshape(X_norm_smooth_red_stat, [N_partial_wt, N_blocks_stat]);
            X_mom = max(X_norm_smooth_red_blocks_mom);
            X_for = max(X_norm_smooth_red_blocks_for);
            X_stat = max(X_norm_smooth_red_blocks_stat);
        % save variables  
            Cm_max{T} = sort(X_mom);
            Cf_max{T} = sort(X_for);
            Cm_stat{T} = sort(X_stat);
            N{T} = size(Cm_max{T},2);
            N_stat{T} = size(Cm_stat{T},2);
        end
end

%% Clear variables

clearvars -except Cm_max Cf_max Cm_stat N N_stat T_fs T_3600

