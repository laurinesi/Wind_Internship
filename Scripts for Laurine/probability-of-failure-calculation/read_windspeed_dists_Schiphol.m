% load the xpt - tables. 

% statistical uncertainties for each wind directions
V_cov_G=0.05;%[0.04];
V_cov_GEV=0.10;%[0.09];

V_directory_G = '..\wind-speeds\fitted-dists';
V_filenames_G = 'V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt';     
           
V_directory_GEV = '..\wind-speeds\fitted-dists';
V_filenames_GEV = 'V_Schiphol_maxima_conf1_Fx_GEV_MLE_sec1.xpt';     
           
% probability of occurence wind speed different sections KLOPT DIT NOG WEL?
%load ('\\tsn.tno.nl\data\SV\sv-069776\Analysis\P_Failure calculation\stochast V\Fx 50-yearly extremes\configuration1\p_sec_conf1');           
