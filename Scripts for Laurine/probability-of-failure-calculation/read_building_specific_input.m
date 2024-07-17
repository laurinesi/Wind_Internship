% in deze file worden de juiste files geladen. Je moet wel zelf aangeven
% welke sections je op overdruk wil laten testen en welke op onderdruk. Je
% moet ook aangeven zelf hoe dat .. ontworpen dient te worden. % 1=max 0=min

% roughness factor
Vcr = 0.15; 
ce_mean_to_specified = 0.8; 
% structural strength
VR = 0.1 ; % (staat vast)


if extremes == 1
[R_mean_M, R_std_M] = parameters_LN(20.3462,0.0998);
%[R_mean_M, R_std_M] = parameters_LN(20.2713,0.05);

C_directory_G = '..\response-coefficients\fitted-dists';
C_filenames_G = 'C_Refbuilding_mom_Fx_G_MLE_20.xpt';   
C_directory_GEV = '..\response-coefficients\fitted-dists';
C_filenames_GEV = 'C_Refbuilding_mom_Fx_GEV_MLE_20.xpt';

crd = 1.1320; % Roughness factor at h_ref = 0.66h

elseif extremes == 2
[R_mean_Q, R_std_Q] = parameters_LN(16.1704,0.0998); 
%[R_mean_Q, R_std_Q] = parameters_LN(16.0925,0.05); 

C_directory_G = '..\response-coefficients\fitted-dists';
C_filenames_G = 'C_Refbuilding_for_Fx_G_MLE_20.xpt';   
C_directory_GEV = '..\response-coefficients\fitted-dists';
C_filenames_GEV = 'C_Refbuilding_for_Fx_GEV_MLE_20.xpt';

crd = 1.0655; % Roughness factor at h_ref = 0.5h
end


% roughness factor

crd_squared = crd^2;
crmean_squared = ce_mean_to_specified * crd_squared;
crstd_squared = Vcr*crmean_squared ;  %;






                