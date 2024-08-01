% Compare used dists (Schiphol) with the new dists of Alex (Cabauw)

%% settings

A_dist = 1; % [1 2 3] ~ [GEV GP GW]
F_dist = 2; % [1 2] ~ [G GEV]

%% GEV Cabauw
if A_dist == 1
    
    % loading struct df for GEV
    addpath('..\wind-speeds\Alex-dists\')
    load('GEV_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GEV';

    % initialization
    n = size(proba,1);
    Fx = zeros(n,1);

    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end
    
    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    A = [original, Fx];
    df.Fx = Fx;
    
    % saving data as a txt file
    writematrix(A,'V_Cabauw_maxima_tail_Fx_GEV','Delimiter','\t','FileType','text') 

end

%% GP Cabauw
if A_dist == 2

    % loading struct df for GP
    addpath('..\wind-speeds\Alex-dists\')
    load('GP_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GP';

    % initialization
    n = size(proba,1);
    Fx = zeros(n,1);

    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end
    
    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    A = [original, Fx];
    df.Fx = Fx;
    
    % saving data as a txt file
    writematrix(A,'V_Cabauw_maxima_tail_Fx_GP','Delimiter','\t','FileType','text') 

end

%% GW Cabauw
if A_dist == 3

    % loading struct df for GW
    addpath('..\wind-speeds\Alex-dists\')
    load('GW_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GW';

    % initialization
    n = size(proba,1);
    Fx = zeros(n,1);

    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end
    
    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    A = [original, Fx];
    df.Fx = Fx;
    
    % saving data as a txt file
    writematrix(A,'V_Cabauw_maxima_tail_Fx_GW','Delimiter','\t','FileType','text') 
    
end

%% G Schiphol
if F_dist == 1

    addpath('..\wind-speeds\fitted-dists\');
    %load('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt');
    A_S = dlmread('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt','\t');

    % 
    SdistType = 'G';

    % Define the start and end indices
    startIdx = 8324;
    endIdx = 60142;
    
    % Calculate the number of rows to copy
    numRowsToCopy = endIdx - startIdx + 1;
    
    % Preallocate A_S_crop
    A_S_crop = zeros(numRowsToCopy, size(A_S, 2));
    
    % Copy the data from A_S to A_S_crop
    A_S_crop(:, :) = A_S(startIdx:endIdx, :);

end

%% GEV Schiphol
if F_dist == 2

    addpath('..\wind-speeds\fitted-dists\');
    %load('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt');
    A_S = dlmread('V_Schiphol_maxima_conf1_Fx_GEV_MLE_sec1.xpt','\t');
    %A_S = xptread('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt');
    
    %
    SdistType = 'GEV';

    % Define the start and end indices
    startIdx = 6346;
    endIdx = 60142;
    
    % Calculate the number of rows to copy
    numRowsToCopy = endIdx - startIdx + 1;
    
    % Preallocate A_S_crop
    A_S_crop = zeros(numRowsToCopy, size(A_S, 2));
    
    % Copy the data from A_S to A_S_crop
    A_S_crop(:, :) = A_S(startIdx:endIdx, :);

end

%% plotting cdf

C_dist = append('Cabauw', ' ', CdistType);
S_dist = append('Schiphol', ' ', SdistType);

figure
plot(original, Fx,'-','color','k') % Cabauw
hold on
plot(A_S_crop(:,1), A_S_crop(:,2),'--','color','k')  % Schiphol
xlabel('Wind Speed (m/s)')
ylabel('F_X')                      
grid minor
legend(C_dist, S_dist)
legend box off
title('CDFs of wind speeds (tail)')

