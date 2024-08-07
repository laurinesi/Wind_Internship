% Compare used dists (Schiphol) with the new dists of Alex (Cabauw)
% for wind speeds

%% settings

C_dist = 1; % [1 2 3] ~ [GEV GP GW]
F_dist = 2; % [1 2] ~ [G GEV]
S_dist = 3; % [1 2 3] ~ [GEV GP GW]
noise = 2; % [1 2] ~ [Gaussian Uniform]

%% GEV Cabauw
if C_dist == 1
    
    % loading struct df for GEV
    addpath('..\wind-speeds\Alex-dists\')
    load('GEV_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GEV';

    % initialization
    n = size(proba,1);
    Fx = zeros(n,1);

    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end
    
    A = [original, Fx];
    df.Fx = Fx;
    
    % saving data as a txt file
    writematrix(A,'..\wind-speeds\Alex-dists\V_Cabauw_maxima_tail_Fx_GEV','Delimiter','\t','FileType','text') 

end

%% GP Cabauw
if C_dist == 2

    % loading struct df for GP
    addpath('..\wind-speeds\Alex-dists\')
    load('GP_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GP';

    % initialization
    n = size(proba,1);
    Fx = zeros(n,1);

    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end
    
    A = [original, Fx];
    df.Fx = Fx;
    
    % saving data as a txt file
    writematrix(A,'..\wind-speeds\Alex-dists\V_Cabauw_maxima_tail_Fx_GP','Delimiter','\t','FileType','text') 

end

%% GW Cabauw
if C_dist == 3

    % loading struct df for GW
    addpath('..\wind-speeds\Alex-dists\')
    load('GW_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GW';

    % initialization
    n = size(proba,1);
    Fx = zeros(n,1);

    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end
    
    A = [original, Fx];
    df.Fx = Fx;
    
    % saving data as a txt file
    writematrix(A,'..\wind-speeds\Alex-dists\V_Cabauw_maxima_tail_Fx_GW','Delimiter','\t','FileType','text') 
    
end

%% G Schiphol, yearly maxima of hourly mean wind speeds (FH)
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

%% GEV Schiphol, yearly maxima of hourly mean wind speeds (FH)
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

%% Schiphol, 10min mean wind speeds

addpath('..\wind-speeds\fitted-dists\')
data_Schiphol = readmatrix("dataSchiphol.csv"); % measurements
model_Schiphol = readmatrix("RACMO_Schiphol.csv"); % model

% Initialisation
FFGaussNoise = zeros(height(data_Schiphol),1);
FFUnifNoise = zeros(height(data_Schiphol),1);

data(:,1) = data_Schiphol(:,2); % Year (starting the new year on july 1)
data(:,2) = data_Schiphol(:,3); % FF : Mean wind speed (m/s) during the 10-minute period preceding the time of observation at 10m (measurement KNMI)
if noise == 1
    data(:,3) = FFGaussNoise(:,1); % FF with gaussian noise
    data_Schiphol = array2table(data, "VariableNames",{'Year', 'FF', 'FFGaussNoise'});
elseif noise == 2
    data(:,3) = FFUnifNoise(:,1); % FF with uniform noise
    data_Schiphol = array2table(data, "VariableNames",{'Year', 'FF', 'FFUnifNoise'});
end

model(:,1) = model_Schiphol(:,2); % Year
model(:,2) = model_Schiphol(:,3); % F010 : Wind Speed at 10m Height
model(:,3) = model_Schiphol(:,14); % Maximum 10m Wind Speed including Gusts

model_Schiphol = array2table(model, "VariableNames",{'Year', 'F010', 'wgmax'});

% saving data as a txt file
writetable(model_Schiphol,'..\wind-speeds\fitted-dists\model_Schiphol','Delimiter','\t','FileType','text')

%% Gaussain noise
if noise == 1
     
    % Calculate the standard deviation of the FF column and scale it by 0.1
        % We use the function std so no MLE ?
    sigma = std(data_Schiphol.FF) * 0.1;
    
    % Add Gaussian noise
    FFGaussNoise = data_Schiphol.FF + randn(length(data_Schiphol.FF), 1) * sigma;
    
    % Set any negative values in the FFGaussNoise column to zero
    FFGaussNoise(FFGaussNoise < 0) = 0;
    
    % Store the modified data back in the table
    data_Schiphol.FFGaussNoise = FFGaussNoise;
    
    % Plot the CDFs
    figure;
    plot(sort(data_Schiphol.FFGaussNoise), linspace(0, 1, length(data_Schiphol.FFGaussNoise)), 'r');
    hold on;
    plot(sort(data_Schiphol.FF), linspace(0, 1, length(data_Schiphol.FF)), 'k');
    xlabel('Wind speeds (m/s)');
    ylabel('F_X(x)');
    title('ECDFs (Schiphol measurements)');
    legend('FFGaussNoise', 'FF');
    grid minor
    legend Box off
    
    % % plot the histogram of Schiphol measurements
    % figure
    % histogram(data_Schiphol.FF,150)
    % xlabel('Wind speeds (m/s)')
    % ylabel('Frequency')
    % title('Histogram of winds speeds data at Schiphol')

    
    % saving data as a txt file
    writetable(data_Schiphol,'..\wind-speeds\fitted-dists\data_Schiphol','Delimiter','\t','FileType','text') 

end

%% Uniform noise
if noise == 2

    % Compute the maximum difference between sorted FF elements
    I = max(diff(sort(data_Schiphol.FF)));
    
    % Add uniform noise
    FFUnifNoise = data_Schiphol.FF + (rand(length(data_Schiphol.FF), 1) - 0.5) * I;
    
    % Set any negative values in the FFUnifNoise column to zero
    FFUnifNoise(FFUnifNoise < 0) = 0;
    
    % Store the modified data back in the table
    data_Schiphol.FFUnifNoise = FFUnifNoise;
    
    % Plot the CDFs
    figure;
    plot(sort(data_Schiphol.FFUnifNoise), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'r');
    hold on;
    plot(sort(data_Schiphol.FF), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'k');
    xlabel('Wind speeds (m/s)');
    ylabel('F_X(x)');
    title('ECDFs (Schiphol measurements)');
    legend('FFUnifNoise', 'FF');
    grid minor
    legend Box off

    % saving data as a txt file
    writetable(data_Schiphol,'..\wind-speeds\fitted-dists\data_Schiphol','Delimiter','\t','FileType','text') 

end

%% GEV Schiphol 10min
if S_dist == 1

    % loading struct GEV_est
    addpath('..\wind-speeds\Alex-dists\')
    load('GEV_Schiphol_est.mat');
    
    probaS = GEV_est.proba;
    originalS = GEV_est.original;
    S10distType = 'GEV';
    
    % initialization
    n = size(probaS,1);
    FxS = zeros(n,1);
    
    for i = 1:n
        FxS(i,1) = 1 - probaS(i,1);
    end
    
    B = [originalS, FxS];
    df.Fx = FxS;
    
    % saving data as a txt file
    writematrix(B,'..\wind-speeds\Alex-dists\V_Schiphol_maxima_tail_Fx_GEV','Delimiter','\t','FileType','text') 

end

%% GP Schiphol 10min
if S_dist == 2

    % loading struct GP_est
    addpath('..\wind-speeds\Alex-dists\')
    load('GP_Schiphol_est.mat');

    probaS = GP_est.proba;
    originalS = GP_est.original;
    S10distType = 'GP';
    
    % initialization
    n = size(probaS,1);
    FxS = zeros(n,1);
    
    for i = 1:n
        FxS(i,1) = 1 - probaS(i,1);
    end
    
    B = [originalS, FxS];
    df.Fx = FxS;
    
    % saving data as a txt file
    writematrix(B,'..\wind-speeds\Alex-dists\V_Schiphol_maxima_tail_Fx_GP','Delimiter','\t','FileType','text') 

end

%% GW Schiphol 10min
if S_dist == 3

    % loading struct GW_est
    addpath('..\wind-speeds\Alex-dists\')
    load('GW_Schiphol_est.mat');

    probaS = GW_est.proba;
    originalS = GW_est.original;
    S10distType = 'GW';
    
    % initialization
    n = size(probaS,1);
    FxS = zeros(n,1);
    
    for i = 1:n
        FxS(i,1) = 1 - probaS(i,1);
    end
    
    B = [originalS, FxS];
    df.Fx = FxS;
    
    % saving data as a txt file
    writematrix(B,'..\wind-speeds\Alex-dists\V_Schiphol_maxima_tail_Fx_GW','Delimiter','\t','FileType','text') 

end

%% plotting cdf

C_name = append('Cabauw', ' ', CdistType);
S_name = append('Schiphol 1h', ' ', SdistType);
S10_name = append('Schiphol 10min', ' ', S10distType);

figure
plot(original, Fx,'-','color','k') % Cabauw
hold on
plot(A_S_crop(:,1), A_S_crop(:,2),'--','color','k')  % Schiphol (hourly mean)
hold on
plot(originalS, FxS,'-.','color','k') % Schiphol (10min mean)
xlabel('Wind Speed (m/s)')
ylabel('F_X')                      
grid minor
legend(C_name, S_name, S10_name)
legend box off
title('CDFs of wind speeds (tail)')

