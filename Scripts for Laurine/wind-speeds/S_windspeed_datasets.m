% read datasets from measurements and model (RACMO) at Schiphol
% add noise (Gaussian or Uniform) to smooth the data as Alex did
% plot ECDF of original and noisy data
% save new data in txt files : model_Schiphol, data_Schiphol
% - Column 1 is the year
% - Column 2 is the wind speed (m/s)
% - Column 3 is the gusts / smoothed dataset of wind speed

addpath('..\wind-speeds\datasets\')
data_Schiphol = readmatrix("dataSchiphol.csv"); % measurements
model_Schiphol = readmatrix("RACMO_Schiphol.csv"); % model

% settings
noise = 2; % [1 2] ~ [Gaussian Uniform]

% Initialization
FFGaussNoise = zeros(height(data_Schiphol),1);
FFUnifNoise = zeros(height(data_Schiphol),1);
data = [];

% Definition of tables
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
model(:,2) = model_Schiphol(:,3); % F010 : Wind Speed at 10m Height (3-hr point values)
model(:,3) = model_Schiphol(:,14); % Maximum 10m Wind Speed including Gusts (3-hr maximum values)

model_Schiphol = array2table(model, "VariableNames",{'Year', 'F010', 'wgmax'});

% saving data as a txt file
writetable(model_Schiphol,'model_Schiphol','Delimiter','\t','FileType','text')

%% Gaussain noise
if noise == 1
     
    % Calculate the standard deviation of the FF column and scale it by 0.1
    sigma = std(data_Schiphol.FF) * 0.1;
    
    % Add Gaussian noise
    FFGaussNoise = data_Schiphol.FF + randn(length(data_Schiphol.FF), 1) * sigma;
    
    % Set any negative values in the FFGaussNoise column to zero
    FFGaussNoise(FFGaussNoise < 0) = 0;
    
    % Store the modified data back in the table
    data_Schiphol.FFGaussNoise = FFGaussNoise;
    
    % Plot the CDFs for Gaussian noise
    figure;
    plot(sort(data_Schiphol.FFGaussNoise), linspace(0, 1, length(data_Schiphol.FFGaussNoise)), 'r');
    hold on;
    plot(sort(data_Schiphol.FF), linspace(0, 1, length(data_Schiphol.FF)), 'k');
    xlabel('Wind speeds (m/s)');
    ylabel('F_X(x)');
    title('Comparison between original and noisy data');
    legend('FFGaussNoise', 'FF');
    legend Box off
    xlim([0, 20]); % Limiting the x-axis to 20 m/s
    
    % saving data as a txt file
    writetable(data_Schiphol,'data_Schiphol','Delimiter','\t','FileType','text') 

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
    
    % Plot the CDFs for Uniform noise
    figure;
    plot(sort(data_Schiphol.FFUnifNoise), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'r');
    hold on;
    plot(sort(data_Schiphol.FF), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'k');
    xlabel('Wind speeds (m/s)');
    ylabel('F_X(x)');
    title('Comparison between original and noisy data');
    legend('FFUnifNoise', 'FF');
    legend Box off
    xlim([0, 20]); % Limiting the x-axis to 20 m/s

    % saving data as a txt file
    writetable(data_Schiphol,'data_Schiphol','Delimiter','\t','FileType','text') 

end

%% Compare measurements and model

% Plot the empirical CDFs comparing measurements and model
figure;
if noise == 1
    plot(sort(data_Schiphol.FFGaussNoise), linspace(0, 1, length(data_Schiphol.FFGaussNoise)), 'k');
elseif noise == 2
    plot(sort(data_Schiphol.FFUnifNoise), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'k');
end
hold on
plot(sort(model_Schiphol.F010), linspace(0, 1, length(model_Schiphol.F010)), '--r');
hold on
plot(sort(model_Schiphol.wgmax), linspace(0, 1, length(model_Schiphol.wgmax)), '--b');
xlabel('Wind speeds (m/s)');
ylabel('F_X(x)');
title('Empirical CDFs of each whole dataset at Schiphol');
legend('Measurements', 'RACMO', 'including gusts',Location='east');
legend Box off
xlim([0, 27]); % Limiting the x-axis to 27 m/s


% Histogram of Schiphol measurements
% figure
% histogram(data_Schiphol.FF,150)
% xlabel('Wind speeds (m/s)')
% ylabel('Frequency')
% title('Histogram of winds speeds data at Schiphol')

