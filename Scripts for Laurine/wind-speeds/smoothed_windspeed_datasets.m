% read datasets from measurements and model (RACMO)
% add noise (Gaussian or Uniform) to smooth the data
% plot ECDF of original and noisy data
% save new data in txt files : model_Schiphol, data_Schiphol

addpath('..\wind-speeds\')
data_Schiphol = readmatrix("dataSchiphol.csv"); % measurements
model_Schiphol = readmatrix("RACMO_Schiphol.csv"); % model

%% settings
noise = 2; % [1 2] ~ [Gaussian Uniform]

% Initialization
FFGaussNoise = zeros(height(data_Schiphol),1);
FFUnifNoise = zeros(height(data_Schiphol),1);

% Definition of tables
data(:,1) = data_Schiphol(:,2); % Year (starting the new year on july 1)
data(:,2) = data_Schiphol(:,3); % FF : Hourly mean wind speed (m/s) during the 10-minute period preceding the time of observation at 10m (measurement KNMI)
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
    
    % Plot the empirical CDFs
    figure;
    plot(sort(data_Schiphol.FFGaussNoise), linspace(0, 1, length(data_Schiphol.FFGaussNoise)), 'r');
    hold on;
    plot(sort(data_Schiphol.FF), linspace(0, 1, length(data_Schiphol.FF)), 'k');
    xlabel('Wind speeds (m/s)');
    ylabel('F_X(x)');
    title('ECDFs (Schiphol measurements)');
    legend('FFGaussNoise', 'FF');
    %grid minor
    legend Box off
    
    % % plot the histogram of Schiphol measurements
    % figure
    % histogram(data_Schiphol.FF,150)
    % xlabel('Wind speeds (m/s)')
    % ylabel('Frequency')
    % title('Histogram of winds speeds data at Schiphol')

    
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
    
    % Plot the empirical CDFs
    figure;
    plot(sort(data_Schiphol.FFUnifNoise), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'r');
    hold on;
    plot(sort(data_Schiphol.FF), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'k');
    xlabel('Wind speeds (m/s)');
    ylabel('F_X(x)');
    title('ECDFs (Schiphol measurements)');
    legend('FFUnifNoise', 'FF');
    %grid minor
    legend Box off

    % saving data as a txt file
    writetable(data_Schiphol,'data_Schiphol','Delimiter','\t','FileType','text') 

end

%% Compare measurements and model

% Plot the empirical CDFs
figure;
if noise == 1
    plot(sort(data_Schiphol.FFGaussNoise), linspace(0, 1, length(data_Schiphol.FFGaussNoise)), 'k');
elseif noise == 2
    plot(sort(data_Schiphol.FFUnifNoise), linspace(0, 1, length(data_Schiphol.FFUnifNoise)), 'k');
end
hold on
plot(sort(model_Schiphol.F010), linspace(0, 1, length(model_Schiphol.F010)), 'b');
hold on
plot(sort(model_Schiphol.wgmax), linspace(0, 1, length(model_Schiphol.wgmax)), 'Color',  "#0072BD");
xlabel('Wind speeds (m/s)');
ylabel('F_X(x)');
title('Empirical CDFs of each whole dataset at Schiphol');
legend('Measurements', 'RACMO', 'including gusts');
legend Box off

