clear all
close all
clc

addpath('..\wind-speeds\tools\')
run("C_windspeed_datasets.m") % Cabauw measure & model (RACMO, KNW)

% working with the measurement data at Cabauw
dataset = measurements_Cabauw;

%% dataset where wind speed is > 5
wind_speed_columns = {'F010', 'F020', 'F040', 'F080', 'F140', 'F200'};
binary_matrices = zeros(height(dataset), numel(wind_speed_columns));
for i = 1:numel(wind_speed_columns)
    binary_matrices(:, i) = dataset.(wind_speed_columns{i}) > 5;
end
row_index = all(binary_matrices, 2);
filtered_dataset = dataset(row_index,:);

%% get 30% highest wind speeds
valid_idx = find(~isnan(dataset.F010));
valid_dataset = dataset(valid_idx, :);
reorder_dataset = sortrows(valid_dataset,"F010","descend");
new_idx = find(~isnan(reorder_dataset.F010));

index_30 = 30*length(reorder_dataset.Year)/100;
dataset_30 = reorder_dataset(1:index_30, :);

% separation time = 12h
dataset_30 = sortrows(dataset_30, "DateTime", "descend");
separation_hours = 12;
selected_rows = dataset_30(1, :);

for i = 2:height(dataset_30)
    current_time = dataset_30.DateTime(i);
    time_diffs = hours(current_time - selected_rows.DateTime);
    if all(abs(time_diffs) >= separation_hours)
        selected_rows = [selected_rows; dataset_30(i, :)];
    end
end
% sample fraction = 10%
selected_rows = sortrows(selected_rows, "F010", "descend");
top_10_percent_count = round(0.1 * height(selected_rows));
dataset_10 = selected_rows(1:top_10_percent_count, :);


%% get the yearly maximum wind speeds
% Initialization
sigma_analysis_z0 = table;

% get annual maximum
[max_values, max_datetime] = BM_select(dataset);
% visualization
figure
plot(max_values(:,1),max_values(:,2),'Marker','x','LineStyle','none')
xlabel('Year');
ylabel('Wind speeds (m/s)');
title('Yearly maximum wind speeds');

% find line of max_datetime in dataset
for i = 1:length(max_datetime)
    idx(i,:) = find(dataset.DateTime==max_datetime(i));
end
% store annual storms
annual_storms = dataset(idx,:);

% writetable(annual_storms, 'annual_storms.xlsx');
% writetable(dataset_30, 'PoT_30Peaks.xlsx');
% writetable(dataset_10, 'PoT_10Peaks.xlsx');


%% Plot of the measured wind direction for different heights at Cabauw for the selected peaks
% Extract the heights and their corresponding direction columns
heights = [10, 20, 40, 80, 140, 200];
direction_cols = {'D010', 'D020', 'D040', 'D080', 'D140', 'D200'};
datetime_col = 'DateTime';
% Extract datetime and wind direction data
datetimes = annual_storms.(datetime_col);
wind_directions = zeros(height(annual_storms), length(direction_cols));
for i = 1:length(direction_cols)
    wind_directions(:, i) = annual_storms.(direction_cols{i});
end

figure;
hold on;
colors = lines(length(heights));
for i = 1:length(heights)
    plot(datetimes, wind_directions(:, i), '-o', 'Color', colors(i, :), 'DisplayName', sprintf('%d m', heights(i)));
end
xlabel('Year');
ylabel('Wind direction (°)');
title('Wind directions at different heights during annual storms');
legend Location northwest
grid on;
xtickangle(45);

%% Plot of the measured wind direction for different heights at Cabauw for the 30% selected peaks
% Extract the heights and their corresponding direction columns
heights = [10, 20, 40, 80, 140, 200];
direction_cols = {'D010', 'D020', 'D040', 'D080', 'D140', 'D200'};
datetime_col = 'DateTime';
% Extract datetime and wind direction data
datetimes = dataset_30.(datetime_col);
wind_directions = zeros(height(dataset_30), length(direction_cols));
for i = 1:length(direction_cols)
    wind_directions(:, i) = dataset_30.(direction_cols{i});
end

figure;
hold on;
colors = lines(length(heights));
for i = 1:length(heights)
    plot(datetimes, wind_directions(:, i), 'o', 'Color', colors(i, :), 'DisplayName', sprintf('%d m', heights(i)));
end
xlabel('Measurement timestamp');
ylabel('Wind direction (°)');
title('Wind directions at different heights for 30% max peaks');
legend Location northwest
grid on;
xtickangle(45);

%% Plot of the measured wind direction for different heights at Cabauw for the 10% selected peaks
% Extract the heights and their corresponding direction columns
heights = [10, 20, 40, 80, 140, 200];
direction_cols = {'D010', 'D020', 'D040', 'D080', 'D140', 'D200'};
datetime_col = 'DateTime';
% Extract datetime and wind direction data
datetimes = dataset_10.(datetime_col);
wind_directions = zeros(height(dataset_10), length(direction_cols));
for i = 1:length(direction_cols)
    wind_directions(:, i) = dataset_10.(direction_cols{i});
end

figure;
hold on;
colors = lines(length(heights));
for i = 1:length(heights)
    plot(datetimes, wind_directions(:, i), 'o', 'Color', colors(i, :), 'DisplayName', sprintf('%d m', heights(i)));
end
xlabel('Measurement timestamp');
ylabel('Wind direction (°)');
title('Wind directions at different heights for 10% max peaks');
legend Location northwest
grid on;
xtickangle(45);


%% Sigma-analysis z0
% Computed local roughness lenght at Cabauw with sigma-analysis (KNMI)
sigma_analysis_z0.direction = [0;5;25;45;65;85;105;125;145;165;185;205;225;245;265;285;305;325;345;361];
sigma_analysis_z0.z0 = [0.138;0.125;0.149;0.151;0.179;0.191;0.183;0.090;0.072;0.070;0.067;0.076;0.089;0.136;0.139;0.13;0.192;0.150;0.138;0.138];

ws_measure_BM = [annual_storms.F010, annual_storms.F020, annual_storms.F040, annual_storms.F080, annual_storms.F140, annual_storms.F200];
vb = annual_storms.F010;  % Wind speed at 10m for each storm
wind_direction = annual_storms.D010;  % Wind direction at 10m for each storm
heights = [10, 20, 40, 80, 140, 200];

% Loop through each storm
for i = 1:length(vb)
    z0_index = find(sigma_analysis_z0.direction <= wind_direction(i) & wind_direction(i) < [sigma_analysis_z0.direction(2:end); 361], 1, 'first');
    sigma_z0(i) = sigma_analysis_z0.z0(z0_index);
    
    % Calculate mean wind speed at each height using the simplified formula
    % (where z0 = z0,ref)
    for j = 1:length(heights)
        z = heights(j);
        a = 1/log(10/sigma_z0(i));
        cr(i,j) = a*log(z/sigma_z0(i));
        vm(i,j) = cr(i,j)*vb(i);  % co = 1, so omitted here
    end
end


%% Yearly storms at Cabauw (measurements) with mean and max
mean_ws = mean(ws_measure_BM(:,:));
max_ws = max(ws_measure_BM(:,:));

heights = [10,20,40,80,140,200];
heights_log = log(heights);

figure
plot(mean_ws,heights_log,'k-',DisplayName='mean')
hold on
plot(max_ws,heights_log,'k--',DisplayName='max')
hold on
p = plot(ws_measure_BM(:,:),heights_log,LineStyle="none",Marker=".",MarkerSize=11);
legend([p(1)],'measurements')
hold on
coeffs = polyfit(heights_log,mean_ws,1);
data_fit = polyval(coeffs, heights_log);
plot(data_fit,heights_log, '-', Color='#9b9b9b',DisplayName='linear fit mean');
title('Yearly storms at Cabauw (measurements)')
xlabel('Wind speed (m/s)')
ylabel('ln(Height)')
legend Location northwest
% set(gca, 'YScale', 'log')

% Assuming logarithmic wind profile
figure
hold on
num_storms = size(vm, 1);
colors = lines(num_storms);
for i = 1:num_storms
    plot(vm(i,:), heights, 'Color', colors(i, :), 'DisplayName', ['Storm ', num2str(i)]);
    plot(ws_measure_BM(i,:), heights, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 11, 'Color', colors(i, :), 'HandleVisibility', 'off');
end
xlabel('Wind Velocity (m/s)')
ylabel('Height (m)')
title('Measurements vs computed vm for every storm')
legend('Location', 'bestoutside')
set(gca, 'YScale', 'log')
grid minor
hold off

% For the first storm (2001)
figure
plot(vm(1,:),heights,DisplayName='log profile')
hold on
% Measurements
plot(ws_measure_BM(1,:),heights,DisplayName='measurements', LineStyle='none',Marker='.',MarkerSize=11)
% hold on
% plot(return_values_measure_model,heights)
xlabel('Wind Velocity (m/s)')
ylabel('Height (m)')
title('Measurements vs computed vm for storm 1')
legend Location northwest
set(gca, 'YScale', 'log')
grid minor

%% 20year extremes
figure
plot(ws_measure_BM(:,:), heights, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 11, 'Color', colors(i, :));
hold on
plot(vm_20,z,'k-',LineWidth=1) % vm_20 computated in compare_vm.mat
xlabel('Mean Wind Velocity (m/s)')
ylabel('Height (m)')
title(['20-year extreme, z0 = ', num2str(z0)])
xlim([16.7435 29.9029])
ylim([10 200])
set(gca, 'YScale', 'log')
legend('v_{20,Cabauw}','Measurements')
legend Location northwest
grid minor


%% Log-linear regression model : MEAN

data = ws_measure_BM(:,2:6);
heights = [20, 40, 80, 140, 200]; % only take measurements made h >= 20m
heights_log = log(heights);
mean_data = mean(ws_measure_BM(:,2:6));

% Log-linear regression model (mean of data)
figure
hold on
for i = 1:size(data, 2)
    p = scatter(data(:, i), repmat(heights(i), size(data, 1), 1), 'Marker', '.', 'MarkerEdgeColor', '#0072BD');
    legend([p(1)],'Measurements')
end

% Fit a linear model for log-height vs mean wind speed
coeffs = polyfit(mean_data, heights_log, 1);
heights_log_fit = polyval(coeffs, mean_data);
z0_log = coeffs(1) * 0 + coeffs(2);  % ln(z0) corresponds to wind speed = 0
z0_mean_ws = exp(z0_log);
fprintf('Estimated roughness length (z0): %.4f m\n', z0_mean_ws);

% Extrapolation below the lowest wind speed to z0
mean_data_ext = linspace(min(mean_data), 0, 100);
heights_log_ext = polyval(coeffs, mean_data_ext);
heights_ext = exp(heights_log_ext);

% Plot the regression line and its extrapolation
plot(mean_data, exp(heights_log_fit), '-', 'Color', '#9b9b9b', 'DisplayName', 'Log-linear fit');
plot(mean_data_ext, heights_ext, '--', 'Color', '#9b9b9b', 'DisplayName', 'Extrapolation');
xlabel('Wind speed (m/s)');
ylabel('Height (m)');
title('Log-Linear plot of wind speed vs height - mean ws');
legend Location northwest
grid on
ylim([z0_log 200])
set(gca, 'YScale', 'log');

%% Log-linear regression model : EVERY STORM
data = [];
heights_log_ext = [];
data_fit = [];
data_fit_ext = [];

for i = 1:length(ws_measure_BM)
    data(i,:) = ws_measure_BM(i,:);
    heights = [10, 20, 40, 80, 140, 200];
    heights_log = log(heights);
    
    % Log-linear regression model (extrapolation)
    % mean_wind_speed = a⋅ln(height) + b
    coeffs(i,:) = polyfit(heights_log,data(i,:),1);
    a(i) = coeffs(i,1);
    b(i) = coeffs(i,2);
    
    % using y = ax + b, where wind speed (y) = 0:
    z0_log(i) = -b(i)/a(i);
    z0_ext(i) = exp(z0_log(i));
    
    % Display the roughness length
    fprintf('Estimated roughness length (z0): %.4f m\n', z0_ext(i));
    
    % Extrapolate to find the height where wind speed = 0
    heights_log_ext(i,:) = linspace(min(heights_log), z0_log(i), 50);
    data_fit(i,:) = polyval(coeffs(i,:), heights_log);
    data_fit_ext(i,:) = polyval(coeffs(i,:),heights_log_ext(i,:));
end

figure
p1 = plot(data,exp(heights_log),LineStyle="none",Marker=".",MarkerSize=11,MarkerEdgeColor="#0072BD");
hold on;
p2 = plot(data_fit,exp(heights_log), '-', Color='#9b9b9b');
hold on
p3 = plot(data_fit_ext,exp(heights_log_ext), '.', Color='#9b9b9b');
legend([p1(1) p2(1) p3(1)],'Measurements','log-linear fit','extrapolation')

ylabel('Height (m)');
xlabel('Wind speed (m/s)');
legend Location northwest;
title('Log-Linear plot with extrapolation');
grid on
xlim([0 34.027])
ylim([min(z0_log) 200])
set(gca, 'YScale', 'log');

mean_z0 = mean(z0_ext);
fprintf('Mean of the estimated roughness length: %.4f m\n', mean_z0);

%% NON-LINEAR REGRESSION MODEL (WITH DISPLACEMENT HEIGHT)
% Define initial guesses for a, b, and d
initial_guess = [0.1, 0.1, 5]; % [a, b, d]
heights = [10, 20, 40, 80, 140, 200];

% Initialize arrays for storing results
num_storms = size(ws_measure_BM, 1);
z0_values = zeros(num_storms, 1);
z0_uncertainties = zeros(num_storms, 1);
d_values = zeros(num_storms, 1);
d_uncertainties = zeros(num_storms, 1);

figure;
hold on;
color_map = lines(num_storms); % Generate a colormap for distinct storm colors
for i = 1:num_storms
    data = ws_measure_BM(i, :);

    % Define the log-linear model with displacement height
    log_linear_model = @(params, h) params(1) * log(h - params(3)) + params(2);

    % Nonlinear fitting using lsqcurvefit
    options = optimset('Display', 'off');
    lb = [-Inf, -Inf, 0]; % Lower bounds: d must be non-negative
    ub = [Inf, Inf, min(heights)]; % Upper bounds: d must be smaller than min(heights)
    [params, resnorm, residual, exitflag, output, lambda, jacobian] = lsqcurvefit(... 
        log_linear_model, initial_guess, heights, data, lb, ub, options);

    % Extract coefficients
    a = params(1);
    b = params(2);
    d = params(3);

    % Calculate roughness length z0
    z0_log = -b / a + d;
    z0_ext_with_d = exp(z0_log);

    % Calculate uncertainties
    if issparse(jacobian)
        jacobian = full(jacobian);
    end
    cov_matrix = inv(jacobian' * jacobian);
    uncertainty = sqrt(diag(cov_matrix));
    partial_a = b / a^2;
    partial_b = -1 / a;
    partial_d = 1;
    delta_z0_log = sqrt((partial_a * uncertainty(1))^2 + (partial_b * uncertainty(2))^2 + (partial_d * uncertainty(3))^2);
    delta_z0 = z0_ext_with_d * delta_z0_log;

    % Store results
    z0_values(i) = z0_ext_with_d;
    z0_uncertainties(i) = delta_z0;
    d_values(i) = d;
    d_uncertainties(i) = uncertainty(3);

    % Plot measurements
    plot(data, heights,'o','MarkerSize',8,'MarkerEdgeColor',color_map(i, :),'DisplayName',sprintf('Storm %d - Measurements', i));

    % Constrain extrapolation range
    max_wind_speed = max(data);
    wind_speed_ext = linspace(0, max_wind_speed, 100);
    heights_log_ext = (wind_speed_ext - b) / a; % Solve for log(heights - d)
    heights_ext = exp(heights_log_ext) + d; % Convert back to original height scale

    plot(wind_speed_ext, heights_ext, '-', 'Color', color_map(i, :), 'DisplayName', sprintf('Storm %d - Extrapolation', i));
end
xlabel('Wind speed (m/s)');
ylabel('Height (m)');
set(gca, 'YScale', 'log');
grid on;
legend('Location', 'northwest');
title('Log-Linear Plot with Displacement Height for All Storms');

% Create a table to display the results
results_table = table((1:num_storms)', d_values, d_uncertainties, z0_values, z0_uncertainties, 'VariableNames', {'Storm', 'd (Displacement Height)', 'Uncertainty in d', 'z0 (Roughness Length)', 'Uncertainty in z0'});
disp(results_table);
% writetable(results_table, 'results.csv')


%% difference between z0 extrapolated with measurements and z0 from sigma-analysis
z0_eurocodes = [0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];

% relative difference
model_error = z0_ext / sigma_z0;
Eurocodes_error = z0_ext / z0_eurocodes;

disp(model_error)
disp(Eurocodes_error)

% absolut difference
model_error = sigma_z0 - z0_ext;
Eurocodes_error = z0_eurocodes - z0_ext;

% histograms
figure;
hold on;
histogram(model_error, floor(sqrt(length(z0_ext))), 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Model Error');
histogram(Eurocodes_error, floor(sqrt(length(z0_ext))), 'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Eurocodes Error');
xlabel('Error');
ylabel('Frequency');
legend Location northwest;
title('Overlayed Histogram: Model Error vs Eurocodes Error');
grid on;

% boxplots
figure;
boxplot([model_error', Eurocodes_error'], 'Labels', {'Model Error', 'Eurocodes Error'});
ylabel('Error');
title('Boxplot Comparison of Errors');
grid on;

% plot of z0/storm
mean_z0_sigma = mean(sigma_z0);
mean_z0_ext = mean(z0_ext);
figure
x = linspace(2001,2019,19); 
plot(x,sigma_z0,LineStyle="none",Marker="o")
hold on
plot(x,z0_ext,LineStyle="none",Marker="o")
hold on
yline(mean_z0_sigma,'--')
hold on
yline(mean_z0_ext,'--')
xlabel('Annual storm (year)')
ylabel('Roughness length (m)')
legend('sigma-analysis','extrapolated')


%% Roughness Length by Wind Direction Range
% Define initial guesses for a, b, and d
initial_guess = [0.1, 0.1, 5]; % [a, b, d]
heights = [10, 20, 40, 80, 140, 200];
direction_ranges = [5, 25; 25, 45; 45, 65; 65, 85; 85, 105; 105, 125; ...
                    125, 145; 145, 165; 165, 185; 185, 205; ...
                    205, 225; 225, 245; 245, 265; 265, 285; ...
                    285, 305; 305, 325; 325, 345; 345, 360];

num_directions = size(direction_ranges, 1);
used_z0_by_category = cell(num_directions, 1); % sigma-analysis
z0_ext_by_category = cell(num_directions, 1);  % regression model : mean_wind_speed = a⋅ln(height) + b
z0_computed_by_category = cell(num_directions, 1); % regression model : mean_wind_speed = a⋅ln(height - d) + b

% Loop through each storm
for i = 1:length(vb)
    for j = 1:num_directions
        if wind_direction(i) >= direction_ranges(j, 1) && wind_direction(i) < direction_ranges(j, 2)
            used_z0_by_category{j} = [used_z0_by_category{j}, sigma_z0(i)];
            z0_ext_by_category{j} = [z0_ext_by_category{j}, z0_ext(i)];
            z0_computed_by_category{j} = [z0_computed_by_category{j}, z0_values(i)]; 
            break;
        end
    end
end

direction_labels = {'5-24°', '25-44°', '45-64°', '65-84°', '85-104°', ...
                    '105-124°', '125-144°', '145-164°', '165-184°', ...
                    '185-204°', '205-224°', '225-244°', '245-264°', ...
                    '265-284°', '285-304°', '305-324°', '325-344°', '345-4°'};
x_positions = 1:num_directions;

figure;
hold on;
for j = 1:num_directions
    % Scatter for sigma_z0
    scatter(repelem(x_positions(j), length(used_z0_by_category{j})), used_z0_by_category{j}, ...
        'o', 'DisplayName', 'used\_z0', 'MarkerEdgeColor', 'b');
    
    % Scatter for z0_ext
    scatter(repelem(x_positions(j), length(z0_ext_by_category{j})), z0_ext_by_category{j}, ...
        'x', 'DisplayName', 'z0\_ext', 'MarkerEdgeColor', 'r');

    % Scatter for z0_ext_with_d
    scatter(repelem(x_positions(j), length(z0_computed_by_category{j})), z0_computed_by_category{j}, ...
        '+', 'DisplayName', 'z0\_ext', 'MarkerEdgeColor', 'k');
    
    % Add Eurocodes_z0 value (constant at 0.2) for each direction
    plot([j - 0.2, j + 0.2], [0.2, 0.2], 'k-', 'LineWidth', 1, 'DisplayName', 'Eurocodes\_z0');
end
xticks(x_positions);
xticklabels(direction_labels);
xlabel('Wind Direction Range (°)');
ylabel('Roughness Length (m)');
title('Roughness Length by Wind Direction Range');
legend({'sigma-analysis', 'extrapolated','extrapolated with d', 'Eurocodes'}, 'Location', 'southeast');
grid on;
set(gca, 'YScale', 'log')
xlim([10, 17]);




