clear all
close all
clc

addpath('..\wind-speeds\tools\')
run("C_windspeed_datasets.m") % Cabauw measure & model (RACMO, KNW)

% working with the measurement data at Cabauw
dataset = measurements_Cabauw;

%% get 30% highest wind speeds
valid_idx = find(~isnan(dataset.F010));
valid_dataset = dataset(valid_idx, :);
reorder_dataset = sortrows(valid_dataset,"F010","descend");
new_idx = find(~isnan(reorder_dataset.F010));

index_30 = 30*length(reorder_dataset.Year)/100;
dataset_30 = reorder_dataset(1:index_30, :);

%% get the yearly maximum wind speeds
% Initialization
idx = [];
z0 = [];
a = [];
cr = [];
vm = [];
direction_z0 = table;

dataset = measurements_Cabauw;

% get annual maximum
[max_values, max_datetime] = BM_select(dataset);
% find line of max_datetime in dataset
for i = 1:length(max_datetime)
    idx(i,:) = find(dataset.DateTime==max_datetime(i));
end
% store annual storms
annual_storms = dataset(idx,:);


% Measured local roughness lenght at Cabauw
direction_z0.direction = [0;5;25;45;65;85;105;125;145;165;185;205;225;245;265;285;305;325;345;361];
direction_z0.z0 = [0.138;0.125;0.149;0.151;0.179;0.191;0.183;0.090;0.072;0.070;0.067;0.076;0.089;0.136;0.139;0.13;0.192;0.150;0.138;0.138];

ws_measure_BM = [annual_storms.F010, annual_storms.F020, annual_storms.F040, annual_storms.F080, annual_storms.F140, annual_storms.F200];
vb = annual_storms.F010;  % Base wind speed at 10m for each storm
wind_direction = annual_storms.D010;  % Wind direction at 10m for each storm
heights = [10, 20, 40, 80, 140, 200];


% Loop through each storm
for i = 1:length(vb)
    z0_index = find(direction_z0.direction <= wind_direction(i) & wind_direction(i) < [direction_z0.direction(2:end); 361], 1, 'first');
    z0(i) = direction_z0.z0(z0_index);
    
    % Calculate mean wind speed at each height using the simplified formula
    for j = 1:length(heights)
        z = heights(j);
        a(i) = 1/log(10/z0(i));
        cr(i,j) = a(i)*log(z/z0(i));
        vm(i,j) = cr(i,j)*vb(i);  % co = 1, so omitted here
    end
end

% Yearly storms at Cabauw (measurements) with mean and max
mean_ws = mean(ws_measure_BM(:,:));
max_ws = max(ws_measure_BM(:,:));
figure
plot(mean_ws,heights,'k-')
hold on
plot(max_ws,heights,'k--')
hold on
plot(ws_measure_BM(:,:),heights,LineStyle="none",Marker="x")
title('Yearly storms at Cabauw (measurements)')
legend('mean','max')
set(gca, 'YScale', 'log')


%% Log-linear regression model

data = ws_measure_BM(1,:);
heights = [10, 20, 40, 80, 140, 200];
heights_log = log10(heights);

% Log-linear regression model (only data)
figure
scatter(heights_log,data)
lsline

% Log-linear regression model (extrapolation)
% mean_wind_speed = a⋅log_10(height) + b
coeffs = polyfit(heights_log,data,1);
a = coeffs(1);
b = coeffs(2);

% using y = ax + b, where wind speed (y) = 0:
z0_log = -b/a;
z0 = 10^z0_log;
% Display the roughness length
fprintf('Estimated roughness length (z0): %.4f m\n', z0);

% Extrapolate to find the height where wind speed = 0
heights_log_ext = linspace(min(heights_log), z0_log, 100);
data_fit = polyval(coeffs, heights_log);
data_fit_ext = polyval(coeffs,heights_log_ext);

figure;
plot(heights_log, data, 'o', 'DisplayName', 'data');
hold on;
plot(heights_log, data_fit, '-', Color='#9b9b9b', DisplayName='log-linear fit');
hold on
plot(heights_log_ext, data_fit_ext, '--', Color='#9b9b9b', DisplayName='extrapolation');
xlabel('log10(Height)');
ylabel('Wind speed (m/s)');
legend Location northwest;
title('Log-Linear plot of wind speed vs height with extrapolation');
grid on



%%
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
% title('Measurements vs computed vm for every storm')
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
% title('Measurements vs computed vm for storm 1')
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

