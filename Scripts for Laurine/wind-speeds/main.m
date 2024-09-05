% probabilistic description of wind loads

addpath('..\wind-speeds\tools\')

% Datasets
run("S_windspeed_datasets.m") % Schiphol measure & model (RACMO)
% run("C_windspeed_datasets.m") % Cabauw measure & model (RACMO, KNW)

%% GEV - BM - MLE - Bootstrap

% tail index (shape) : GEV MLE on weather model dataset (BM)

% get the yearly maximum wind speeds
dataset = model_Schiphol;
[max_values] = BM_select(dataset);

[parmhat, ~,se] = gevfit2(max_values(:,2));

tail = parmhat(1);
scale_model = parmhat(2);
location_model = parmhat(3);
se_tail = se(1);
se_scale = se(2);
se_loc = se(3);

disp(['MLE GEV model - tail index: ', num2str(tail)]);
disp(['MLE GEV model - scale: ', num2str(scale_model)]);
disp(['MLE GEV model - location: ', num2str(location_model)]);

% scale & location : bootstrap on measurements dataset (BM)

dataset = data_Schiphol;
[max_values] = BM_select(dataset);

[parmhat, parmci, se] = gevfit2(max_values(:,2));
tail_measure = parmhat(1);
scale_measure = parmhat(2);
location_measure = parmhat(3);
std_tail_measure = se(1);
std_scale_measure = se(2);
std_location_measure = se(3);

disp(['MLE GEV - shape: ', num2str(tail_measure)]);
disp(['MLE GEV - scale: ', num2str(scale_measure)]);
disp(['MLE GEV - location: ', num2str(location_measure)]);

population = max_values(:,2); % original population
n = 50; % number of bootstrap samples

% need to fix the tail in gevfit for bootstrap
[GEVparameters] = bootstrap(population, n, tail);


% Calculate statistics for each parameter
mean_shape = tail;
std_shape = se_tail;

% 'mean' = ML estimates
% if we take infinite bootstrap samples, the mean of the parameters should
% approximate the ML estimates
mean_scale = scale_measure;
mean_location = location_measure;

% take parmci to plot 95% confidence instead of mean +/- std
std_scale = std(GEVparameters.scale);
std_location = std(GEVparameters.location);

% 95% confidence 
tail_confidence = parmci(:,1);
scale_confidence =  parmci(:,2);
location_confidence = parmci(:,3);

% Display calculated statistics
fprintf('shape parameters - mu: %.4f, sigma: %.4f\n', mean_shape, std_shape);
fprintf('scale parameters - mu: %.4f, sigma: %.4f\n', mean_scale, std_scale);
fprintf('location parameters - mu: %.4f, sigma: %.4f\n', mean_location, std_location);


%%
% Plot normal probability plot for each parameter
% figure;
% normplot(tail);
% title('Normal Probability Plot (shape)');

figure;
normplot(GEVparameters.scale);
title('Normal Probability Plot (scale)');

figure;
normplot(GEVparameters.location);
title('Normal Probability Plot (location)');

% Plot PDF of each parameter assuming normal distribution
x_shape = linspace(mean_shape - 4*std_shape, mean_shape + 4*std_shape, 100); % Range for shape parameter
x_scale = linspace(mean_scale - 4*std_scale, mean_scale + 4*std_scale, 100); % Range for scale parameter
x_location = linspace(mean_location - 4*std_location, mean_location + 4*std_location, 100); % Range for location parameter

pdf_shape = normpdf(x_shape, mean_shape, std_shape); % PDF of shape parameter
pdf_scale = normpdf(x_scale, mean_scale, std_scale); % PDF of scale parameter
pdf_location = normpdf(x_location, mean_location, std_location); % PDF of location parameter

% shape
figure;
plot(x_shape, pdf_shape, 'Color',	"#A2142F", 'LineWidth', 1.5);
title('PDF of Shape Parameter');
xlabel('Shape');
ylabel('Density');
xlim([min(x_shape), max(x_shape)])
hold on;
% Draw a vertical line at the mean
plot([mean_shape, mean_shape], [0, max(pdf_shape)], 'k--', 'LineWidth', 1);
% Shade the region representing ±1 std deviation
x_fill = [mean_shape-std_shape, linspace(mean_shape-std_shape, mean_shape+std_shape, 100), mean_shape+std_shape];
y_fill = [0, normpdf(linspace(mean_shape-std_shape, mean_shape+std_shape, 100), mean_shape, std_shape), 0];
fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;

% scale
figure;
plot(x_scale, pdf_scale, 'Color',	"#0072BD", 'LineWidth', 1.5);
title('PDF of Scale Parameter');
xlabel('Scale');
ylabel('Density');
xlim([min(x_scale), max(x_scale)])
hold on;
% Draw a vertical line at the mean
plot([mean_scale, mean_scale], [0, max(pdf_scale)], 'k--', 'LineWidth', 1);
    % Draw a vertical line at scale estimate from weather model
    plot([scale_model, scale_model], [0, max(pdf_scale)], 'k--', 'LineWidth', 1);
    % Draw a vertical line at scale estimate
    plot([scale_measure, scale_measure], [0, max(pdf_scale)], 'k--', 'LineWidth', 1);
%Shade the region representing ±1 std deviation
x_fill = [mean_scale-std_scale, linspace(mean_scale-std_scale, mean_scale+std_scale), mean_scale+std_scale];
y_fill = [0, normpdf(linspace(mean_scale-std_scale, mean_scale+std_scale), mean_scale, std_scale), 0];
fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;

% location
figure;
plot(x_location, pdf_location, 'Color',	"#77AC30", 'LineWidth', 1.5);
title('PDF of Location Parameter');
xlabel('Location');
ylabel('Density');
xlim([min(x_location), max(x_location)])
hold on;
% Draw a vertical line at the mean
plot([mean_location, mean_location], [0, max(pdf_location)], 'k--', 'LineWidth', 1);
%     % Draw a vertical line at the mean
%     plot([location_model, location_model], [0, max(pdf_location)], 'k--', 'LineWidth', 1);
%     % Draw a vertical line at the mean
%     plot([location_measure, location_measure], [0, max(pdf_location)], 'k--', 'LineWidth', 1);
% Shade the region representing ±1 std deviation
x_fill = [mean_location-std_location, linspace(mean_location-std_location, mean_location+std_location), mean_location+std_location];
y_fill = [0, normpdf(linspace(mean_location-std_location, mean_location+std_location), mean_location, std_location), 0];
fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;

%% table to save data

% Assume these are your variables
tail_model = {[num2str(tail) ', sd=' num2str(se_tail)]};
tail_measurements = {[num2str(tail_measure) ', sd=' num2str(std_tail_measure)]};
tail_model_measurements = {[num2str(tail) ', sd=' num2str(std_shape)]};

scale_model = {[num2str(scale_model) ', sd=' num2str(se_scale)]};
scale_measurements = {[num2str(scale_measure) ', sd=' num2str(std_scale_measure)]};
scale_model_measurements = {[num2str(mean_scale) ', sd=' num2str(std_scale)]};

location_model = {[num2str(location_model) ', sd=' num2str(se_loc)]};
location_measurements = {[num2str(location_measure) ', sd=' num2str(std_location_measure)]};
location_model_measurements = {[num2str(mean_location) ', sd=' num2str(std_location)]};

% Create the cell array by storing your variables
table = {
    tail_model, tail_measurements, tail_model_measurements;
    scale_model, scale_measurements, scale_model_measurements;
    location_model, location_measurements, location_model_measurements
    };


% Create a table with appropriate column names
results = array2table(table, 'VariableNames', {'model', 'measurements', 'model + measurements'}, 'RowNames',{'tail index', 'scale', 'location'});

writetable(results, 'results.txt', 'WriteRowNames', true, 'Delimiter', '\t');

writetable(results, 'results.xlsx', 'WriteRowNames', true);




%% GP - PoT - MLE - Analytical unc


%% GW - PoT - iHill(i) - Analytical unc


%% log-GW - PoT - iHill(i) - Analytical unc

