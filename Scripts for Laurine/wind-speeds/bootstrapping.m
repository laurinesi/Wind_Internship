% Bootstrapping on values obtained with the Block Maxima method and then
% fit a Generalized Extreme Value (GEV) distribution to estimate the
% distribution of the GEV fit parameters (location, scale, and shape).


addpath('..\wind-speeds\')
BM_Schiphol = readmatrix("Schiphol_BM.txt"); % measurements


% Step 1 : population = dataset of BM values assumed independent (74 years for measurements at Schiphol)
population = BM_Schiphol(:,2);

% Step 2 : resample nb of time
n = 10; % number of bootstrap samples

% Initialization
samples = zeros(length(population), n);
parmhat = zeros(n, 3);  % [k, sigma, mu] parameters of GEV

% Step 3 : generate bootstrap samples and fit GEV distribution to each
for i = 1:n
    % Generate a bootstrap sample with replacement
    samples(:, i) = datasample(population, length(population), 'Replace', true);
    
    % Fit GEV distribution to the bootstrap sample
    % gevfit(X) returns maximum likelihood estimates
    parmhat(i, :) = gevfit(samples(:, i));
end

GEVparameters = array2table(parmhat, "VariableNames",{'shape', 'scale', 'location'});

% Step 4 : display the GEV parameters for each bootstrap sample
disp('GEV Parameters for each bootstrap sample :');
disp(['Result for ', num2str(n), ' estimations'])
disp(['ndraw = ', num2str(length(population)), 'years'])
disp(GEVparameters);

% Step 5 : compute mean and standard deviation of each statistic 
% and draw their distribution (assumed normal)

% Calculate statistics for each parameter
mean_shape = mean(parmhat(:, 1));
std_shape = std(parmhat(:, 1));
mean_scale = mean(parmhat(:, 2));
std_scale = std(parmhat(:, 2));
mean_location = mean(parmhat(:, 3));
std_location = std(parmhat(:, 3));

% Display calculated statistics
fprintf('shape parameters - mu: %.4f, sigma: %.4f\n', mean_shape, std_shape);
fprintf('scale parameters - mu: %.4f, sigma: %.4f\n', mean_scale, std_scale);
fprintf('location parameters - mu: %.4f, sigma: %.4f\n', mean_location, std_location);

% Plot normal probability plot for each parameter
figure;
normplot(parmhat(:, 1));
title('Normal Probability Plot (shape)');

figure;
normplot(parmhat(:, 2));
title('Normal Probability Plot (scale)');

figure;
normplot(parmhat(:, 3));
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
% Shade the region representing ±1 std deviation
x_fill = [mean_scale-std_scale, linspace(mean_scale-std_scale, mean_scale+std_scale, 100), mean_scale+std_scale];
y_fill = [0, normpdf(linspace(mean_scale-std_scale, mean_scale+std_scale, 100), mean_scale, std_scale), 0];
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
% Shade the region representing ±1 std deviation
x_fill = [mean_location-std_location, linspace(mean_location-std_location, mean_location+std_location, 100), mean_location+std_location];
y_fill = [0, normpdf(linspace(mean_location-std_location, mean_location+std_location, 100), mean_location, std_location), 0];
fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold off;


