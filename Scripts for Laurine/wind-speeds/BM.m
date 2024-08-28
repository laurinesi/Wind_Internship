% exploration of Schiphol dataset

addpath('..\wind-speeds\')
data_Schiphol = readmatrix("data_Schiphol.txt"); % measurements
model_Schiphol = readmatrix("model_Schiphol.txt"); % model

% Block Maxima = max wind speed of each year

years = unique(data_Schiphol(:,1)); % Find unique years
max_values = zeros(length(years), 2); % Preallocate for speed

for i = 1:length(years)
    year = years(i);
    % Find the maximum value for the current year
    max_value = max(data_Schiphol(data_Schiphol(:,1) == year, 2));
    max_values(i,:) = [year, max_value]; % Store the year and max value
end

% Create a table with appropriate column names
FF_BM = array2table(max_values, 'VariableNames', {'Year', 'FF'});

% Save the table as a CSV file
writetable(FF_BM, 'Schiphol_BM.csv');

% Save the table as a text file
writetable(FF_BM, 'Schiphol_BM.txt', 'Delimiter', '\t');

% BM on all dataset (with the 1st July = new year)
figure
plot(data_Schiphol(:,1), data_Schiphol(:,2),'Marker','o')
hold on
plot(max_values(:,1),max_values(:,2),'Marker','o','LineStyle','none','Color','r','MarkerFaceColor','r')
xlabel('Year');
ylabel('Wind speeds (m/s)');
title('Measurements at Schiphol');
yticks(0:4:30)
xlim([1950 2023])

% BM values
figure
plot(max_values(:,1),max_values(:,2),'Marker','x','LineStyle','none')
xlabel('Year');
ylabel('Wind speeds (m/s)');
title('Yearly maximum wind speeds');
yticks(15:2:27)

% histogram of all the measurements
% figure
% histogram(data_Schiphol(:,2),150)
% hold on
%plot une  ligne vertcial à 14.4 et à 27.5 pour montrer ce qu'on utilise
%apres (les maximums)


%% Distribution fitting & parameter estimation

% Sort data for ECDF
sorted_data = sort(max_values(:,2));
n = length(max_values);
ecdf = (1:n) / n;

% MLE

% Normal MLE
phat = mle(max_values(:,2),'Distribution','Normal');
mean_norm = phat(1);
std_norm = phat(2);

% Gumbel MLE
parmhat = evfit(max_values(:,2));
location_gum = parmhat(1);
scale_gum = parmhat(2);

% Uniform MLE
phat = mle(max_values(:,2),'Distribution','Uniform');
min_unif = phat(1);
max_unif = phat(2);

% Log-normal MLE
phat = mle(max_values(:,2),'Distribution','Lognormal');
mean_lognorm = phat(1);
std_lognorm = phat(2);

% Rayleigh MLE
phat = mle(max_values(:,2),'Distribution','Rayleigh');
scale_ray = phat(1);

% GEV MLE
phat = mle(max_values(:,2),'Distribution','Generalized Extreme Value');
shape_gev = phat(1);
scale_gev = phat(2);
location_gev = phat(3);

fprintf('MLE normal - mu: %.4f, sigma: %.4f\n', mean_norm, std_norm);


%% Normal plot

% Generate x values for plotting
x = linspace(12, 30, 1000);

% Plot ECDF
figure;
stairs(sorted_data, ecdf, 'Marker', 'o', 'LineStyle', 'none', 'DisplayName', 'ECDF', LineWidth=1);
hold on;

% Plot CDFs of fitted distributions

% Normal CDF
cdf_normal = normcdf(x, mean_norm, std_norm);
plot(x, cdf_normal, 'DisplayName', 'Normal CDF', 'LineStyle', '-', LineWidth=1);

% Gumbel CDF
% cdf_gumbel = exp(-exp(-alpha*(x-u)));
cdf_gumbel = evcdf(x, location_gum, scale_gum);
plot(x, cdf_gumbel, 'DisplayName', 'Gumbel CDF', 'LineStyle', '--', LineWidth=1);

% Uniform CDF
cdf_uniform = unifcdf(x, min_unif, max_unif);
plot(x, cdf_uniform, 'DisplayName', 'Uniform CDF', 'LineStyle', '--', LineWidth=1);

% Lognormal CDF
cdf_lognorm = logncdf(x, mean_lognorm, std_lognorm);
plot(x, cdf_lognorm, 'DisplayName', 'Lognormal CDF', 'LineStyle', '--', LineWidth=1);

% Rayleigh CDF
cdf_rayleigh = raylcdf(x, scale_ray);
plot(x, cdf_rayleigh, 'DisplayName', 'Rayleigh CDF', 'LineStyle', '--', LineWidth=1);

% GEV CDF
cdf_gev = gevcdf(x,shape_gev,scale_gev,location_gev);
plot(x, cdf_gev, 'DisplayName', 'GEV CDF', 'LineStyle', '--', LineWidth=1);


% Labels and legend
xlabel('x');
ylabel('P(X <= x)');
title('ECDF and Fitted CDFs');
legend('show', Location='southeast');
legend box off;
grid on;
hold off;


%% Gumbel plot

% Generate x values for plotting
x = linspace(12, 30, 1000);

% Plot ECDF
figure;
stairs(sorted_data, ecdf, 'Marker', 'o', 'LineStyle', 'none', 'DisplayName', 'ECDF', LineWidth=1);
hold on;

% Plot CDFs of fitted distributions

% Normal CDF
cdf_normal = normcdf(x, mean_norm, std_norm);
log_cdf_normal = -log(-log(cdf_normal));
plot(x, log_cdf_normal, 'DisplayName', 'Normal CDF', 'LineStyle', '-', LineWidth=1);

% Gumbel CDF
cdf_gumbel = evcdf(x, location_gum, scale_gum);
log_cdf_gumbel = -log(-log(cdf_gumbel));
plot(x, log_cdf_gumbel, 'DisplayName', 'Gumbel CDF', 'LineStyle', '--', LineWidth=1);

% Uniform CDF
cdf_uniform = unifcdf(x, min_unif, max_unif);
plot(x, cdf_uniform, 'DisplayName', 'Uniform CDF', 'LineStyle', '--', LineWidth=1);

% Lognormal CDF
cdf_lognorm = logncdf(x, mean_lognorm, std_lognorm);
log_cdf_lognorm = -log(-log(cdf_lognorm));
plot(x, log_cdf_lognorm, 'DisplayName', 'Lognormal CDF', 'LineStyle', '--', LineWidth=1);

% Rayleigh CDF
cdf_rayleigh = raylcdf(x, scale_ray);
log_cdf_rayleigh = -log(-log(cdf_rayleigh));
plot(x, log_cdf_rayleigh, 'DisplayName', 'Rayleigh CDF', 'LineStyle', '--', LineWidth=1);


% Labels and legend
xlabel('x');
ylabel('-log(-log(P(X <= x)))');
title('ECDF and Fitted CDFs');
legend('show', Location='northwest');
legend box off;
grid on;
hold off;


