% probabilistic description of wind loads

addpath('..\wind-speeds\tools\')

% Datasets
run("S_windspeed_datasets.m") % Schiphol measure & model (RACMO)
% run("C_windspeed_datasets.m") % Cabauw measure & model (RACMO, KNW)

%settings
dataType = 2;  % [1 2] ~ [Model Measurements]
fixedTail = 2;  % [1 2] ~ [no yes]

%% GEV - BM - MLE - Bootstrap

% tail index (shape) : GEV MLE on weather model dataset (BM)

if dataType == 1
    % get the yearly maximum wind speeds
    dataset = model_Schiphol;
    [max_values] = BM_select(dataset);
    
    % fit a gev distribution
    [parmhat, parmci, se] = gevfit2(max_values(:,2));
    
    tail = parmhat(1);
    scale_model = parmhat(2);
    location_model = parmhat(3);
    se_tail = se(1);
    se_scale = se(2);
    se_loc = se(3);
    
    % settings
    population = max_values(:,2); % original population
    n = 1000; % number of bootstrap samples
    
    % bootstrap
    [GEVparameters] = bootstrap(population, n);
    
    % check normality of bootstrap estimates and compute quantiles for 95%
    [Qt, Qs, Ql] = check_normality(GEVparameters, parmhat, 1, 1);

    % computation COV
    cov_shape = -(std(GEVparameters.tail)/mean(GEVparameters.tail));
    cov_scale = std(GEVparameters.scale)/mean(GEVparameters.scale);
    cov_location = std(GEVparameters.location)/mean(GEVparameters.location);

% difference with/without bootstrapping 
    % boxplot of tail estimate using bootstrap method
    figure
    boxplot(GEVparameters.tail)
    hold on
    plot(tail,Marker="_",LineWidth=1,Color='k', DisplayName='MLE estimate')
    legend box on

%     % Create error bar plot
%     figure;
%     errorbar(1, tail, se_tail, 'o', 'MarkerSize', 10, 'LineWidth', 2);
%     xlim([0 20]);   % Set x-axis limits to keep it simple
%     ylim([tail - 2*se_tail, tail + 2*se_tail]); % Set y-axis limits
%     xticks(1);     % Only one x-tick for the single data point
%     xticklabels({'Data Point'});
%     ylabel('Value');
%     title('Single Value with Standard Error');
    
    % % MLE
    % disp(['MLE - tail index: ', num2str(tail)]);
    % disp(['MLE - scale: ', num2str(scale_model)]);
    % disp(['MLE - location: ', num2str(location_model)]);
    % % MLE Bootstrap
    % fprintf('Bootstrap - tail - mu: %.4f, sigma: %.4f\n', mean(GEVparameters.tail), std(GEVparameters.tail));
    % fprintf('Bootstrap - scale - mu: %.4f, sigma: %.4f\n', mean(GEVparameters.scale), std(GEVparameters.scale));
    % fprintf('Bootstrap - location - mu: %.4f, sigma: %.4f\n',mean(GEVparameters.location), std(GEVparameters.location));


    % Comparison of different methods for parameter estimation
    tail_mle = [num2str(tail) ', se=' num2str(se_tail)];
    tail_bootstrap = [num2str(mean(GEVparameters.tail)) ', sd=' num2str(std(GEVparameters.tail))];
    scale_mle = [num2str(scale_model) ', se=' num2str(se_scale)];
    scale_bootstrap = [num2str(mean(GEVparameters.scale)) ', sd=' num2str(std(GEVparameters.scale))];
    location_mle = [num2str(location_model) ', se=' num2str(se_loc)];
    location_bootstrap = [num2str(mean(GEVparameters.location)) ', sd=' num2str(std(GEVparameters.location))];
    
    table = {
        tail_mle, tail_bootstrap;
        scale_mle, scale_bootstrap;
        location_mle, location_bootstrap
        };
    
    MLE_boot = array2table(table, 'VariableNames', {'MLE', 'MLE - Bootstrap'}, 'RowNames',{'tail index', 'scale', 'location'});
    % writetable(MLE_boot, 'pestimation_model.xlsx', 'WriteRowNames', true);
    disp(MLE_boot)
    
    
    % Comparison of different methods to compute confidence interval
    ci_tail_mle = [num2str(parmci(1,1)) ', ' num2str(parmci(2,1))];
    ci_tail_bootstrap = [num2str(Qt(1)) ', ' num2str(Qt(2))];
    ci_scale_mle = [num2str(parmci(1,2)) ', ' num2str(parmci(2,2))];
    ci_scale_bootstrap = [num2str(Qs(1)) ', ' num2str(Qs(2))];
    ci_location_mle = [num2str(parmci(1,3)) ', ' num2str(parmci(2,3))];
    ci_location_bootstrap = [num2str(Ql(1)) ', ' num2str(Ql(2))];
    
    table = {
        ci_tail_mle, ci_tail_bootstrap;
        ci_scale_mle, ci_scale_bootstrap;
        ci_location_mle, ci_location_bootstrap
        };
    
    comp_ci = array2table(table, 'VariableNames', {'MLE', 'MLE - Bootstrap'}, 'RowNames',{'ci - tail index', ' ci - scale', 'ci - location'});
    % writetable(comp_ci, 'confint_model.xlsx', 'WriteRowNames', true);
    disp(comp_ci)

    % disp COV
    disp(['COV - shape: ', num2str(cov_shape)]);
    disp(['COV - scale: ', num2str(cov_scale)]);
    disp(['COV - location: ', num2str(cov_location)]);

end

%% scale & location : bootstrap on measurements dataset (BM)

if dataType == 2

    % Compute tail from weather model
    dataset = model_Schiphol;
    [max_values] = BM_select(dataset);
    [parmhat, ~, se] = gevfit2(max_values(:,2));
    tail = parmhat(1);
    se_tail = se(1);

    population = max_values(:,2); % original population
    n = 1000; % number of bootstrap samples

    % bootstrap
    [GEVparameters] = bootstrap(population, n);
  %  GEVparameters.tail = boot_tail;


    % Compute scale and location from measurements
    dataset = data_Schiphol;
    [max_values] = BM_select(dataset);
    
    [parmhat, parmci, se] = gevfit2(max_values(:,2));
    tail_measure = parmhat(1);
    scale_measure = parmhat(2);
    location_measure = parmhat(3);
    se_tail_measure = se(1);
    se_scale_measure = se(2);
    se_location_measure = se(3);
    
    population = max_values(:,2); % original population
    n = 1000; % number of bootstrap samples

    if fixedTail == 1
    
        % parameter estimation using bootstrapping
        [GEVparameters] = bootstrap(population, n);
        
        % check normality of bootstrap estimates and compute quantiles for 95%
        [Qt, Qs, Ql] = check_normality(GEVparameters, parmhat, 1, 1);

        % Calculate statistics for each parameter
        mean_shape = mean(GEVparameters.tail);
        std_shape = std(GEVparameters.tail);

        % computation COV
    cov_shape = -(std(GEVparameters.tail)/mean(GEVparameters.tail));
    cov_scale = std(GEVparameters.scale)/mean(GEVparameters.scale);
    cov_location = std(GEVparameters.location)/mean(GEVparameters.location);


        % Comparison of different methods for parameter estimation
        tail_mle = [num2str(tail_measure) ', se=' num2str(se_location_measure)];
        tail_bootstrap = [num2str(mean_shape) ', sd=' num2str(std_shape)];
        scale_mle = [num2str(scale_measure) ', se=' num2str(se_scale_measure)];
        scale_bootstrap = [num2str(mean(GEVparameters.scale)) ', sd=' num2str(std(GEVparameters.scale))];
        location_mle = [num2str(location_measure) ', se=' num2str(se_location_measure)];
        location_bootstrap = [num2str(mean(GEVparameters.location)) ', sd=' num2str(std(GEVparameters.location))];
        
        table = {
            tail_mle, tail_bootstrap;
            scale_mle, scale_bootstrap;
            location_mle, location_bootstrap
            };
        
        MLE_boot = array2table(table, 'VariableNames', {'MLE', 'MLE - Bootstrap'}, 'RowNames',{'tail index', 'scale', 'location'});
        % writetable(MLE_boot, 'pestimation_model.xlsx', 'WriteRowNames', true);
        disp(MLE_boot)
        
        
        % Comparison of different methods to compute confidence interval
        ci_tail_mle = [num2str(parmci(1,1)) ', ' num2str(parmci(2,1))];
        ci_tail_bootstrap = [num2str(Qt(1)) ', ' num2str(Qt(2))];
        ci_scale_mle = [num2str(parmci(1,2)) ', ' num2str(parmci(2,2))];
        ci_scale_bootstrap = [num2str(Qs(1)) ', ' num2str(Qs(2))];
        ci_location_mle = [num2str(parmci(1,3)) ', ' num2str(parmci(2,3))];
        ci_location_bootstrap = [num2str(Ql(1)) ', ' num2str(Ql(2))];
        
        table = {
            ci_tail_mle, ci_tail_bootstrap;
            ci_scale_mle, ci_scale_bootstrap;
            ci_location_mle, ci_location_bootstrap
            };
        
        comp_ci = array2table(table, 'VariableNames', {'MLE', 'MLE - Bootstrap'}, 'RowNames',{'ci - tail index', ' ci - scale', 'ci - location'});
        % writetable(comp_ci, 'confint_model.xlsx', 'WriteRowNames', true);
        disp(comp_ci)

        % disp COV
        disp(['COV - shape: ', num2str(cov_shape)]);
        disp(['COV - scale: ', num2str(cov_scale)]);
        disp(['COV - location: ', num2str(cov_location)]);

    elseif fixedTail == 2

        % parameter estimation using bootstrapping with tail fixed
        [GEVparameters] = bootstrap(population, n, tail);
        
        % check normality of bootstrap estimates and compute quantiles for 95%
        [Qt, Qs, Ql] = check_normality(GEVparameters, parmhat, 1, 1);

        % Calculate statistics for each parameter
        mean_shape = tail;
        std_shape = se_tail;

        % computation COV
        cov_scale = std(GEVparameters.scale)/mean(GEVparameters.scale);
        cov_location = std(GEVparameters.location)/mean(GEVparameters.location);
        
        % disp COV
        disp(['COV - scale: ', num2str(cov_scale)]);
        disp(['COV - location: ', num2str(cov_location)]);

    end

    mean_scale = mean(GEVparameters.scale);
    mean_location = mean(GEVparameters.location);
    std_scale = std(GEVparameters.scale);
    std_location = std(GEVparameters.location);
    
    % 95% confidence 
    tail_confidence = parmci(:,1);
    scale_confidence =  parmci(:,2);
    location_confidence = parmci(:,3);
    
%     % Estimates MLE
%     disp(['MLE GEV - shape: ', num2str(tail_measure)]);
%     disp(['MLE GEV - scale: ', num2str(scale_measure)]);
%     disp(['MLE GEV - location: ', num2str(location_measure)]);
%     % Estimates MLE Bootstrap
%     fprintf('shape parameters - mu: %.4f, sigma: %.4f\n', mean_shape, std_shape);
%     fprintf('scale parameters - mu: %.4f, sigma: %.4f\n', mean_scale, std_scale);
%     fprintf('location parameters - mu: %.4f, sigma: %.4f\n', mean_location, std_location);



%     % Plot PDF of each parameter assuming normal distribution
%     x_shape = linspace(mean_shape - 4*std_shape, mean_shape + 4*std_shape, 100); % Range for shape parameter
%     x_scale = linspace(mean_scale - 4*std_scale, mean_scale + 4*std_scale, 100); % Range for scale parameter
%     x_location = linspace(mean_location - 4*std_location, mean_location + 4*std_location, 100); % Range for location parameter
%     
%     pdf_shape = normpdf(x_shape, mean_shape, std_shape); % PDF of shape parameter
%     pdf_scale = normpdf(x_scale, mean_scale, std_scale); % PDF of scale parameter
%     pdf_location = normpdf(x_location, mean_location, std_location); % PDF of location parameter
%     
%     % shape
%     figure;
%     plot(x_shape, pdf_shape, 'Color',	"#A2142F", 'LineWidth', 1.5);
%     title('PDF of Shape Parameter');
%     xlabel('Shape');
%     ylabel('Density');
%     xlim([min(x_shape), max(x_shape)])
%     hold on;
%     % Draw a vertical line at the mean
%     plot([mean_shape, mean_shape], [0, max(pdf_shape)], 'k--', 'LineWidth', 1);
%     % Shade the region representing ±1 std deviation
%     x_fill = [mean_shape-std_shape, linspace(mean_shape-std_shape, mean_shape+std_shape, 100), mean_shape+std_shape];
%     y_fill = [0, normpdf(linspace(mean_shape-std_shape, mean_shape+std_shape, 100), mean_shape, std_shape), 0];
%     fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     hold off;
%     
%     % scale
%     figure;
%     plot(x_scale, pdf_scale, 'Color',	"#0072BD", 'LineWidth', 1.5);
%     title('PDF of Scale Parameter');
%     xlabel('Scale');
%     ylabel('Density');
%     xlim([min(x_scale), max(x_scale)])
%     hold on;
%     % Draw a vertical line at the mean
%     plot([mean_scale, mean_scale], [0, max(pdf_scale)], 'k--', 'LineWidth', 1);
%     %Shade the region representing ±1 std deviation
%     x_fill = [mean_scale-std_scale, linspace(mean_scale-std_scale, mean_scale+std_scale), mean_scale+std_scale];
%     y_fill = [0, normpdf(linspace(mean_scale-std_scale, mean_scale+std_scale), mean_scale, std_scale), 0];
%     fill(x_fill, y_fill, 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     hold off;
%     
%     % location
%     figure;
%     plot(x_location, pdf_location, 'Color',	"#77AC30", 'LineWidth', 1.5);
%     title('PDF of Location Parameter');
%     xlabel('Location');
%     ylabel('Density');
%     xlim([min(x_location), max(x_location)])
%     hold on;
%     % Draw a vertical line at the mean
%     plot([mean_location, mean_location], [0, max(pdf_location)], 'k--', 'LineWidth', 1);
%     % Shade the region representing ±1 std deviation
%     x_fill = [mean_location-std_location, linspace(mean_location-std_location, mean_location+std_location), mean_location+std_location];
%     y_fill = [0, normpdf(linspace(mean_location-std_location, mean_location+std_location), mean_location, std_location), 0];
%     fill(x_fill, y_fill, 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%     hold off;

end


%% table to save data

if exist("tail_model", "var") && exist("tail_measure", "var") && exist("tail_model_measurements","var")

% Variables
tail_model = {[num2str(tail) ', sd=' num2str(se_tail)]};
tail_measurements = {[num2str(tail_measure) ', sd=' num2str(se_tail_measure)]};
tail_model_measurements = {[num2str(tail) ', sd=' num2str(std_shape)]};

scale_model = {[num2str(scale_model) ', sd=' num2str(se_scale)]};
scale_measurements = {[num2str(scale_measure) ', sd=' num2str(se_scale_measure)]};
scale_model_measurements = {[num2str(mean_scale) ', sd=' num2str(std_scale)]};

location_model = {[num2str(location_model) ', sd=' num2str(se_loc)]};
location_measurements = {[num2str(location_measure) ', sd=' num2str(se_location_measure)]};
location_model_measurements = {[num2str(mean_location) ', sd=' num2str(std_location)]};

% Create the cell array by storing the variables
table = {
    tail_model, tail_measurements, tail_model_measurements;
    scale_model, scale_measurements, scale_model_measurements;
    location_model, location_measurements, location_model_measurements
    };


% Create a table
results = array2table(table, 'VariableNames', {'model', 'measurements', 'model + measurements'}, 'RowNames',{'tail index', 'scale', 'location'});
writetable(results, 'results.txt', 'WriteRowNames', true, 'Delimiter', '\t');
writetable(results, 'results.xlsx', 'WriteRowNames', true);
end


%% return period

% Define the return periods (in years)
return_periods = [10, 50, 100, 1000, 10000];

% Calculate wind speeds associated with return periods
wind_speeds = zeros(size(return_periods));

for i = 1:length(return_periods)
    T = return_periods(i);
    % Calculate the exceedance probability for the return period
    P = 1 / T;
    
    % Inverse CDF of the GEV distribution
    % We need to solve for x in 1 - F(x) = P, which translates to F(x) = 1 - P
    prob = 1 - P;
    
    % Calculate the wind speed corresponding to the probability
    wind_speeds(i) = gevinv(prob, mean(GEVparameters.tail), mean(GEVparameters.scale), mean(GEVparameters.location));
end

% Display the results
disp('Return Period (years)     Wind Speed (m/s)');
disp([return_periods', wind_speeds']);

figure
plot(return_periods, wind_speeds)
grid on
set(gca, 'XScale', 'log')

%% GP - PoT - MLE - Analytical unc


%% GW - PoT - iHill(i) - Analytical unc


%% log-GW - PoT - iHill(i) - Analytical unc

