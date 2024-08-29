% probabilistic description of wind loads

addpath('..\wind-speeds\tools\')

% Datasets
run("S_windspeed_datasets.m") % Schiphol measure & model (RACMO)
% run("C_windspeed_datasets.m") % Cabauw measure & model (RACMO, KNW)

%% GEV - BM - MLE - Bootstrap

% tail index (shape) : GEV MLE on weather model dataset (BM)

dataset = model_Schiphol;
[max_values] = BM_select(dataset);

parmhat = gevfit(max_values(:,2));
tail = parmhat(1);

disp(['MLE GEV - tail index: ', num2str(tail)]);

% scale & location : bootstrap on measurements dataset (BM)

dataset = data_Schiphol;
[max_values] = BM_select(dataset);

population = max_values(:,2); % original population
n = 10; % number of bootstrap samples

% need to fix the tail in gevfit for bootstrap
[GEVparameters] = bootstrap(population, n);

%% GP - PoT - MLE - Analytical unc


%% GW - PoT - iHill(i) - Analytical unc


%% log-GW - PoT - iHill(i) - Analytical unc