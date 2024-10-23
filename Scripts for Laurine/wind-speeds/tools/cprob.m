clear all;
clc;

%% Transformation of Vb from Eurocode (50-years extreme) to 20-years extreme

addpath('..\..\wind-speeds\')
table_data_Cabauw = readtable("measure_Cabauw.txt");
% run("C_windspeed_datasets.m") % to create the txt file 'measure_Cabauw.txt'

% get the yearly maximum wind speeds
dataset = table_data_Cabauw;
[max_values] = BM_select(dataset);


% Compute the cumulative probability for these yearly maxima: F_N_max() function
[F_N_max] = F_N_max(max_values(:,2));

% Take the cumulative distribution value for the maximum wind speed and fill in for p
[max_value, max_index] = max(max_values(:,2));
p = F_N_max(max_index,1);

vb_0 = 24.5;
z = 10;
z0 = 0.132;
z0_ref = 0.05;
% Compute roughness factor
cr = 0.19*(z0/z0_ref)^0.07 * log(z/z0);

K = 0.281;
n = 0.5;
A = 1-K*log(-log(1-p));
B = 1-K*log(-log(0.98));
c_prob = (A/B)^n;

vb_20y = vb_0*c_prob*cr(1);

% Wind speed you obtain with cprob*vb*cr can be compared to the maximum yearly maximum in the data
disp(['Maximum value: ', num2str(max_value)]);
disp(['Basic wind speed (20-years extreme): ', num2str(vb_20y)]);

