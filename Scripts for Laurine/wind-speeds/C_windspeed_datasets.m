% Review windspeed datasets at Cabauw

addpath('..\wind-speeds\datasets\')

data_Cabauw = readtable("Cabauw_measure.csv"); % measurements
data_RACMO_Cabauw = readtable("fulldata.csv"); % RACMO
data_KNW_Cabauw = readtable("modeldata.csv"); % KNW

%% Initialization
data_10 = table;
RACMO_10 = table;
KNW_10 = table;


% Definition of tables
% Measurements
year_column = data_Cabauw.Year;
filtered_idx = year(year_column) > 2000 & year(year_column) < 2020;  % Measurements (only 2000-2020)
measurements_Cabauw = data_Cabauw(filtered_idx,:);
measurements_Cabauw.Year = year(measurements_Cabauw.Year);
data_10.Year = measurements_Cabauw.Year;
data_10.F010 = measurements_Cabauw.F010;


% RACMO
year_columnR = data_RACMO_Cabauw.Year;
filtered_idxR = year_columnR > 974 & year_columnR < 2014;  % Measurements (only 2000-2020)
filtered_RACMO_Cabauw = data_RACMO_Cabauw(filtered_idxR,:);
RACMO_Cabauw = movevars(filtered_RACMO_Cabauw,"Year","Before","fh050");
RACMO_10.Year = RACMO_Cabauw.Year; % Year
RACMO_10.w10m = RACMO_Cabauw.w10m; % w10m : 10-m Wind Speed (m/s)


% KNW
year_column_KNW = data_KNW_Cabauw.Year;
year_values = year(year_column_KNW);
KNW_Cabauw = data_KNW_Cabauw;
KNW_Cabauw.Year = year_values;
KNW_10.Year = KNW_Cabauw.Year; % Year
KNW_10.F010 = KNW_Cabauw.F010; % F010 : Wind Speed at 10m Height


% saving data as a txt file : [years data_at_10m]
writetable(data_10,'measure_Cabauw','Delimiter','\t','FileType','text')
writetable(RACMO_10,'RACMO_Cabauw','Delimiter','\t','FileType','text')
writetable(KNW_10,'KNW_Cabauw','Delimiter','\t','FileType','text')

clearvars -except measurements_Cabauw RACMO_Cabauw KNW_Cabauw
