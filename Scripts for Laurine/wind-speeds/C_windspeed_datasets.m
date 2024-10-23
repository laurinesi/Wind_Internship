% Review windspeed datasets at Cabauw

addpath('..\wind-speeds\datasets\')

data_Cabauw = readtable("Cabauw_measure.csv"); % measurements
RACMO_Cabauw = readtable("fulldata.csv"); % RACMO
KNW_Cabauw = readtable("modeldata.csv"); % KNW

%% Initialization
data = [];
RACMO = [];
KNW = [];


% Definition of tables
datetime_column = data_Cabauw.Year;
year_values = year(datetime_column);
mask = (year_values > 2000) & (year_values < 2020);
filtered_years = year_values(mask);
data(:,1) = filtered_years; % Year (starting the new year on july 1)
filtered_data = data_Cabauw.F010(mask);
data(:,2) = filtered_data; % F010 : 10-m Wind Speed (m/s)

for i = 1:size(RACMO_Cabauw.Year)
    RACMO(i,1) = RACMO_Cabauw.Year(i) + 1000; % Year
end
RACMO(:,2) = RACMO_Cabauw.w10m; % w10m : 10-m Wind Speed (m/s)

datetime_column = KNW_Cabauw.Year;
year_values = year(datetime_column);
KNW(:, 1) = year_values; % Year
KNW(:,2) = KNW_Cabauw.F010; % F010 : Wind Speed at 10m Height


% saving data as a txt file
table_data_Cabauw = array2table(data, "VariableNames",{'Year', 'F010'});
writetable(table_data_Cabauw,'measure_Cabauw','Delimiter','\t','FileType','text')

table_RACMO_Cabauw = array2table(RACMO, "VariableNames",{'Year', 'w10m'});
writetable(table_RACMO_Cabauw,'RACMO_Cabauw','Delimiter','\t','FileType','text')

table_KNW_Cabauw = array2table(KNW, "VariableNames",{'Year', 'F010'});
writetable(table_KNW_Cabauw,'KNW_Cabauw','Delimiter','\t','FileType','text')