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
data(:,1) = data_Cabauw.DateTime; % DateTime
data(:,2) = data_Cabauw.Year; % Year (starting the new year on july 1)
data(:,3) = data_Cabauw.F010; % F010 : 10-m Wind Speed (m/s)

RACMO(:,1) = RACMO_Cabauw.DateTime; % DateTime
RACMO(:,2) = RACMO_Cabauw.Year; % Year
RACMO(:,3) = RACMO_Cabauw.w10m; % w10m : 10-m Wind Speed (m/s)

KNW(:,1) = KNW_Cabauw.DateTime; % DateTime
KNW(:,2) = KNW_Cabauw.Year; % Year
KNW(:,3) = KNW_Cabauw.F010; % F010 : Wind Speed at 10m Height


% saving data as a txt file
data_Cabauw = array2table(data, "VariableNames",{'DateTime', 'Year', 'F010'});
writetable(data_Cabauw,'measure_Cabauw','Delimiter','\t','FileType','text')

RACMO_Cabauw = array2table(RACMO, "VariableNames",{'DateTime', 'Year', 'w10m'});
writetable(RACMO_Cabauw,'RACMO_Cabauw','Delimiter','\t','FileType','text')

KNW_Cabauw = array2table(KNW, "VariableNames",{'DateTime', 'Year', 'F010'});
writetable(KNW_Cabauw,'KNW_Cabauw','Delimiter','\t','FileType','text')