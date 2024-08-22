% Review windspeed datasets at Cabauw

addpath('..\wind-speeds\datasets\')

% % Check if the matrices have already been loaded
% if exist("data_Cabauw","var") == 0
%     data_Cabauw = readmatrix("Cabauw_measure.csv"); % measurements
% end
% 
% if exist("RACMO_Cabauw","var") == 0
%     RACMO_Cabauw = readmatrix("fulldata.csv"); % RACMO
% end
% 
% if exist("KNW_Cabauw","var") == 0
%     KNW_Cabauw = readmatrix("modeldata.csv"); % KNW
% end

data_Cabauw = readmatrix("Cabauw_measure.csv"); % measurements
RACMO_Cabauw = readmatrix("fulldata.csv"); % RACMO
KNW_Cabauw = readmatrix("modeldata.csv"); % KNW


% Initialization
data = [];
RACMO = [];
KNW = [];


% Definition of tables
data(:,1) = data_Cabauw(:,1); % DateTime
data(:,2) = data_Cabauw(:,2); % Year (starting the new year on july 1)
data(:,3) = data_Cabauw(:,3); % F010 : 10-m Wind Speed (m/s)

RACMO(:,1) = RACMO_Cabauw(:,1); % DateTime
RACMO(:,2) = RACMO_Cabauw(:,3); % Year
RACMO(:,3) = RACMO_Cabauw(:,19); % w10m : 10-m Wind Speed (m/s)

KNW(:,1) = KNW_Cabauw(:,1); % DateTime
KNW(:,2) = KNW_Cabauw(:,2); % Year
KNW(:,3) = KNW_Cabauw(:,3); % F010 : Wind Speed at 10m Height


% saving data as a txt file
data_Cabauw = array2table(RACMO, "VariableNames",{'DateTime', 'Year', 'F010'});
writetable(data,'measure_Cabauw','Delimiter','\t','FileType','text')

RACMO_Cabauw = array2table(RACMO, "VariableNames",{'DateTime', 'Year', 'w10m'});
writetable(RACMO,'RACMO_Cabauw','Delimiter','\t','FileType','text')

KNW_Cabauw = array2table(RACMO, "VariableNames",{'DateTime', 'Year', 'F010'});
writetable(KNW,'KNW_Cabauw','Delimiter','\t','FileType','text')

