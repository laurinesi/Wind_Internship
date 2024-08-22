% exploration of Schiphol dataset

% run("smoothed_windspeed_datasets.m")

addpath('..\wind-speeds\')
data_Schiphol = readmatrix("data_Schiphol.txt"); % measurements
model_Schiphol = readmatrix("model_Schiphol.txt"); % model

% Block Maxima = ne prendre que la valeur maximale de chaque année

% Assuming data_Schiphol is an Nx2 matrix where:
% - Column 1 is the year
% - Column 2 is the value (e.g., wind speed)

years = unique(data_Schiphol(:,1)); % Find unique years
max_values = zeros(length(years), 2); % Preallocate for speed

for i = 1:length(years)
    year = years(i);
    % Find the maximum value for the current year
    max_value = max(data_Schiphol(data_Schiphol(:,1) == year, 2));
    max_values(i,:) = [year, max_value]; % Store the year and max value
end

% Display the block maxima for each year
disp(max_values);

% Create a table with appropriate column names
FF_BM = array2table(max_values, 'VariableNames', {'Year', 'FF_BM'});

% Save the table as a CSV file
writetable(FF_BM, 'Schiphol_BM.csv');

% Save the table as a text file
writetable(FF_BM, 'Schiphol_BM.txt', 'Delimiter', '\t');

figure
% plot measurements/years
plot(data_Schiphol(:,1), data_Schiphol(:,2),'Marker','x')
hold on
plot(max_values(:,1),max_values(:,2),'Marker','x','Color','r','LineStyle','none')
yticks(0:4:30)

figure
plot(max_values(:,1),max_values(:,2),'Marker','x','LineStyle','none')
yticks(15:2:27)

