%% plot return values at different heights for each model and possibily interpolate return values(heights wanted)
% normal if RACMO plot is shifted because scale and location are biased

load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_10m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_20m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_40m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_80m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_140m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_200m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_model_10m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_model_100m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_model_150m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_model_200m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_model_250m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_model_300m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_model_10m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_model_20m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_model_40m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_model_80m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_model_140m.mat')
load('..\..\..\Alex''s code\Alex-dists\GW_Cabauw_measure_model_200m.mat')

height_measure = [10,20,40,80,140,200];
height_model = [10,100,150,200,250,300];

return_values_measure(:,1) = GW_10.original;
return_values_measure(:,2) = GW_20.original;
return_values_measure(:,3) = GW_40.original;
return_values_measure(:,4) = GW_80.original;
return_values_measure(:,5) = GW_140.original;
return_values_measure(:,6) = GW_200.original;

return_values_model(:,1) = GW_m10.model_est;
return_values_model(:,2) = GW_m100.model_est;
return_values_model(:,3) = GW_m150.model_est;
return_values_model(:,4) = GW_m200.model_est;
return_values_model(:,5) = GW_m250.model_est;
return_values_model(:,6) = GW_m300.model_est;

return_values_measure_model(:,1) = GW_mm10.model_est(17);
return_values_measure_model(:,2) = GW_mm20.model_est(17);
return_values_measure_model(:,3) = GW_mm40.model_est(17);
return_values_measure_model(:,4) = GW_mm80.model_est(17);
return_values_measure_model(:,5) = GW_mm140.model_est(17);
return_values_measure_model(:,6) = GW_mm200.model_est(17);
%%
% pour Curve Fitter App
data(:,1) = return_values_measure(1,:);
data_model(:,1) = return_values_model(1,:);
data_measure_model(:,1) = return_values_measure_model(1,:);


%% Evolution of mean wind speed vb -> vm (Eurocodes)

% Parameters
z = [10, 20, 40, 50, 80, 100, 140, 150, 200, 250, 300];

% TERRAIN CATEGORY: roughness lengh for different terrain category
z0_0 = 0.003;
zi_0 = 0.01;
zii_0 = 0.05;
ziii_0 = 0.3;
ziv_0 = 1;

z0 = 0.132; % roughness lenght of current study location
z0_ref = 0.05; % roughness lenght of vb,0

% WIND ZONES: basic wind velocity based on Eurocodes
vbi_0 = 29.5;
vbii_0 = 27;        % Schiphol
vbiii_0 = 24.5;     % Cabauw

% Compute roughness factor
cr = 0.19*(z0/z0_ref)^0.07 * log(z/z0);

% Orography factor
co = 1;

% Compute mean wind speed
vm1 = cr*co*vbi_0;
vm2 = cr*co*vbii_0;
vm3 = cr*co*vbiii_0;

% % Plot Vb/Height
% figure
% % plot vm based on Eurocode basic wind speed
% plot(vm1,z,LineWidth=1)
% hold on
% plot(vm2,z,LineWidth=1)
% hold on
% plot(vm3,z,LineWidth=1)
% hold on
% plot(vm1,z,'LineStyle','none','Marker','o', 'MarkerSize',4, 'MarkerFaceColor','w','MarkerEdgeColor','w')
% hold on
% plot(vm2,z,'LineStyle','none','Marker','o', 'MarkerSize',4, 'MarkerFaceColor','w','MarkerEdgeColor','w')
% hold on
% plot(vm3,z,'LineStyle','none','Marker','o', 'MarkerSize',4, 'MarkerFaceColor','w','MarkerEdgeColor','w')
% title(['Eurocodes - 50-year extreme, r = ', num2str(z0)]); 
% xlabel('Mean Wind Velocity (m/s)');
% ylabel('Height (m)');
% legend('Wind Zone I','Wind Zone II','Wind Zone III', 'Location','northwest')
% legend box off
% grid minor

%% Fit exponential curve
% Sample data (replace with your actual data)
x = return_values_measure_model(1,:);
y = height_measure;

% Fit the data using an exponential model
f = fit(x', y', 'exp1');  % 'exp1' refers to y = a * exp(b * x)

% Display the fitting parameters
disp(f);

% Generate x values, including values going down to 0
extended_x = linspace(0, min(x), 100);  % Generate 100 points from x = 0 to max(x)

% Evaluate the fitted model for the extended x range
extended_y = f(extended_x);

% Plot the original data and the extended fit
figure;
plot(x, y, 'o', 'DisplayName', 'Data');  % Original data points
hold on;
plot(extended_x, extended_y, '-r', 'DisplayName', 'Exponential Fit');  % Extended fit curve
hold off;

% Labels and legend
set(gca, 'YScale', 'log')
title('Exponential Curve Fit with Extrapolation to x = 0');
xlabel('x');
ylabel('y');
legend show;
grid on;

%%
% Sample data (replace with your actual data)
x = vm3;
y = z;

% Fit the data using an exponential model
f = fit(x', y', 'exp1');  % 'exp1' refers to y = a * exp(b * x)

% Display the fitting parameters
disp(f);

% Generate x values, including values going down to 0
extended_xE = linspace(0, min(x), 100);  % Generate 100 points from x = 0 to max(x)

% Evaluate the fitted model for the extended x range
extended_yE = f(extended_xE);

% Plot the original data and the extended fit
figure;
plot(x, y, 'o', 'DisplayName', 'Data');  % Original data points
hold on;
plot(extended_xE, extended_yE, '-r', 'DisplayName', 'Exponential Fit');  % Extended fit curve
hold off;

% Labels and legend
set(gca, 'YScale', 'log')
title('Exponential Curve Fit with Extrapolation to x = 0');
xlabel('x');
ylabel('y');
legend show;
grid on;


%% Find the roughness based on log wind profile
figure
% plot(return_values_measure(1,:),height_measure,'r--',LineWidth=1,DisplayName='V_m measurements')
% hold on
% plot(return_values_model(1,:),height_model,'b--',LineWidth=1,DisplayName='V_m model')
% hold on
plot(return_values_measure_model(1,:),height_measure,'g-',LineWidth=1,DisplayName='V_m measurements + model')
hold on
plot(vm3,z,'k-',LineWidth=1, DisplayName=' V_m Eurocode')
% hold on
% plot(extended_x, extended_y, 'k--', 'DisplayName', 'Exponential Fit');
% hold on
% plot(extended_xE, extended_yE, 'k--', 'DisplayName', 'Exponential Fit');
% xline(34,'k--',LineWidth=0.7)
% xline(38,'k--',LineWidth=0.7)
% yline(60,'k--',LineWidth=0.7)
xlabel('Mean Wind Velocity (m/s)')
ylabel('Height (m)')
title(['r = ', num2str(z0)])
xlim([20 36])
ylim([10 200])
set(gca, 'YScale', 'log')
% legend off
legend Location northwest
grid minor

