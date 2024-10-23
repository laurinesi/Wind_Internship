% Monte Carlo method to compute failure probability

%% Random continuous distributions

% Parameters for resistance R and load S
mu_R = 1200;
sigma_R = 180;
% mu_S = 8;
% sigma_S = 1.5;

% Normal distributions
R = makedist('Normal', 'mu', mu_R, 'sigma', sigma_R);
% S = makedist('Normal', 'mu', mu_S, 'sigma', sigma_S);

%% Return values

location = 'Cabauw';

% Load Alex data for GW dist
if location == "Cabauw"
    load('../../wind-speeds/Alex-dists/GW_Cabauw_est.mat');
elseif location == "Schiphol"
    load('../../wind-speeds/Alex-dists/GW_Schiphol_est.mat');
end

return_period = df.return;  % Return period
original = df.original;  % Return values measurements only
model_est = df.model_est;  % Return values RACMO tail
proba = df.proba;  % Probability of exceedance

% Find the index where the return period equals 50
index = find(return_period == 50);

% % Probability of exceedance/wind velocity
% figure
% plot(model_est,proba)
% xlabel('Wind speed (m/s)')
% ylabel('P(X > x)')

p = 1 - proba;
% CDF of derived ws with RACMO tail
figure
plot(model_est, p)
xlabel('Wind speed (m/s)')
ylabel('P(X <= x)')

% Plot return periods/return values 
figure
plot(return_period,model_est, 'b-', LineWidth=1.5)
hold on
plot(return_period, original, 'r-', LineWidth=1.5)
hold on
plot(df.points(:,4), df.points(:,1), Marker="+", LineStyle="none", Color='k')
hold on 
xline(50, 'k--')
xlabel('Return period (year)')
ylabel('Wind speed (m/s)')
legend('RACMO tail','Measurements tail','Measurements observed', 'Location','northwest')
legend Box off
grid on
set(gca, 'XScale', 'log')

%% Computation of wind load using Alex results

% ce = 1;
% cf = 1;
% cs_d = 0.85;
% A = 1;

rho = 1.25;
vb = model_est;

% basic pressure
qb = 0.5*rho*vb.^2;

% wind load
S = qb; % *ce*cf*cs_d*A;

mu_S = mean(S);
sigma_S = std(S);

S = makedist('Normal', 'mu', mu_S, 'sigma', sigma_S);
% x_pdf = linspace(min(mu_S - 3*sigma_S), max(mu_S + 3*sigma_S), 100);
% pdf_S = pdf(S, x_pdf);
% figure
% plot(x_pdf,pdf_S, 'k-')

%% Plot PDF

x_pdf = linspace(min(mu_R - 3*sigma_R, mu_S - 3*sigma_S), max(mu_R + 3*sigma_R, mu_S + 3*sigma_S), 100);
pdf_R = pdf(R, x_pdf);
pdf_S = pdf(S, x_pdf);

figure;
plot(x_pdf, pdf_R, 'b-', 'LineWidth', 2); % PDF of R
hold on;
plot(x_pdf, pdf_S, 'g-', 'LineWidth', 2); % PDF of S
xlabel('Value');
ylabel('PDF');
title('PDF of R and S');
legend('R (Resistance)', 'S (Load)', 'Location', 'Best');
grid minor;
hold off;

%% Generate grid
x = linspace(mu_R - 3*sigma_R, mu_R + 3*sigma_R, 100);
y = linspace(mu_S - 3*sigma_S, mu_S + 3*sigma_S, 100);
[X, Y] = meshgrid(x, y);

% Compute pdf
f_R = pdf(R, X);
f_S = pdf(S, Y);

% Joint pdf
joint_pdf = f_R .* f_S;

%% Plot 2D
figure;
contour(X, Y, joint_pdf, 'LineWidth', 1.5);
hold on;
plot(x, x, 'r--', 'LineWidth', 2);  % Line between safe zone and failure zone
xlabel('R (Resistance)');
ylabel('S (Load)');
title('Joint pdf (2D)');
grid on;

%% Plot 3D
figure;
surf(X, Y, joint_pdf);
hold on;
plot3(x, x, zeros(size(x)), 'r--', 'LineWidth', 2);
xlabel('R (Resistance)');
ylabel('S (Load)');
zlabel('Density');
title('Joint pdf (3D)');
shading interp;
grid on;

%% Computation of failure probability (P(R < S))

% analytical method
integrand = @(r, s) pdf(R, r) .* pdf(S, s) .* (r < s);
failure_prob_analytical = integral2(integrand, -Inf, Inf, -Inf, Inf);

fprintf('Failure probability (analytical) : %.4f\n', failure_prob_analytical);

%% Monte Carlo

N = 1e4;  % number of simulations
R_samples = random(R, N, 1);  % Generate samples of R
S_samples = random(S, N, 1);  % Generate samples of S

% Count number of times where R < S
failures = sum(R_samples < S_samples);

% Computation of failure probability
failure_prob_mc = failures / N;

fprintf('Failure probability (Monte Carlo) : %.4f\n', failure_prob_mc);

%% Reliability index

mu_Z = mu_R - mu_S;
sigma_Z_squarred = sigma_R^2 + sigma_S^2; % if R and s independent (otherwise : sigma_R^2 + sigma_S^2 - 2*correlation_RS*sigma_R*sigma_S
sigma_Z = sqrt(sigma_Z_squarred);

beta_target = 3.8;

beta = mu_Z / sigma_Z;

fprintf('Reliability index : %.4f\n', beta);
