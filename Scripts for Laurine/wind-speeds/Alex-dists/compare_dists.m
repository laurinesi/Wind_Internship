% Compare fitted dists (Schiphol 1h mean) with the new dists of Alex (Cabauw & Schiphol 10min mean)
% for wind speeds

%% settings

C_dist = 1; % [1 2 3] ~ [GEV GP GW]
F_dist = 2; % [1 2] ~ [G GEV]
S_dist = 1; % [1 2 3] ~ [GEV GP GW]

set(groot,'defaultLineLineWidth',0.7)

% Initialization
Fx = [];
V_Cabauw = [];
P_Cabauw = [];
V_Schiphol = [];
P_Schiphol = [];

%% GEV Cabauw
if C_dist == 1
    
    % loading struct df for GEV
    addpath('..\wind-speeds\Alex-dists\')
    load('GEV_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GEV';

    % initialization
    n = size(proba,1);

    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end

    index = 9;
    % add values before the tail (can't fit GEV on those)
    for j = 1:index
        V_Cabauw(j,1) = df.points(j,1);
        P_Cabauw(j,1) = 1 - 1/df.points(j,2);
    end
    P_alldata = reshape([Fx;P_Cabauw],[],1);
    Fx = sort(P_alldata);
    V_alldata = reshape([original;V_Cabauw],[],1);
    original = sort(V_alldata);
    
    A = [original, Fx];
    
    % saving data as a txt file
    writematrix(A,'..\wind-speeds\Alex-dists\V_Cabauw_maxima_tail_Fx_GEV','Delimiter','\t','FileType','text') 

end

%% GP Cabauw
if C_dist == 2

    % loading struct df for GP
    addpath('..\wind-speeds\Alex-dists\')
    load('GP_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GP';

    % initialization
    n = size(proba,1);

    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end

    index = 2100;
    % whole dataset
    for j = 1:index
        V_Cabauw(j,1) = df.points(j,1);
        P_Cabauw(j,1) = df.points(j,2);
    end
    P_alldata = reshape([Fx;P_Cabauw],[],1);
    Fx = sort(P_alldata);
    V_alldata = reshape([original;V_Cabauw],[],1);
    original = sort(V_alldata);
    
    A = [original, Fx];
    
    % saving data as a txt file
    writematrix(A,'..\wind-speeds\Alex-dists\V_Cabauw_maxima_tail_Fx_GP','Delimiter','\t','FileType','text') 

end

%% GW Cabauw
if C_dist == 3

    % loading struct df for GW
    addpath('..\wind-speeds\Alex-dists\')
    load('GW_Cabauw_est.mat');

    proba = df.proba;
    original = df.original;
    CdistType = 'GW';

    % initialization
    n = size(proba,1);

    % new field : Fx = probability P(X <= x)
    % P(X <= x) = 1 - df.proba
    for i = 1:n
        Fx(i,1) = 1 - proba(i,1);
    end

    index = 2100;
    % whole dataset
    for j = 1:index
        V_Cabauw(j,1) = df.points(j,1);
        P_Cabauw(j,1) = df.points(j,2);
    end
    P_alldata = reshape([Fx;P_Cabauw],[],1);
    Fx = sort(P_alldata);
    V_alldata = reshape([original;V_Cabauw],[],1);
    original = sort(V_alldata);
    
    A = [original, Fx];
    
    % saving data as a txt file
    writematrix(A,'..\wind-speeds\Alex-dists\V_Cabauw_maxima_tail_Fx_GW','Delimiter','\t','FileType','text') 
    
end

%% G Schiphol, yearly maxima of hourly mean wind speeds (FH)
if F_dist == 1

    addpath('..\wind-speeds\fitted-dists\');
    %load('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt');
    A_S = dlmread('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt','\t');

    % shift F_T (50years maxima proba) to F_t (1year maxima proba)
    A_S(:,2) = A_S(:,2).^(1/50);

    % 
    SdistType = 'G';

    % saving data as a txt file
    writematrix(A_S,'V_Schiphol_maxima_conf2_Fx_G_MLE_sec1','Delimiter','\t','FileType','text')

end

%% GEV Schiphol, yearly maxima of hourly mean wind speeds (FH)
if F_dist == 2

    addpath('..\wind-speeds\fitted-dists\');
    %load('V_Schiphol_maxima_conf1_Fx_G_MLE_sec1.xpt');
    A_S = dlmread('V_Schiphol_maxima_conf1_Fx_GEV_MLE_sec1.xpt','\t');

    % shift F_T (50years maxima proba) to F_t (1year maxima proba)
    A_S(:,2) = A_S(:,2).^(1/50);
    
    %
    SdistType = 'GEV';

    % saving data as a txt file
    writematrix(A_S,'V_Schiphol_maxima_conf2_Fx_GEV_MLE_sec1','Delimiter','\t','FileType','text')

end

%% GEV Schiphol, 10min mean wind speeds (FF)
% BM MLE
if S_dist == 1

    % loading struct GEV_est
    addpath('..\wind-speeds\Alex-dists\')
    load('GEV_Schiphol_est.mat');
    
    proba_S = GEV_est.proba;
    original_S = GEV_est.original;
    S10distType = 'GEV';

    % initialization
    n = size(proba_S,1);
    Fx_S = zeros(n,1);
    
    for i = 1:n
        Fx_S(i,1) = 1 - proba_S(i,1);
    end
    
    indx = 34;
    % add values before the tail (can't fit GEV on those)
    for j = 1:indx
        V_Schiphol(j,1) = GEV_est.points(j,1);
        P_Schiphol(j,1) = 1 - 1/GEV_est.points(j,2);
    end
    x = reshape([Fx_S;P_Schiphol],[],1);
    Fx_S = sort(x);
    y = reshape([original_S;V_Schiphol],[],1);
    original_S = sort(y);

    % faire la distinction entre fin de courbe fit GEV et debut de courbe
    % cdf avec les valeurs
    
    B = [original_S, Fx_S];
    
    % saving data as a txt file
    writematrix(B,'..\wind-speeds\Alex-dists\V_Schiphol_maxima_tail_Fx_GEV','Delimiter','\t','FileType','text') 

end

%% GP Schiphol 10min
% POT MLE
if S_dist == 2

    % loading struct GP_est
    addpath('..\wind-speeds\Alex-dists\')
    load('GP_Schiphol_est.mat');

    proba_S = GP_est.proba;
    original_S = GP_est.original;
    S10distType = 'GP';
    
    % initialization
    n = size(proba_S,1);
    Fx_S = zeros(n,1);
    
    for i = 1:n
        Fx_S(i,1) = 1 - proba_S(i,1);
    end

    indx = 7348;
    % whole dataset
    for j = 1:indx
        V_Schiphol(j,1) = GP_est.points(j,1);
        P_Schiphol(j,1) = GP_est.points(j,2);
    end
    x = reshape([Fx_S;P_Schiphol],[],1);
    Fx_S = sort(x);
    y = reshape([original_S;V_Schiphol],[],1);
    original_S = sort(y);
    
    B = [original_S, Fx_S];
    
    % saving data as a txt file
    writematrix(B,'..\wind-speeds\Alex-dists\V_Schiphol_maxima_tail_Fx_GP','Delimiter','\t','FileType','text') 

end

%% GW Schiphol 10min
% POT iHilli
if S_dist == 3

    % loading struct GW_est
    addpath('..\wind-speeds\Alex-dists\')
    load('GW_Schiphol_est.mat');

    proba_S = GW_est.proba;
    original_S = GW_est.original;
    S10distType = 'GW';
    
    % initialization
    n = size(proba_S,1);
    Fx_S = zeros(n,1);
    
    for i = 1:n
        Fx_S(i,1) = 1 - proba_S(i,1);
    end

    indx = 7353;
    % whole dataset
    for j = 1:indx
        V_Schiphol(j,1) = GW_est.points(j,1);
        P_Schiphol(j,1) = GW_est.points(j,2);
    end
    x = reshape([Fx_S;P_Schiphol],[],1);
    Fx_S = sort(x);
    y = reshape([original_S;V_Schiphol],[],1);
    original_S = sort(y);
    
    B = [original_S, Fx_S];
    
    % saving data as a txt file
    writematrix(B,'..\wind-speeds\Alex-dists\V_Schiphol_maxima_tail_Fx_GW','Delimiter','\t','FileType','text') 

end


%% CDFs

C_name = append('Cabauw', ' ', CdistType);
S_name = append('Schiphol 1h', ' ', SdistType);
S10_name = append('Schiphol', ' ', S10distType);

figure
% Plot Cabauw tail
plot(A(index+1:end,1), A(index+1:end,2),'Color','r')
hold on
% Plot Cabauw .points
plot(A(1:index+1,1), A(1:index+1,2), '--','color','r')

% Plot Schiphol 1h
plot(A_S(:,1), A_S(:,2),'color','k')  
hold on

% Plot Schiphol 10min tail
plot(B(indx+1:end,1),B(indx+1:end,2),'color','b')
hold on
% Plot Schiphol 10min .points
plot(B(1:indx+1,1),B(1:indx+1,2), '--','color','b')

% Set x-axis to logarithmic scale
set(gca, 'XScale', 'log');

xlabel('Wind Speed (m/s)')
ylabel('P(X <= x)')                      
grid minor
legend(C_name, 'Cabauw data', S_name, S10_name, 'Schiphol data', Location='northwest')
legend box off
title('CDFs of whole datasets')
xlim([0, 30]); % Limiting the x-axis to 30 m/s


%% CDFs (gumbel plot)

C_name = append('Cabauw', ' ', CdistType);
S_name = append('Schiphol 1h', ' ', SdistType);
S10_name = append('Schiphol', ' ', S10distType);

figure
% Plot Cabauw tail
plot(log(A(index+1:end,1)), -log(-log(A(index+1:end,2))),'Color','r')
hold on
% Plot Cabauw .points
plot(log(A(1:index+1,1)), -log(-log(A(1:index+1,2))), '--','color','r')

% Plot Schiphol 1h
plot(log(A_S(:,1)), -log(-log(A_S(:,2))),'Color','k')  
hold on

% Plot Schiphol 10min tail
plot(log(B(indx+1:end,1)),-log(-log(B(indx+1:end,2))),'Color','b')
hold on
% Plot Schiphol 10min .points
plot(log(B(1:indx+1,1)),-log(-log(B(1:indx+1,2))), '--','color','b')

% Set x-axis to logarithmic scale
set(gca, 'XScale', 'log');

xlabel('Wind Speed (m/s)')
ylabel('-log(-log(P(X <= x))')
grid minor
legend(C_name, 'Cabauw data', S_name, S10_name, 'Schiphol data', Location='northwest')
legend box off
title('Gumbel plot')
ylim([-2, 13]);

