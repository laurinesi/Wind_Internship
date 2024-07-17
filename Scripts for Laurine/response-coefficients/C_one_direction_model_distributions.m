%% explanation script
% ** Derive model parameters
% ** Plot distribution functions for both chosen block duration and different block durations
% ** Derive fractile values 

% 1: Run C_one_direction_processing. The first 250 datapoints are skipped
% because of irregular behaviour in FEM response. 
% 2: Change settings and preferences
% 3: Saves data automatically if required

%% settings and preferences
% calculation settings - changable 
          extremes = 1; % [1 2] ~ [moment force] 
          sections = 1;
          bootstrapping = 1; % [1 2] ~ [no yes] ! Change b accordingly, if no bootstrapping b = 1 !
          Nreps = 10; % bootstrap repetitions
          n = 100000; % number of array elements in Fx
          % fracile values
          P_cook = 0.78;
          P_code = 0.95;
          P_fractiles_hour = [P_cook; P_code];
          P_fractiles = P_cook;
      
% plot settings - changable
          plot_set = 3;  % [1 2 3] ~ [normaldomain logdomain gumbeldomain]
          plot_tT_extremes = 1; % [1 2 3] ~ [t-yearly-extremes  T-yearly extremes  t+T-yearly extremes]
          plotFs = 20; % (Chosen) block duration for which to plot fitted distribution function
          
          z = 20; % x-domain, for xpt files 20 is chosen   % 4 is chosen for model distribution figures
          b = 1; % bootstrap-data which are plotted
          fontsize = 16;
% saving settins: [0 1] ~ [no-save save] 
          save_Fx_G_MLE = 1;
          save_Fx_GEV_MLE = 1;
          plot_scatter = 1; %[0 means no 1 means yes] Plot original sample data in distribution plots
          plot_figures_model = 1;
                save_figures_model = 0; % should the plots be saved?
          plotblocks = 0; % [plot fitted distributions for both shifted t-extremes as original hourly extremes for T_fs]
                save_plotblocks = 0;
          if plotblocks == 1
              extremes = 3; % extremes of static moment coefficients
          else
          end
          plot_pdf = 0;
          moments_and_forces = 0; % 0 = not run both at least once; 1 = both for moments and forces the script has been run at least once. 
          saving_directory = 'figures';
% plot settings - layout
          markers = {'d','s','^','o','s', 'p', 'h'};  
          colors = {[0 0 1],[1 0 0],[0 0 0],[0 1 0],[1 0 1],[1 1 0]};
          colors2 = {[0.8 0.8 0.8],[0.8 0.8 0.8],[0.6 0.6 0.6],[1 1 1]};

          % required ticks
          if plot_set == 1
                ytick = [0, 0.1, 0.2 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
          elseif plot_set == 2
                ytick = [0.001 0.01, 0.1,0.2, 0.4, 0.6, 0.8, 1];
          elseif plot_set == 3
                %ytick = [0.000001^10, 0.5, 0.78, 0.9, 0.95, 0.99, 0.99999, 1];
                ytick = [0.000001^10, 0.5, 0.9, 0.99, 0.99999, 1];
          end



%%  RUN SCRIPT 
addpath('..\general\')
run('initialize_variables') 


%% simplify analysis data

        
for T = T_fs
            X_t_analysis_mom{T}(1,:) = Cm_max{T};
            X_t_analysis_for{T}(1,:) = Cf_max{T};
            X_t_analysis_stat{T}(1,:) = Cm_stat{T};
 
end
        
%% apply bootstrap settings
        if bootstrapping == 1        
                % do nothing  
        else if bootstrapping == 2  % nonparametric bootstrap - resampling            
                           
                for i = T_fs
                    id{i} = sort(ceil(rand(Nreps,N{i})*N{i}),2); 
                    id_stat{i} = sort(ceil(rand(Nreps,N_stat{i})*N_stat{i}),2); 
                    X_t_analysis_mom{i} = [sort(X_t_analysis_mom{i}(1,:)); Cm_max{i}(id{i})];    % first row original data 
                    X_t_analysis_for{i} = [sort(X_t_analysis_for{i}(1,:)); Cf_max{i}(id{i})];    % first row original data                            
                    X_t_analysis_stat{i} = [sort(X_t_analysis_stat{i}(1,:)); Cm_stat{i}(id_stat{i})];    % first row original data                            

                end
            end        
        end

clearvars id id_stat

%%  Derivation model parameters
if extremes == 1      
    for i =  T_fs          

        for j = 1 : size(X_t_analysis_mom{i},1)  % for all bootstrapped samples, might also be just 1            
            % current sample data
                    X_t_analysis_mom{i}(j,:)=sort(X_t_analysis_mom{i}(j,:));  %  X = current (bootstrapped) sample
            % moments
                    [mean_X_mom(j,i) , std_X_mom(j,i), a_X_mom(j,i) ] = moments_sample_data(X_t_analysis_mom{i}(j,:));
            % empirical distribution function T_fs
                    F_N_Tfs_mom{i}(j,:) = F_N_max(X_t_analysis_mom{i}(j,:)); 
            % empirical distribution function T3600
                    F_N_T3600_mom{i}(j,:) = F_N_Tfs_mom{i}(j,:).^(T_3600/i);
            % model parameters
                    [alpha1_MLE_mom(j,i) , u1_MLE_mom(j,i), nlogL_G_MLE_mom(1,i)] = MLE_Gumbel_max(X_t_analysis_mom{i}(j,:)); % Gumbel-MLE
                    parGEVmom{j}(:,i) = gevfit(X_t_analysis_mom{i}(j,:));
            % settings analytic distribution functions Tfs
                    xmin_max = min(X_t_analysis_mom{i}(j,:)); 
                    xmax_max = max(X_t_analysis_mom{i}(j,:))+z*abs(std_X_mom(j,i));
                    x_mom{i}(j,:) = linspace(xmin_max,xmax_max, n);
            % analytical distribution function Tfs
                     Fx_GEV_MLE_Tfs_mom{i}(j,:) = gevcdf(x_mom{i}(j,:),parGEVmom{j}(1,i),parGEVmom{j}(2,i),parGEVmom{j}(3,i));
                     Fx_G_MLE_Tfs_mom{i}(j,:) = cdf_Gumbel_max(alpha1_MLE_mom(j,i) , u1_MLE_mom(j,i), x_mom{i}(j,:)); % Gumbel-MLE
                             
            % analytical distribution function T3600
                    Fx_GEV_MLE_T3600_mom{i}(j,:) = Fx_GEV_MLE_Tfs_mom{i}(j,:).^(T_3600/i);
                    Fx_G_MLE_T3600_mom{i}(j,:) = Fx_G_MLE_Tfs_mom{i}(j,:).^(T_3600/i);

              % DESIGN FRACTILES 10S EXTREMES   
                    x_fractiles_GEV_MLE_mom_Tfs{i}(j,:) = fractile_values_empirical_Fx(x_mom{i}(j,:),Fx_GEV_MLE_Tfs_mom{i}(j,:), P_fractiles.^(plotFs/T_3600));
                    x_fractiles_G_MLE_mom_Tfs{i}(j,:) = fractile_values_empirical_Fx(x_mom{i}(j,:),Fx_G_MLE_Tfs_mom{i}(j,:), P_fractiles.^(plotFs/T_3600));

              % DESIGN FRACTILES HOURLY EXTREMES 
                   x_fractiles_GEV_MLE_mom{i}(j,:) = fractile_values_empirical_Fx(x_mom{i}(j,:),Fx_GEV_MLE_T3600_mom{i}(j,:), P_fractiles);
                   cdfG = interp1(x_mom{i}(j,:)',Fx_GEV_MLE_T3600_mom{i}(j,:)',0.97);
                   x_fractiles_G_MLE_mom{i}(j,:) = fractile_values_empirical_Fx(x_mom{i}(j,:),Fx_G_MLE_T3600_mom{i}(j,:), P_fractiles);

              % save distribution functions

                    if save_Fx_G_MLE == 1
                    str= [ 'C_Refbuilding_mom_Fx_G_MLE_' num2str(i) ' = transpose([x_mom{i}(1,:); Fx_G_MLE_T3600_mom{i}(1,:)]);']; 
                    eval(str);  
                    str= [ 'dlmwrite(' char(39) 'fitted-dists\C_Refbuilding_mom_Fx_G_MLE_' num2str(i) '.xpt' char(39) ',C_Refbuilding_mom_Fx_G_MLE_' num2str(i) ',' char(39) 'delimiter' char(39) ',' char(39) '\t' char(39) ',' char(39) 'precision' char(39) ',16)'  ];  
                    eval(str)            
                    end    
                    if save_Fx_GEV_MLE == 1
                    str= [ 'C_Refbuilding_mom_Fx_GEV_MLE_' num2str(i) ' = transpose([x_mom{i}(1,:); Fx_GEV_MLE_T3600_mom{i}(1,:)]);']; 
                    eval(str);  
                    str= [ 'dlmwrite(' char(39) 'fitted-dists\C_Refbuilding_mom_Fx_GEV_MLE_' num2str(i) '.xpt' char(39) ',C_Refbuilding_mom_Fx_GEV_MLE_' num2str(i) ',' char(39) 'delimiter' char(39) ',' char(39) '\t' char(39) ',' char(39) 'precision' char(39) ',16)'  ];  
                    eval(str)         
                    end   
        end
    end

  
elseif extremes == 2    
        for i =  T_fs          
        F_N_Tfs_for{i} = [];
        x_for{i} = [];
        for j = 1 : size(X_t_analysis_for{i},1)  % for all bootstrapped samples, might also be just 1            
            % current sample data
                    X_t_analysis_for{i}(j,:)=sort(X_t_analysis_for{i}(j,:));  %  X = current (bootstrapped) sample
            % moments
                    [mean_X_for(j,i) , std_X_for(j,i), a_X_for(j,i) ] = moments_sample_data(X_t_analysis_for{i}(j,:));
            % empirical distribution function T_fs
                    F_N_Tfs_for{i}(j,:) = F_N_max(X_t_analysis_for{i}(j,:)); 
            % empirical distribution function T3600
                    F_N_T3600_for{i}(j,:) = F_N_Tfs_for{i}(j,:).^(T_3600/i);
            % model parameters
                    [alpha1_MLE_for(j,i) , u1_MLE_for(j,i), nlogL_G_MLE_for(1,i)] = MLE_Gumbel_max(X_t_analysis_for{i}(j,:)); % Gumbel-MLE
                    parGEVfor{j}(:,i) = gevfit(X_t_analysis_for{i}(j,:));
            % settings analytic distribution functions Tfs
                    xmin_max = min(X_t_analysis_for{i}(j,:)) ; 
                    xmax_max = max(X_t_analysis_for{i}(j,:))+z*abs(std_X_for(j,i));
                    x_for{i}(j,:) = linspace(xmin_max,xmax_max, n);
            % analytical distribution function Tfs
                     Fx_GEV_MLE_Tfs_for{i}(j,:) = gevcdf(x_for{i}(j,:),parGEVfor{j}(1,i),parGEVfor{j}(2,i),parGEVfor{j}(3,i));
                     Fx_G_MLE_Tfs_for{i}(j,:) = cdf_Gumbel_max(alpha1_MLE_for(j,i) , u1_MLE_for(j,i), x_for{i}(j,:)); % Gumbel-MLE
                             
            % analytical distribution function T3600
                    Fx_GEV_MLE_T3600_for{i}(j,:) = Fx_GEV_MLE_Tfs_for{i}(j,:).^(T_3600/i);
                    Fx_G_MLE_T3600_for{i}(j,:) = Fx_G_MLE_Tfs_for{i}(j,:).^(T_3600/i);

              % DESIGN FRACTILES tS EXTREMES   
                    x_fractiles_GEV_MLE_for_Tfs{i}(j,:) = fractile_values_empirical_Fx(x_for{i}(j,:),Fx_GEV_MLE_Tfs_for{i}(j,:), P_fractiles.^(plotFs/T_3600));
                    x_fractiles_G_MLE_for_Tfs{i}(j,:) = fractile_values_empirical_Fx(x_for{i}(j,:),Fx_G_MLE_Tfs_for{i}(j,:), P_fractiles.^(plotFs/T_3600));

              % DESIGN FRACTILES HOURLY EXTREMES 
                   x_fractiles_GEV_MLE_for{i}(j,:) = fractile_values_empirical_Fx(x_for{i}(j,:),Fx_GEV_MLE_T3600_for{i}(j,:), P_fractiles);
                   x_fractiles_G_MLE_for{i}(j,:) = fractile_values_empirical_Fx(x_for{i}(j,:),Fx_G_MLE_T3600_for{i}(j,:), P_fractiles);

              % save distribution functions

                    if save_Fx_G_MLE == 1
                    str= [ 'C_Refbuilding_for_Fx_G_MLE_' num2str(i) ' = transpose([x_for{i}(1,:); Fx_G_MLE_T3600_for{i}(1,:)]);']; 
                    eval(str);  
                    str= [ 'dlmwrite(' char(39) 'fitted-dists\C_Refbuilding_for_Fx_G_MLE_' num2str(i) '.xpt' char(39) ',C_Refbuilding_for_Fx_G_MLE_' num2str(i) ',' char(39) 'delimiter' char(39) ',' char(39) '\t' char(39) ',' char(39) 'precision' char(39) ',16)'  ];  
                    eval(str)            
                    end    
                    if save_Fx_GEV_MLE == 1
                    str= [ 'C_Refbuilding_for_Fx_GEV_MLE_' num2str(i) ' = transpose([x_for{i}(1,:); Fx_GEV_MLE_T3600_for{i}(1,:)]);']; 
                    eval(str);  
                    str= [ 'dlmwrite(' char(39) 'fitted-dists\C_Refbuilding_for_Fx_GEV_MLE_' num2str(i) '.xpt' char(39) ',C_Refbuilding_for_Fx_GEV_MLE_' num2str(i) ',' char(39) 'delimiter' char(39) ',' char(39) '\t' char(39) ',' char(39) 'precision' char(39) ',16)'  ];  
                    eval(str)         
                    end   
        end
        end
        
elseif extremes == 3      
    for i =  T_fs          
        F_N_Tfs_stat{i} = [];
        x_mom{i} = [];
        for j = 1 : size(X_t_analysis_stat{i},1)  % for all bootstrapped samples, might also be just 1            
            % current sample data
                    X_t_analysis_stat{i}(j,:)=sort(X_t_analysis_stat{i}(j,:));  %  X = current (bootstrapped) sample
            % moments
                    [mean_X_stat(j,i) , std_X_stat(j,i), a_X_stat(j,i) ] = moments_sample_data(X_t_analysis_stat{i}(j,:));
            % empirical distribution function T_fs
                    F_N_Tfs_stat{i}(j,:) = F_N_max(X_t_analysis_stat{i}(j,:)); 
            % empirical distribution function T3600
                    F_N_T3600_stat{i}(j,:) = F_N_Tfs_stat{i}(j,:).^(T_3600/i);
            % model parameters
                    [alpha1_MLE_stat(j,i) , u1_MLE_stat(j,i), nlogL_G_MLE_stat(1,i)] = MLE_Gumbel_max(X_t_analysis_stat{i}(j,:)); % Gumbel-MLE
                    parGEVstat{j}(:,i) = gevfit(X_t_analysis_stat{i}(j,:),[],'MaxIter',500);
            % settings analytic distribution functions Tfs
                    xmin_max = min(X_t_analysis_stat{i}(j,:)); 
                    xmax_max = max(X_t_analysis_stat{i}(j,:))+z*abs(std_X_stat(j,i));
                    x_stat{i}(j,:) = linspace(xmin_max,xmax_max, n);
            % analytical distribution function Tfs
                     Fx_GEV_MLE_Tfs_stat{i}(j,:) = gevcdf(x_stat{i}(j,:),parGEVstat{j}(1,i),parGEVstat{j}(2,i),parGEVstat{j}(3,i));
                     Fx_G_MLE_Tfs_stat{i}(j,:) = cdf_Gumbel_max(alpha1_MLE_stat(j,i) , u1_MLE_stat(j,i), x_stat{i}(j,:)); % Gumbel-MLE
                             
            % analytical distribution function T3600
                    Fx_GEV_MLE_T3600_stat{i}(j,:) = Fx_GEV_MLE_Tfs_stat{i}(j,:).^(T_3600/i);
                    Fx_G_MLE_T3600_stat{i}(j,:) = Fx_G_MLE_Tfs_stat{i}(j,:).^(T_3600/i);

              % DESIGN FRACTILES tS EXTREMES   
                    x_fractiles_GEV_MLE_stat_Tfs{i}(j,:) = fractile_values_empirical_Fx(x_stat{i}(j,:),Fx_GEV_MLE_Tfs_stat{i}(j,:), P_fractiles.^(plotFs/T_3600));
                    x_fractiles_G_MLE_stat_Tfs{i}(j,:) = fractile_values_empirical_Fx(x_stat{i}(j,:),Fx_G_MLE_Tfs_stat{i}(j,:), P_fractiles.^(plotFs/T_3600));

              % DESIGN FRACTILES HOURLY EXTREMES 
                   x_fractiles_GEV_MLE_stat{i}(j,:) = fractile_values_empirical_Fx(x_stat{i}(j,:),Fx_GEV_MLE_T3600_stat{i}(j,:), P_fractiles);
                   x_fractiles_G_MLE_stat{i}(j,:) = fractile_values_empirical_Fx(x_stat{i}(j,:),Fx_G_MLE_T3600_stat{i}(j,:), P_fractiles);

        end
    end
end


%% Plot fitted distribution functions for sample data for plotFs

if plot_figures_model == 1
close all
    
    if extremes == 1
    k=1;
    h=figure;
                        F_N_Tfs_mom_transformed = [];
                        F_N_T3600_mom_transformed = [];
                        Fx_Tfs_mom = [];
                        Fx_T3600_mom = [];
for i =  plotFs
    if plot_tT_extremes == 1
        P_fractiles = P_fractiles.^(plotFs/T_3600);
    else
    end
                    for j = 1:b; % plot original data                
                  % plot settings
                        xmin{i} = min(x_mom{i}(1,:));
                        xmax{i} = max(X_t_analysis_mom{i}(1,:))+0.5;
                        
                        if plot_tT_extremes == 2
                            ymin = 0 ;
                            ymax = 0.99;
                        elseif plot_tT_extremes == 3
                            ymin = 0 ;
                            ymax = 0.99999; 
                        else
                            ymin = min(F_N_T3600_mom{i}(j,:)) ;
                            ymax = 0.99999;
                        end

                   % plot distribution functions
                        % summary of analyical distribution functions

                        Fx_Tfs_mom{i}{j} = [Fx_G_MLE_Tfs_mom{i}(j,:); Fx_GEV_MLE_Tfs_mom{i}(j,:)];
                        Fx_T3600_mom{i}{j} = [Fx_G_MLE_T3600_mom{i}(j,:); Fx_GEV_MLE_T3600_mom{i}(j,:)];

                    % fractile values
                        x_fractiles_F_N_mom{i} = [];
                        x_fractiles_F_N_mom{i} = fractile_values_empirical_Fx(X_t_analysis_mom{i}(j,:), F_N_T3600_mom{i}(j,:), P_cook);
                        x_cook_F_N_mom(:, i) = x_fractiles_F_N_mom{i}(1);
    

                     % transform distribution functions for different plot settings       
                          if plot_set == 1
                                    F_N_Tfs_mom_transformed{i} = F_N_Tfs_mom{i}(j,:);
                                    F_N_T3600_mom_transformed{i} = F_N_T3600_mom{i}(j,:);
                                    ymin_transformed{i} = ymin;
                                    ymax_transformed{i} = ymax;
                                    Fx_Tfs_mom_transformed{i}{j} =  Fx_Tfs_mom{i}{j}; 
                                    Fx_T3600_mom_transformed{i}{j} =  Fx_T3600_mom{i}{j}; 
                                    P_fractiles_transformed = P_fractiles;                                          
                                    ytick_transf = ytick;

                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}_{me}\leq C)_{t}$';
                                    else 
                                    str = '$P(\hat{c}_{me} \leq C)_{1hr}$';
                                    end
                          elseif plot_set == 2
                                    F_N_Tfs_mom_transformed{i}(j,:) = log(1-F_N_Tfs_mom{i}(j,:));
                                    F_N_T3600_mom_transformed{i}(j,:) = log(1-F_N_T3600_mom{i}(j,:));
                                    ymin_transformed{i} = reallog(1-ymax);
                                    ymax_transformed{i} = reallog(1-ymin);
                                    Fx_Tfs_mom_transformed{i}{j} =  log(1-Fx_Tfs_mom{i}{j}); 
                                    Fx_T3600_mom_transformed{i}{j} =  log(1-Fx_T3600_mom{i}{j}); 
                                    P_fractiles_transformed = log(1-P_fractiles);                                          
                                    ytick_transf = log(1-ytick);

                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}_{me} \req C)_{t}$';
                                    else 
                                    str = '$P(\hat{c}_{me} \req C)_{1hr}$';
                                    end
                          elseif plot_set == 3
                                    F_N_Tfs_mom_transformed{i}(j,:) = -log(-log(F_N_Tfs_mom{i}(j,:)));
                                    F_N_T3600_mom_transformed{i}(j,:) = -log(-log(F_N_T3600_mom{i}(j,:)));
                                    Fx_Tfs_mom_transformed{i}{j} =  -log(-log(Fx_Tfs_mom{i}{j})); 
                                    Fx_T3600_mom_transformed{i}{j} =  -log(-log(Fx_T3600_mom{i}{j})); 
                                    P_fractiles_transformed = -log(-log(P_fractiles));

                                    ymin_transformed{i} = -log(-log(0.000000000000001));
                                    ymax_transformed{i} = -log(-log(ymax));                        
                                    ytick_transf = -log(-log(ytick));
                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}\leq \hat{c}_M)_{t}$';
                                    else 
                                    str = '$P(\hat{c} \leq \hat{c}_M)_{1hr}$';
                                    end
                          end 

                                Fx_G_MLE_Tfs_mom_transformed{i}{j} = Fx_Tfs_mom_transformed{i}{j}(1,:);
                                Fx_GEV_MLE_Tfs_mom_transformed{i}{j} = Fx_Tfs_mom_transformed{i}{j}(2,:);

                                Fx_G_MLE_T3600_mom_transformed{i}{j} = Fx_T3600_mom_transformed{i}{j}(1,:);
                                Fx_GEV_MLE_T3600_mom_transformed{i}{j} = Fx_T3600_mom_transformed{i}{j}(2,:);

                    % plot 

                        if plot_tT_extremes == 1
               
                           if i == plotFs % Only show bootstrapping for the wanted t-extreme distribution functions
                               if plot_scatter == 1
                                   if j == 1 
                               scatter(X_t_analysis_mom{i}(1,:), F_N_Tfs_mom_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               %scatter(X_t_analysis_mom{i}(1,:), F_N_Tfs_mom_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'DisplayName','sample data')
                                   else
                                   end
                               else
                               end
                           hold all
                               if j == 1
                                   h1 = plot(x_mom{i}(j,:), Fx_G_MLE_Tfs_mom_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1.5,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_mom(j,i),2) ',' num2str(u1_MLE_mom(j,i),2) '), t = ' num2str(plotFs,2) 's' ]);
                                   h2 = plot(x_mom{i}(j,:), Fx_GEV_MLE_Tfs_mom_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1.5,  'DisplayName', ['Type III (' num2str(parGEVmom{j}(1,i),2) ',' num2str(parGEVmom{j}(2,i),2) ',' num2str(parGEVmom{j}(3,i),2) '), t = ' num2str(plotFs,2) 's' ]);
                                   %h1 = plot(x_mom{i}(j,:), Fx_G_MLE_Tfs_mom_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1.5,'DisplayName', ['Type I, t = ' num2str(plotFs,2) 's' ]);
                                   %h2 = plot(x_mom{i}(j,:), Fx_GEV_MLE_Tfs_mom_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1.5,  'DisplayName', ['Type III, t = ' num2str(plotFs,2) 's' ]);

                               else
                                    plot(x_mom{i}(j,:), Fx_GEV_MLE_Tfs_mom_transformed{i}{j}, '-', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_mom{i}(j,:), Fx_G_MLE_Tfs_mom_transformed{i}{j}, '-', 'Color', colors2{k}, 'LineWidth', 1, 'HandleVisibility','off');
                               end
                                %[ horizontal_lines, vertical_lines ] = plot_fractile_values( xmin{i}, xmax{i}, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_GEV_MLE_mom_Tfs{i}(1,:),'DisplayName','Cook-Mayne design value');
                                [ horizontal_lines, vertical_lines ] = plot_fractile_values( xmin{i}, xmax{i}, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_GEV_MLE_mom_Tfs{i}(1,:),'HandleVisibility','off');
                                [ horizontal_lines, vertical_lines ] = plot_fractile_values( xmin{i}, xmax{i}, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_G_MLE_mom_Tfs{i}(1,:), 'HandleVisibility','off' );
                                %line([design_C_G design_C_G],[ymin_transformed{i} P_fractiles_transformed], 'Color', [0 0 0]+0.05*10,'DisplayName','Computed design value');
                                %line([design_C_GEV_GEV design_C_GEV_GEV],[ymin_transformed{i} P_fractiles_transformed], 'Color', [0 0 0]+0.05*10, 'HandleVisibility','off');
                                axis([xmin{i} xmax{i} ymin_transformed{i} ymax_transformed{i}])
                           else
          
                           end
                           
                        elseif plot_tT_extremes == 2  
                            % Shifted to hourly extremes
                           if i == plotFs % Only show bootstrapping for the wanted t-extreme distribution functions
                               if plot_scatter == 1
                               scatter(X_t_analysis_mom{i}(1,:), F_N_T3600_mom_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               hold all
                               else
                               end
                           
                               if j == 1
                                   h1 = plot(x_mom{i}(j,:), Fx_G_MLE_T3600_mom_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_mom(j,i),2) ',' num2str(u1_MLE_mom(j,i),2) '), shifted ' num2str(plotFs,2) 's' ]);
                                   h2 = plot(x_mom{i}(j,:), Fx_GEV_MLE_T3600_mom_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVmom{j}(1,i),2) ',' num2str(parGEVmom{j}(2,i),2) ',' num2str(parGEVmom{j}(3,i),2) '), shifted ' num2str(plotFs,2) 's' ]);
                                    
                               else
                                    plot(x_mom{i}(j,:), Fx_GEV_MLE_T3600_mom_transformed{i}{j}, '--', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_mom{i}(j,:), Fx_G_MLE_T3600_mom_transformed{i}{j}, '--', 'Color', colors2{k+1}, 'LineWidth', 1, 'HandleVisibility','off');

                               end
                               axis([min(X_t_analysis_mom{i}(1,:)+0.3) xmax{i} ymin_transformed{i} ymax_transformed{i}])
                           else
          
                           end                          
                        
                        elseif plot_tT_extremes == 3  
                            % Original + shifted to hourly extremes
                           if i == plotFs % Only show bootstrapping for the wanted t-extreme distribution functions
                               if plot_scatter == 1
                               scatter(X_t_analysis_mom{i}(1,:), F_N_T3600_mom_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               hold all
                               scatter(X_t_analysis_mom{i}(1,:), F_N_Tfs_mom_transformed{i}(1,:), 40,'o', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               else
                               end
                           
                               if j == 1
                                    plot(x_mom{i}(j,:), Fx_GEV_MLE_T3600_mom_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVmom{j}(1,i),2) ',' num2str(parGEVmom{j}(2,i),2) ',' num2str(parGEVmom{j}(3,i),2) '), shifted ' num2str(T_fs(k),2) 's' ]);
                                    plot(x_mom{i}(j,:), Fx_G_MLE_T3600_mom_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_mom(j,i),2) ',' num2str(u1_MLE_mom(j,i),2) '), shifted ' num2str(T_fs(k),2) 's' ]);
                                    plot(x_mom{i}(j,:), Fx_GEV_MLE_Tfs_mom_transformed{i}{j}, '--', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVmom{j}(1,i),2) ',' num2str(parGEVmom{j}(2,i),2) ',' num2str(parGEVmom{j}(3,i),2) '),' num2str(T_fs(k),2) 's' ]);
                                    plot(x_mom{i}(j,:), Fx_G_MLE_Tfs_mom_transformed{i}{j}, '--', 'Color', colors{k+1}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_mom(j,i),2) ',' num2str(u1_MLE_mom(j,i),2) '),' num2str(T_fs(k),2) 's' ]);

                               else
                                    plot(x_mom{i}(j,:), Fx_GEV_MLE_T3600_mom_transformed{i}{j}, '--', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_mom{i}(j,:), Fx_G_MLE_T3600_mom_transformed{i}{j}, '--', 'Color', colors2{k+1}, 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_mom{i}(j,:), Fx_GEV_MLE_Tfs_mom_transformed{i}{j}, '--', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_mom{i}(j,:), Fx_G_MLE_Tfs_mom_transformed{i}{j}, '--', 'Color', colors2{k+1}, 'LineWidth', 1, 'HandleVisibility','off');

                               end
                               axis([min(X_t_analysis_mom{i}+0.3) xmax{i} ymin_transformed{i} ymax_transformed{i}])
                           else
          
                           end
                            

                        end               

                    end
                    xtick = round(xmin{i},1):0.2:round(xmax{i},1);
                    k=k+1;
                    set(gca,'YTick',ytick_transf);
                    set(gca,'XTick',xtick);
                    set(gca,'YTickLabel',ytick);
                    set(gca,'fontsize',fontsize);
                    xlabel('Moment coefficient $\hat{c}_M$ [-]', 'Interpreter', 'Latex')
                    ylabel(str,'Interpreter', 'Latex')                      
                    grid minor
                    legend('show', 'Location', 'northwest')
                    uistack(h1, 'top')
                    uistack(h2, 'top')
                    if save_figures_model == 1;
                    set(gcf, 'PaperPositionMode','auto')
                    str=['saveSameSize(gcf,' char(39) 'format' char(39) ',' char(39) 'png' char(39) ',' char(39) 'file' char(39) ',' char(39) 'Reference_building_maxima_sec' num2str(i) '.png' char(39) ')'];
                    eval(str)
                    str=['movefile(' char(39) 'Reference_building_maxima_sec' num2str(i) '.png' char(39) ',' char(39) num2str(saving_directory) char(39) ')'] ;
                    eval(str)
                    end
end



elseif extremes == 2
    k=1;
    h=figure;
                        F_N_Tfs_for_transformed = [];
                        F_N_T3600_for_transformed = [];
                        Fx_Tfs_for = [];
                        Fx_T3600_for = [];
    for i =  plotFs
        if plot_tT_extremes == 1
            P_fractiles = P_fractiles.^(plotFs/T_3600);
        else
        end
                    for j = 1:b; % plot original data                
                  % plot settings
                        xmin{i} = min(x_for{i}(1,:));
                        xmax{i} = max(X_t_analysis_for{i}(1,:))+0.5;
                        
                        if plot_tT_extremes == 2
                            ymin = 0 ;
                            ymax = 0.99;
                        elseif plot_tT_extremes == 3
                            ymin = 0 ;
                            ymax = 0.99999; 
                        else
                            ymin = min(F_N_T3600_for{i}(j,:)) ;
                            ymax = 0.99999;
                        end

                   % plot distribution functions
                        % summary of analyical distribution functions

                        Fx_Tfs_for{i}{j} = [Fx_G_MLE_Tfs_for{i}(j,:); Fx_GEV_MLE_Tfs_for{i}(j,:)];
                        Fx_T3600_for{i}{j} = [Fx_G_MLE_T3600_for{i}(j,:); Fx_GEV_MLE_T3600_for{i}(j,:)];

                    % fractile values
                        x_fractiles_F_N_for{i} = [];
                        x_fractiles_F_N_for{i} = fractile_values_empirical_Fx(X_t_analysis_for{i}(j,:), F_N_T3600_for{i}(j,:), P_cook);
                        x_cook_F_N_for(:, i) = x_fractiles_F_N_for{i}(1);
    


                     % transform distribution functions for different plot settings       
                          if plot_set == 1
                                    F_N_Tfs_for_transformed{i} = F_N_Tfs_for{i}(j,:);
                                    F_N_T3600_for_transformed{i} = F_N_T3600_for{i}(j,:);
                                    ymin_transformed{i} = ymin;
                                    ymax_transformed{i} = ymax;
                                    Fx_Tfs_for_transformed{i}{j} =  Fx_Tfs_for{i}{j}; 
                                    Fx_T3600_for_transformed{i}{j} =  Fx_T3600_for{i}{j}; 
                                    P_fractiles_transformed = P_fractiles;                                          
                                    ytick_transf = ytick;

                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}_{me}\leq C)_{t}$';
                                    else 
                                    str = '$P(\hat{c}_{me} \leq C)_{1hr}$';
                                    end
                          elseif plot_set == 2
                                    F_N_Tfs_for_transformed{i}(j,:) = log(1-F_N_Tfs_for{i}(j,:));
                                    F_N_T3600_for_transformed{i}(j,:) = log(1-F_N_T3600_for{i}(j,:));
                                    ymin_transformed{i} = reallog(1-ymax);
                                    ymax_transformed{i} = reallog(1-ymin);
                                    Fx_Tfs_for_transformed{i}{j} =  log(1-Fx_Tfs_for{i}{j}); 
                                    Fx_T3600_for_transformed{i}{j} =  log(1-Fx_T3600_for{i}{j}); 
                                    P_fractiles_transformed = log(1-P_fractiles);                                          
                                    ytick_transf = log(1-ytick);

                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}_{me} \req C)_{t}$';
                                    else 
                                    str = '$P(\hat{c}_{me} \req C)_{1hr}$';
                                    end
                          elseif plot_set == 3
                                    F_N_Tfs_for_transformed{i}(j,:) = -log(-log(F_N_Tfs_for{i}(j,:)));
                                    F_N_T3600_for_transformed{i}(j,:) = -log(-log(F_N_T3600_for{i}(j,:)));
                                    Fx_Tfs_for_transformed{i}{j} =  -log(-log(Fx_Tfs_for{i}{j})); 
                                    Fx_T3600_for_transformed{i}{j} =  -log(-log(Fx_T3600_for{i}{j})); 
                                    P_fractiles_transformed = -log(-log(P_fractiles));

                                    ymin_transformed{i} = -log(-log(0.000000000000001));
                                    ymax_transformed{i} = -log(-log(ymax));                        
                                    ytick_transf = -log(-log(ytick));
                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}\leq C)_{t}$';
                                    else 
                                    str = '$P(\hat{c} \leq C)_{1hr}$';
                                    end
                          end 

                                Fx_G_MLE_Tfs_for_transformed{i}{j} = Fx_Tfs_for_transformed{i}{j}(1,:);
                                Fx_GEV_MLE_Tfs_for_transformed{i}{j} = Fx_Tfs_for_transformed{i}{j}(2,:);

                                Fx_G_MLE_T3600_for_transformed{i}{j} = Fx_T3600_for_transformed{i}{j}(1,:);
                                Fx_GEV_MLE_T3600_for_transformed{i}{j} = Fx_T3600_for_transformed{i}{j}(2,:);

                    % plot 

                        if plot_tT_extremes == 1
               
                           if i == plotFs % Only show bootstrapping for the wanted t-extreme distribution functions
                               if plot_scatter == 1
                               scatter(X_t_analysis_for{i}(1,:), F_N_Tfs_for_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')

                               else
                               end
                               hold all
                               
                               if j == 1
                                  h1 =  plot(x_for{i}(j,:), Fx_G_MLE_Tfs_for_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1.5,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_for(j,i),2) ',' num2str(u1_MLE_for(j,i),2) '), t = ' num2str(plotFs,2) 's' ]);
                                  h2 =  plot(x_for{i}(j,:), Fx_GEV_MLE_Tfs_for_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1.5,  'DisplayName', ['Type III (' num2str(parGEVfor{j}(1,i),2) ',' num2str(parGEVfor{j}(2,i),2) ',' num2str(parGEVfor{j}(3,i),2) '), t = ' num2str(plotFs,2) 's' ]);
                                    
                               else
                                    plot(x_for{i}(j,:), Fx_GEV_MLE_Tfs_for_transformed{i}{j}, '-', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_for{i}(j,:), Fx_G_MLE_Tfs_for_transformed{i}{j}, '-', 'Color', colors2{k}, 'LineWidth', 1, 'HandleVisibility','off');
                               end

                                [ horizontal_lines, vertical_lines ] = plot_fractile_values( xmin{i}, xmax{i}, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_GEV_MLE_for_Tfs{i}(1,:),'DisplayName','Cook-Mayne design value');
                                [ horizontal_lines, vertical_lines ] = plot_fractile_values( xmin{i}, xmax{i}, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_G_MLE_for_Tfs{i}(1,:), 'HandleVisibility','off' );
                                %line([design_C_G design_C_G],[ymin_transformed{i} P_fractiles_transformed], 'Color', [0 0 0]+0.05*10,'DisplayName','Computed design value');
                                %line([design_C_GEV_GEV design_C_GEV_GEV],[ymin_transformed{i} P_fractiles_transformed], 'Color', [0 0 0]+0.05*10, 'HandleVisibility','off');
                               axis([xmin{i} xmax{i} ymin_transformed{i} ymax_transformed{i}])
                           else
          
                           end
                           
                        elseif plot_tT_extremes == 2  
                            % Shifted to hourly extremes
                           if i == plotFs % Only show bootstrapping for the wanted t-extreme distribution functions
                               if plot_scatter == 1
                               scatter(X_t_analysis_for{i}(1,:), F_N_T3600_for_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               hold all
                               else
                               end
                           
                               if j == 1
                                   h1 = plot(x_for{i}(j,:), Fx_G_MLE_T3600_for_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_for(j,i),2) ',' num2str(u1_MLE_for(j,i),2) '), shifted ' num2str(T_fs(k),2) 's' ]);
                                   h2 = plot(x_for{i}(j,:), Fx_GEV_MLE_T3600_for_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVfor{j}(1,i),2) ',' num2str(parGEVfor{j}(2,i),2) ',' num2str(parGEVfor{j}(3,i),2) '), shifted ' num2str(T_fs(k),2) 's' ]);
                                    
                               else
                                    plot(x_for{i}(j,:), Fx_GEV_MLE_T3600_for_transformed{i}{j}, '--', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_for{i}(j,:), Fx_G_MLE_T3600_for_transformed{i}{j}, '--', 'Color', colors2{k+1}, 'LineWidth', 1, 'HandleVisibility','off');

                               end
                               axis([min(X_t_analysis_for{i}(1,:)+0.3) xmax{i} ymin_transformed{i} ymax_transformed{i}])
                           else
          
                           end                          
                        
                        elseif plot_tT_extremes == 3  
                            % Original + shifted to hourly extremes
                           if i == plotFs % Only show bootstrapping for the wanted t-extreme distribution functions
                               if plot_scatter == 1
                               scatter(X_t_analysis_for{i}(1,:), F_N_T3600_for_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               hold all
                               scatter(X_t_analysis_for{i}(1,:), F_N_Tfs_for_transformed{i}(1,:), 40,'o', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                               else
                               end
                           
                               if j == 1
                                    plot(x_for{i}(j,:), Fx_GEV_MLE_T3600_for_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVfor{j}(1,i),2) ',' num2str(parGEVfor{j}(2,i),2) ',' num2str(parGEVfor{j}(3,i),2) '), shifted ' num2str(T_fs(k),2) 's' ]);
                                    plot(x_for{i}(j,:), Fx_G_MLE_T3600_for_transformed{i}{j}, '-', 'Color', colors{k+1}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_for(j,i),2) ',' num2str(u1_MLE_for(j,i),2) '), shifted ' num2str(T_fs(k),2) 's' ]);
                                    plot(x_for{i}(j,:), Fx_GEV_MLE_Tfs_for_transformed{i}{j}, '--', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVfor{j}(1,i),2) ',' num2str(parGEVfor{j}(2,i),2) ',' num2str(parGEVfor{j}(3,i),2) '),' num2str(T_fs(k),2) 's' ]);
                                    plot(x_for{i}(j,:), Fx_G_MLE_Tfs_for_transformed{i}{j}, '--', 'Color', colors{k+1}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_for(j,i),2) ',' num2str(u1_MLE_for(j,i),2) '),' num2str(T_fs(k),2) 's' ]);

                               else
                                    plot(x_for{i}(j,:), Fx_GEV_MLE_T3600_for_transformed{i}{j}, '--', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_for{i}(j,:), Fx_G_MLE_T3600_for_transformed{i}{j}, '--', 'Color', colors2{k+1}, 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_for{i}(j,:), Fx_GEV_MLE_Tfs_for_transformed{i}{j}, '--', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                                    plot(x_for{i}(j,:), Fx_G_MLE_Tfs_for_transformed{i}{j}, '--', 'Color', colors2{k+1}, 'LineWidth', 1, 'HandleVisibility','off');

                               end
                               axis([min(X_t_analysis_for{i}+0.3) xmax{i} ymin_transformed{i} ymax_transformed{i}])
                           else
          
                           end
                            

                        end               

                    end
                    xtick = round(xmin{i},1):0.2:round(xmax{i},1);
                    k=k+1;
                    set(gca,'YTick',ytick_transf);
                    set(gca,'XTick',xtick);
                    set(gca,'YTickLabel',ytick);
                    xlabel('Dynamic force coefficient C [-]', 'Interpreter', 'Latex')
                    ylabel(str,'Interpreter', 'Latex')                      
                    grid minor
                    legend('show', 'Location', 'northwest')
                    uistack(h1, 'top')
                    uistack(h2, 'top')
                    if save_figures_model ==1;
                    set(gcf, 'PaperPositionMode','auto')
                    str=['saveSameSize(gcf,' char(39) 'format' char(39) ',' char(39) 'png' char(39) ',' char(39) 'file' char(39) ',' char(39) 'Reference_building_maxima_sec' num2str(i) '.png' char(39) ')'];
                    eval(str)
                    str=['movefile(' char(39) 'Reference_building_maxima_sec' num2str(i) '.png' char(39) ',' char(39) num2str(saving_directory) char(39) ')'] ;
                    eval(str)
                    end
    end 
    end
else 
end
%% Plot fitted distribution functions for sample data for T_fs
if plotblocks == 1 
     k=1;
        for i =  T_fs
                    for j = 1:b; % plot original data                
                  % plot settings
                        xmin{i} = min(x_stat{i});
                        xmax{i} = max(X_t_analysis_stat{i})+0.3;
                        
                        if plot_tT_extremes == 2
                            ymin = 0 ;
                            ymax = 0.99;
                        elseif plot_tT_extremes == 3
                            ymin = 0 ;
                            ymax = 0.99999; 
                        else
                            ymin = min(F_N_T3600_stat{i}(j,:)) ;
                            ymax = 0.99999;
                        end

                   % plot distribution functions
                        % summary of analyical distribution functions

                        Fx_Tfs_stat{i}{j} = [Fx_G_MLE_Tfs_stat{i}(j,:); Fx_GEV_MLE_Tfs_stat{i}(j,:)];
                        Fx_T3600_stat{i}{j} = [Fx_G_MLE_T3600_stat{i}(j,:); Fx_GEV_MLE_T3600_stat{i}(j,:)];

                    % fractile values
                        x_fractiles_F_N_stat{i} = [];
                        x_fractiles_F_N_stat{i} = fractile_values_empirical_Fx(X_t_analysis_stat{i}(j,:), F_N_T3600_stat{i}(j,:), P_fractiles);
                        x_cook_F_N_stat(:, i) = x_fractiles_F_N_stat{i}(1);

                     % transform distribution functions for different plot settings       
                          if plot_set == 1
                                    F_N_Tfs_stat_transformed{i} = F_N_Tfs_stat{i}(j,:);
                                    F_N_T3600_stat_transformed{i} = F_N_T3600_stat{i}(j,:);
                                    ymin_transformed{i} = ymin;
                                    ymax_transformed{i} = ymax;
                                    Fx_Tfs_stat_transformed{i}{j} =  Fx_Tfs_stat{i}{j}; 
                                    Fx_T3600_stat_transformed{i}{j} =  Fx_T3600_stat{i}{j}; 
                                    P_fractiles_transformed = P_fractiles;                                          
                                    ytick_transf =ytick;

                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}_{me}\leq C)_{t}$';
                                    else 
                                    str = '$P(\hat{c}_{me} \leq C)_{1hr}$';
                                    end
                          elseif plot_set == 2
                                    F_N_Tfs_stat_transformed{i}(j,:) = log(1-F_N_Tfs_stat{i}(j,:));
                                    F_N_T3600_stat_transformed{i}(j,:) = log(1-F_N_T3600_stat{i}(j,:));
                                    ymin_transformed{i} = reallog(1-ymax);
                                    ymax_transformed{i} = reallog(1-ymin);
                                    Fx_Tfs_stat_transformed{i}{j} =  log(1-Fx_Tfs_stat{i}{j}); 
                                    Fx_T3600_stat_transformed{i}{j} =  log(1-Fx_T3600_stat{i}{j}); 
                                    P_fractiles_transformed = log(1-P_fractiles);                                          
                                    ytick_transf = log(1-ytick);

                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}_{me} \req C)_{t}$';
                                    else 
                                    str = '$P(\hat{c}_{me} \req C)_{1hr}$';
                                    end
                          elseif plot_set == 3
                                    F_N_Tfs_stat_transformed{i}(j,:) = -log(-log(F_N_Tfs_stat{i}(j,:)));
                                    F_N_T3600_stat_transformed{i}(j,:) = -log(-log(F_N_T3600_stat{i}(j,:)));
                                    Fx_Tfs_stat_transformed{i}{j} =  -log(-log(Fx_Tfs_stat{i}{j})); 
                                    Fx_T3600_stat_transformed{i}{j} =  -log(-log(Fx_T3600_stat{i}{j})); 
                                    P_fractiles_transformed = -log(-log(P_fractiles));

                                    ymin_transformed{i} = -log(-log(0.000000000000001));
                                    ymax_transformed{i} = -log(-log(ymax));                        
                                    ytick_transf = -log(-log(ytick));
                                    if plot_tT_extremes == 1
                                    str = '$P(\hat{c}\leq C)_{t}$';
                                    else 
                                    str = '$P(\hat{c} \leq C)_{1hr}$';
                                    end
                          end 

                                Fx_G_MLE_Tfs_stat_transformed{i}{j} = Fx_Tfs_stat_transformed{i}{j}(1,:);
                                Fx_GEV_MLE_Tfs_stat_transformed{i}{j} = Fx_Tfs_stat_transformed{i}{j}(2,:);

                                Fx_G_MLE_T3600_stat_transformed{i}{j} = Fx_T3600_stat_transformed{i}{j}(1,:);
                                Fx_GEV_MLE_T3600_stat_transformed{i}{j} = Fx_T3600_stat_transformed{i}{j}(2,:);
                    end
        end
        for i = T_fs
           for j = 1:b; % plot original data          
               if plot_scatter == 1
                    if i == T_3600
               scatter(X_t_analysis_stat{i}(1,:), F_N_T3600_stat_transformed{i}(1,:), 40,'+', 'MarkerEdgeColor','black', 'HandleVisibility','off')
                    else
                    end
               else
               end
               hold all
               if i == T_3600 % Only show bootstrapping for the hourly extreme distribution functions
                   if j == 1
                        plot(x_stat{i}(j,:), Fx_G_MLE_T3600_stat_transformed{i}{j}, '--', 'Color', colors{k}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_stat(j,i),2) ',' num2str(u1_MLE_stat(j,i),2) '), hourly extremes' ]);
                        plot(x_stat{i}(j,:), Fx_GEV_MLE_T3600_stat_transformed{i}{j}, '-', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVstat{j}(1,i),2) ',' num2str(parGEVstat{j}(2,i),2) ',' num2str(parGEVstat{j}(3,i),2) '), hourly extremes' ]);
                   else
                        plot(x_stat{i}(j,:), Fx_GEV_MLE_T3600_stat_transformed{i}{j}, '-', 'Color', colors2{k} , 'LineWidth', 1, 'HandleVisibility','off');
                        plot(x_stat{i}(j,:), Fx_G_MLE_T3600_stat_transformed{i}{j}, '--', 'Color', colors2{k}, 'LineWidth', 1, 'HandleVisibility','off');
                   end
               else
                   if j == 1
                   plot(x_stat{i}(1,:), Fx_G_MLE_T3600_stat_transformed{i}{1}, '--', 'Color', colors{k}, 'LineWidth', 1,'DisplayName', ['Type I (' num2str(1/alpha1_MLE_stat(j,i),2) ',' num2str(u1_MLE_stat(j,i),2) '),' num2str(T_fs(k),2) 's' ]);
                   plot(x_stat{i}(1,:), Fx_GEV_MLE_T3600_stat_transformed{i}{1}, '-', 'Color', colors{k} , 'LineWidth', 1,  'DisplayName', ['Type III (' num2str(parGEVstat{j}(1,i),2) ',' num2str(parGEVstat{j}(2,i),2) ',' num2str(parGEVstat{j}(3,i),2) '),' num2str(T_fs(k),2) 's' ]);
                   
                   else

                   end
               end
               axis([xmin{T_3600} xmax{T_3600} ymin_transformed{T_3600} ymax_transformed{T_3600}])

               [ horizontal_lines, vertical_lines ] = plot_fractile_values( 0.3, 1.4, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_GEV_MLE_stat{i}(j,:) );
               [ horizontal_lines, vertical_lines ] = plot_fractile_values( 0.3, 1.4, ymin_transformed{i}, P_fractiles_transformed, x_fractiles_G_MLE_stat{i}(j,:) );
           end
 
               k=k+1;
        end
                    set(gca,'YTick',ytick_transf);
                    set(gca,'YTickLabel',ytick);
                    xlabel('Static moment coefficient C [-]', 'Interpreter', 'Latex')
                    ylabel(str,'Interpreter', 'Latex')                      
                    grid minor
                    legend('show', 'Location', 'northeast')
                    if save_plotblocks == 1;
                    set(gcf, 'PaperPositionMode','auto')
                    str=['saveSameSize(gcf,' char(39) 'format' char(39) ',' char(39) 'png' char(39) ',' char(39) 'file' char(39) ',' char(39) 'Plot_blocks' num2str(plot_tT_extremes) '' num2str(plot_set) '.png' char(39) ')'];
                    eval(str)
                    str=['movefile(' char(39) 'Plot_blocks' num2str(plot_tT_extremes) '' num2str(plot_set) '.png' char(39) ',' char(39) num2str(saving_directory) char(39) ')'] ;
                    eval(str)
                    end

else
end

%% Overview of Fx moments compared to other moments. 


if extremes == 1

% initialize
summary_mean_mom = [];
summary_std_mom = [];
summary_skewness_mom = [];
inputparameters_mom = [];
sample_statistics_mom = [];    

% run

[mean_G_MLE_mom, std_G_MLE_mom, skewness_G_MLE_mom] = moments_TypeI(alpha1_MLE_mom(1,plotFs), u1_MLE_mom(1,plotFs));
[mean_GEV_MLE_mom, std_GEV_MLE_mom, skewness_GEV_MLE_mom] = moments_TypeIII(parGEVmom{1}(:,plotFs));


summary_mean_mom = transpose([ mean_X_mom(1,plotFs); mean_G_MLE_mom; mean_GEV_MLE_mom ]);
summary_std_mom = transpose([std_X_mom(1,plotFs) ; std_G_MLE_mom; std_GEV_MLE_mom ]);
summary_skewness_mom = transpose([a_X_mom(1,plotFs) ; skewness_G_MLE_mom; skewness_GEV_MLE_mom ]);                
inputparameters_mom = transpose([   alpha1_MLE_mom(1,plotFs); u1_MLE_mom(1,plotFs); parGEVmom{1}(1,plotFs); parGEVmom{1}(2,plotFs); parGEVmom{1}(3,plotFs)]);
sample_statistics_mom = transpose( [  mean_X_mom(1,plotFs) ; std_X_mom(1,plotFs) ; std_X_mom(1,plotFs)./mean_X_mom(1,plotFs)  ; a_X_mom(1,plotFs) ]);  

% fractile values % cook fractile
j=1;


x_fractiles_G_MLE_mom_summary = x_fractiles_G_MLE_mom{plotFs}(j,1); 
x_fractiles_GEV_MLE_mom_summary = x_fractiles_GEV_MLE_mom{plotFs}(j,1);


x_fractiles_mom = transpose([x_fractiles_G_MLE_mom_summary; x_fractiles_GEV_MLE_mom_summary]);


elseif extremes == 2

% initialize
summary_mean_for = [];
summary_std_for = [];
summary_skewness_for = [];
inputparameters_for = [];
sample_statistics_for = [];    

% run

[mean_G_MLE_for, std_G_MLE_for, skewness_G_MLE_for] = moments_TypeI(alpha1_MLE_for(1,plotFs), u1_MLE_for(1,plotFs));
[mean_GEV_MLE_for, std_GEV_MLE_for, skewness_GEV_MLE_for] = moments_TypeIII(parGEVfor{1}(:,plotFs));


summary_mean_for = transpose([ mean_X_for(1,plotFs); mean_G_MLE_for; mean_GEV_MLE_for ]);
summary_std_for = transpose([std_X_for(1,plotFs) ; std_G_MLE_for; std_GEV_MLE_for ]);
summary_skewness_for = transpose([a_X_for(1,plotFs) ; skewness_G_MLE_for; skewness_GEV_MLE_for ]);                
inputparameters_for = transpose([   alpha1_MLE_for(1,plotFs); u1_MLE_for(1,plotFs); parGEVfor{1}(1,plotFs); parGEVfor{1}(2,plotFs); parGEVfor{1}(3,plotFs)]);
sample_statistics_for = transpose( [  mean_X_for(1,plotFs) ; std_X_for(1,plotFs) ; std_X_for(1,plotFs)./mean_X_for(1,plotFs)  ; a_X_for(1,plotFs) ]);  

% fractile values % cook fractile
j=1;


x_fractiles_G_MLE_for_summary = x_fractiles_G_MLE_for{plotFs}(j,1); 
x_fractiles_GEV_MLE_for_summary = x_fractiles_GEV_MLE_for{plotFs}(j,1);


x_fractiles_for = transpose([x_fractiles_G_MLE_for_summary; x_fractiles_GEV_MLE_for_summary]);

end

if moments_and_forces == 1
    summary_statistics_complete = transpose([summary_mean_mom summary_mean_for ; summary_std_mom summary_std_for ; summary_skewness_mom summary_skewness_for ]) ; 
end
    
%% EPDF coefficients
if plot_pdf == 1
bins = 50;
if extremes == 1
    figure
f = histogram(Cm_max{plotFs},bins,'Normalization','probability') ;
            
%             figure
%             bar(x,f/(sum(f)),'BarWidth', 1,'FaceColor',[1 .8 .8]);
elseif extremes == 2
    figure
f = histogram(Cf_max{plotFs},bins,'Normalization','probability') ;
            
%              figure
%              bar(x,f),'BarWidth', 1,'FaceColor',[1 .8 .8]); 
end
end
%% clear variables

clearvars -except x_fractiles_GEV_MLE_for_Tfs design_C_G design_C_GEV_GEV Cf_max Cm_max Cm_stat F_N_T3600_for F_N_T3600_mom F_N_Tfs_for F_N_Tfs_mom N N_stat T_fs T_3600 Fx_Tfs_mom Fx_T3600_mom Fx_Tfs_for Fx_T3600_for inputparameters_mom inputparameters_for sample_statistics_for sample_statistics_mom summary_mean_for summary_mean_mom summary_skewness_for summary_skewness_mom summary_std_for summary_std_mom x_fractiles_F_N_for x_fractiles_F_N_mom x_fractiles_mom x_fractiles_for summary_statistics_complete u1_MLE_mom mean_G_MLE_mom

