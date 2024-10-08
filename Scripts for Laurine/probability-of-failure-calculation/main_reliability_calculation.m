clear all
close all
clc

%% Settings and input data
%p_sec=[0.06 0.07 0.08 0.05 0.06 0.09 0.13 0.12 0.11 0.08 0.07 0.06];
% choices
        location = 'Schiphol';
        building = 'Reference building';
        configuration_V = 1;
        configuration_C = 1;
%         extremes = 3; 
        both_extremes = 1; 

        prob2b_GEV_GEV = 1;
%         prob2b_G_G = 1;
        prob2b_G_GEV = 1;
        run_all = 1;
% figure set up        
        fontsize = 12;
% input data
        run('initialize_structure_stochastic_variables');
        
% Prob2B settings
        addpath 'Prob2B'
        func = @myZfunction; % will be overwritten later
        %correlations = eye([6 6]);
        correlations =[];
        settings = struct('itmax',1000,'relaxf',0.25,'epsB',0.01,'epsZ',0.01);
        %settings = struct('itmax',100,'relaxf',1.0,'epsB',0.001,'epsZ',0.001,'epsH',0.001,'maxloop',10);

        
%% Run prob2B  G - G 
% if prob2b_G_G == 1
%     
%              
%                 % define parameters of distribution functions
%                 C_filename_G = [C_directory_G '\' C_filenames_G];  
%                 T = [];
%                 T_tmp = dlmread(C_filename_G,'\t');
%                 [p,ia,ic] = unique(T_tmp(:,2));   
%                 T(:,1) = T_tmp(ia,1);
%                 T(:,2) = p;
%                 parameters(7).data = T;
%                 
%                 % here it is determined which C_section corresponds to which limit state
%                 % function and which R-value
%     if extremes == 1
%                         func = @myZfunctionM2;
%                         % strength, blijft hetzelfde voor alle maxima
% 
%     elseif extremes == 2
%                         func = @myZfunctionQ2;
%                         % strength, blijft hetzelfde voor alle maxima
% 
%     end
%                     % roughness factor, blijft hetzelfde in elke limit state function
%                     parameters(8).gemX = crmean_squared;
%                     parameters(8).sigX = crstd_squared;
%                     parameters(8).x = parameters(4).gemX;
%          
%             
%                 V_filename_G = [V_directory_G '\' V_filenames_G];
%                 T = [];
%                 T_tmp = dlmread(V_filename_G,'\t');
%                 [p,ia,ic] = unique(T_tmp(:,2));   
%                 T(:,1) = T_tmp(ia,1);
%                 T(:,2) = p;
%                 parameters(6).data = T;
%                 parameters(10).gemX = 1;
%                 parameters(10).sigX = V_cov_G;
%                 
%         % the run        
%         [results, parametersOut] = Form(func, parameters, correlations, settings);
%         %% 
%         
%          beta_values_G = results.beta;
% 
%         Pf_values_G = results.Pf;
%         alfa_theta_G = parametersOut(1,1).alfa;
%         design_theta_G = parametersOut(1,1).x;
%         alfa_rho_G = parametersOut(2,1).alfa;
%         design_rho_G = parametersOut(2,1).x;
%         alfa_fy_G = parametersOut(3,1).alfa;
%         design_fy_G = parametersOut(3,1).x;
%         alfa_fc_G = parametersOut(4,1).alfa;
%         design_fc_G = parametersOut(4,1).x;
%         alfa_t_G = parametersOut(5,1).alfa;
%         design_t_G = parametersOut(5,1).x;
%         alfa_V_G = parametersOut(6,1).alfa;
%         design_V_G = parametersOut(6,1).x;
%         alfa_C_G = parametersOut(7,1).alfa;
%         design_C_G = parametersOut(7,1).x;
%         alfa_cr_G = parametersOut(8,1).alfa;
%         design_cr_G= parametersOut(8,1).x;
%         alfa_X_G = parametersOut(9,1).alfa;
%         design_X_G= parametersOut(9,1).x;
%         alfa_Sv_G = parametersOut(10,1).alfa;
%         design_Sv_G= parametersOut(10,1).x;
%         alfa_cd_G = parametersOut(11,1).alfa;
%         design_cd_G= parametersOut(11,1).x;
% 
%  
%     if extremes == 1
%     alfas_mean_G_mom = [ alfa_V_G; alfa_C_G; alfa_cr_G; alfa_X_G; alfa_Sv_G ; alfa_cd_G; alfa_theta_G; alfa_fy_G; alfa_fc_G; alfa_t_G];
%     beta_value_G_mom = beta_values_G;
%     design_values_G = [ design_theta_G; design_rho_G; design_fy_G; design_fc_G; design_t_G; design_V_G; design_C_G; design_cr_G; design_X_G; design_Sv_G ; design_cd_G];
% 
%     elseif extremes == 2
%     alfas_mean_G_for = [ alfa_V_G; alfa_C_G; alfa_cr_G; alfa_X_G; alfa_Sv_G ; alfa_cd_G; alfa_theta_G; alfa_fy_G; alfa_fc_G; alfa_t_G];
%     beta_value_G_for = beta_values_G;
% 
%     end
% end        
%% Run prob2B  G - GEV 

if prob2b_G_GEV == 1

    % here it is determined which C_section corresponds to which limit state
    % function and which R-value
    for extremes = 1:2
        run('read_windspeed_dists_Schiphol'); % needs input of Alex later
        run('read_building_specific_input');
        % define parameters of distribution functions
        C_filename_GEV = [C_directory_GEV '\' C_filenames_GEV];  
        T = [];
        T_tmp = dlmread(C_filename_GEV,'\t');
        [p,ia,ic] = unique(T_tmp(:,2));   
        T(:,1) = T_tmp(ia,1);
        T(:,2) = p;
        parameters(7).data = T;
        if extremes == 1
            func = @myZfunctionM2;
            % strength, stays the same for all maxima

        elseif extremes == 2
            func = @myZfunctionQ2;
            % strength, stays the same for all  maxima

        end
            % roughness factor, stays the same in every limit state function
            parameters(8).gemX = crmean_squared;
            parameters(8).sigX = crstd_squared;
            parameters(8).x = parameters(4).gemX;


            V_filename_G = [V_directory_G '\' V_filenames_G];
            T = [];
            T_tmp = dlmread(V_filename_G,'\t');
            [p,ia,ic] = unique(T_tmp(:,2));   
            T(:,1) = T_tmp(ia,1);
            T(:,2) = p;
            parameters(6).data = T;
            parameters(10).gemX = 1;
            parameters(10).sigX = V_cov_G;
            
            for jj = 1:2
                if jj == 1
                    parameters(10).type_v = 'DET     ';
                elseif jj == 2
                    parameters(10).type_v = 'NOR     ';
                end        

            % the run        
                [results, parametersOut] = Form(func, parameters, correlations, settings);
                %% 

                beta_values_G_GEV = results.beta;
                Pf_values_G_GEV = results.Pf;
                alfa_theta_G_GEV = parametersOut(1,1).alfa;
                design_theta_G_GEV = parametersOut(1,1).x;
                alfa_rho_G_GEV = parametersOut(2,1).alfa;
                design_rho_G_GEV = parametersOut(2,1).x;
                alfa_fy_G_GEV = parametersOut(3,1).alfa;
                design_fy_G_GEV = parametersOut(3,1).x;
                alfa_fc_G_GEV = parametersOut(4,1).alfa;
                design_fc_G_GEV = parametersOut(4,1).x;
                alfa_t_G_GEV = parametersOut(5,1).alfa;
                design_t_G_GEV = parametersOut(5,1).x;
                alfa_V_G_GEV = parametersOut(6,1).alfa;
                design_V_G_GEV = parametersOut(6,1).x;
                alfa_C_G_GEV = parametersOut(7,1).alfa;
                design_C_G_GEV = parametersOut(7,1).x;
                alfa_cr_G_GEV = parametersOut(8,1).alfa;
                design_cr_G_GEV= parametersOut(8,1).x;
                alfa_X_G_GEV = parametersOut(9,1).alfa;
                design_X_G_GEV= parametersOut(9,1).x;
                alfa_Sv_G_GEV = parametersOut(10,1).alfa;
                design_Sv_G_GEV= parametersOut(10,1).x;
                alfa_cd_G_GEV = parametersOut(11,1).alfa;
                design_cd_G_GEV= parametersOut(11,1).x;


                if extremes == 1
                alfas_mean_G_GEV_mom{jj} = [ alfa_V_G_GEV; alfa_C_G_GEV; alfa_cr_G_GEV; alfa_X_G_GEV; alfa_Sv_G_GEV ; alfa_cd_G_GEV; alfa_theta_G_GEV; alfa_fy_G_GEV; alfa_fc_G_GEV; alfa_t_G_GEV];
                beta_value_G_GEV_mom{jj} = beta_values_G_GEV;
                design_values_G_GEV_mom{jj} = [ design_theta_G_GEV; design_rho_G_GEV; design_fy_G_GEV; design_fc_G_GEV; design_t_G_GEV; design_V_G_GEV; design_C_G_GEV; design_cr_G_GEV; design_X_G_GEV; design_Sv_G_GEV ; design_cd_G_GEV];
                alfas_mean_G_GEV_squared_mom{jj} = transpose(alfas_mean_G_GEV_mom{jj}.^2);
%                 alfas_mean_GEV_squared_mom = transpose(alfas_mean_GEV_squared_mom);
                elseif extremes == 2
                alfas_mean_G_GEV_for{jj} = [ alfa_V_G_GEV; alfa_C_G_GEV; alfa_cr_G_GEV; alfa_X_G_GEV; alfa_Sv_G_GEV ; alfa_cd_G_GEV; alfa_theta_G_GEV; alfa_fy_G_GEV; alfa_fc_G_GEV; alfa_t_G_GEV];
                beta_value_G_GEV_for{jj} = beta_values_G_GEV;
                design_values_G_GEV_for{jj} = [ design_theta_G_GEV; design_rho_G_GEV; design_fy_G_GEV; design_fc_G_GEV; design_t_G_GEV; design_V_G_GEV; design_C_G_GEV; design_cr_G_GEV; design_X_G_GEV; design_Sv_G_GEV ; design_cd_G_GEV];
                alfas_mean_G_GEV_squared_for{jj} = transpose(alfas_mean_G_GEV_for{jj}.^2);
%                 alfas_mean_G_GEV_squared_for = transpose(alfas_mean_G_GEV_squared_for);
                end
                
            end
    end
end


%% Run Prob2B GEV - GEV 

if prob2b_GEV_GEV == 1
    for extremes = 1:2
        run('read_windspeed_dists_Schiphol');
        run('read_building_specific_input');
        % define parameters of distribution functions
        C_filename_GEV = [C_directory_GEV '\' C_filenames_GEV];
        T = [];
        T_tmp = dlmread(C_filename_GEV,'\t');
        [p,ia,ic] = unique(T_tmp(:,2));   
        T(:,1) = T_tmp(ia,1);
        T(:,2) = p;
        parameters(7).data = T;
        if extremes == 1
            func = @myZfunctionM2;
            % strength, blijft hetzelfde voor alle maxima

        elseif extremes == 2
            func = @myZfunctionQ2;
            % strength, blijft hetzelfde voor alle maxima

        end
      
        % roughness factor, blijft hetzelfde in elke limit state function
        parameters(8).gemX = crmean_squared;
        parameters(8).sigX = crstd_squared;
        parameters(8).x = parameters(4).gemX;

      
                
       
        V_filename_GEV = [V_directory_GEV '\' V_filenames_GEV];
        T = [];
        T_tmp = dlmread(V_filename_GEV,'\t');
        [p,ia,ic] = unique(T_tmp(:,2));   
        T(:,1) = T_tmp(ia,1);
        T(:,2) = p;
        parameters(6).data = T;                
        parameters(10).gemX = 1;
        parameters(10).sigX = V_cov_GEV;
                
        for jj = 1:2
            if jj == 1
                parameters(10).type_v = 'DET     ';
            elseif jj == 2
                parameters(10).type_v = 'NOR     ';
            end        
            [results, parametersOut] = Form(func, parameters, correlations, settings);
            %% 

            beta_values_GEV = results.beta;
            Pf_values_GEV = results.Pf;
            alfa_theta_GEV_GEV = parametersOut(1,1).alfa;
            design_theta_GEV_GEV = parametersOut(1,1).x;
            alfa_rho_GEV_GEV = parametersOut(2,1).alfa;
            design_rho_GEV_GEV = parametersOut(2,1).x;
            alfa_fy_GEV_GEV = parametersOut(3,1).alfa;
            design_fy_GEV_GEV = parametersOut(3,1).x;
            alfa_fc_GEV_GEV = parametersOut(4,1).alfa;
            design_fc_GEV_GEV = parametersOut(4,1).x;
            alfa_t_GEV_GEV = parametersOut(5,1).alfa;
            design_t_GEV_GEV = parametersOut(5,1).x;
            alfa_V_GEV_GEV = parametersOut(6,1).alfa;
            design_V_GEV_GEV = parametersOut(6,1).x;
            alfa_C_GEV_GEV = parametersOut(7,1).alfa;
            design_C_GEV_GEV = parametersOut(7,1).x;
            alfa_cr_GEV_GEV = parametersOut(8,1).alfa;
            design_cr_GEV_GEV= parametersOut(8,1).x;
            alfa_X_GEV_GEV = parametersOut(9,1).alfa;
            design_X_GEV_GEV= parametersOut(9,1).x;
            alfa_Sv_GEV_GEV = parametersOut(10,1).alfa;
            design_Sv_GEV_GEV= parametersOut(10,1).x;
            alfa_cd_GEV_GEV = parametersOut(11,1).alfa;
            design_cd_GEV_GEV= parametersOut(11,1).x;


        if extremes == 1
        alfas_mean_GEV_mom{jj} = [ alfa_V_GEV_GEV; alfa_C_GEV_GEV; alfa_cr_GEV_GEV; alfa_X_GEV_GEV; alfa_Sv_GEV_GEV ; alfa_cd_GEV_GEV; alfa_theta_GEV_GEV; alfa_fy_GEV_GEV; alfa_fc_GEV_GEV; alfa_t_GEV_GEV];
        beta_value_GEV_mom{jj} = beta_values_GEV;
        design_values_GEV_mom{jj} = [ design_theta_GEV_GEV; design_rho_GEV_GEV; design_fy_GEV_GEV; design_fc_GEV_GEV; design_t_GEV_GEV; design_V_GEV_GEV; design_C_GEV_GEV; design_cr_GEV_GEV; design_X_GEV_GEV; design_Sv_GEV_GEV ; design_cd_GEV_GEV];
        alfas_mean_GEV_squared_mom{jj} = transpose(alfas_mean_GEV_mom{jj}.^2);
%         alfas_mean_G_GEV_squared_mom = transpose(alfas_mean_G_GEV_squared_mom);
        elseif extremes == 2
        alfas_mean_GEV_for{jj} = [ alfa_V_GEV_GEV; alfa_C_GEV_GEV; alfa_cr_GEV_GEV; alfa_X_GEV_GEV; alfa_Sv_GEV_GEV ; alfa_cd_GEV_GEV; alfa_theta_GEV_GEV; alfa_fy_GEV_GEV; alfa_fc_GEV_GEV; alfa_t_GEV_GEV];
        beta_value_GEV_for{jj} = beta_values_GEV;
        design_values_GEV_for{jj} = [ design_theta_GEV_GEV; design_rho_GEV_GEV; design_fy_GEV_GEV; design_fc_GEV_GEV; design_t_GEV_GEV; design_V_GEV_GEV; design_C_GEV_GEV; design_cr_GEV_GEV; design_X_GEV_GEV; design_Sv_GEV_GEV ; design_cd_GEV_GEV];
        alfas_mean_GEV_squared_for{jj} = transpose(alfas_mean_GEV_for{jj}.^2);
%         alfas_mean_GEV_squared_for = transpose(alfas_mean_GEV_squared_for);
        end

        end
    end
        
        
end


%% 

% alfas_mean_G_squared_mom = alfas_mean_G_mom.^2;
% alfas_mean_G_squared_mom = transpose(alfas_mean_G_squared_mom);



% alfas_mean_G_squared_for = alfas_mean_G_for.^2;
% alfas_mean_G_squared_for = transpose(alfas_mean_G_squared_for);



% 
if both_extremes == 1
    if run_all == 1
        %     
%         alfas_mean = [alfas_mean_G_GEV_squared_mom{1}; alfas_mean_G_GEV_squared_mom{2}; alfas_mean_G_GEV_squared_for{1}; alfas_mean_G_GEV_squared_for{2}; alfas_mean_GEV_squared_mom{1}; alfas_mean_GEV_squared_mom{2}; alfas_mean_GEV_squared_for{1}; alfas_mean_GEV_squared_for{2}];
        %alfas_mean = [alfas_mean_G_squared_mom; alfas_mean_G_GEV_squared_mom; alfas_mean_GEV_squared_mom];
        alfas_mean = [alfas_mean_G_GEV_squared_mom{2}; alfas_mean_G_GEV_squared_for{2}; alfas_mean_GEV_squared_mom{2}; alfas_mean_GEV_squared_for{2}];


        % alfas_std = [std_alfa_G_mom; std_alfa_GEV_mom; std_alfa_G_for; std_alfa_GEV_for];  
        beta_values = [beta_value_G_GEV_mom{1}  beta_value_G_GEV_mom{2}; beta_value_G_GEV_for{1}  beta_value_G_GEV_for{2}; beta_value_GEV_mom{1}  beta_value_GEV_mom{2}; beta_value_GEV_for{1}  beta_value_GEV_for{2}];

        close all

        ax2 = [1 , 2 , 3 , 4 ];
        %ax2 = [1 , 2, 3];
        x = {'(1) M' '(1) Q' '(2) M' '(2) Q'};
        figure
        hold all
        bar_alfas = bar(ax2,alfas_mean,0.8,'stacked','FaceColor','flat');
%         b = bar(y,'FaceColor','flat');
        for k = 1:size(alfas_mean,2)
            bar_alfas(k).CData = k;
        end
            ylabel('Squared sensitivity factors $\alpha^{2}$','FontSize',12,'Interpreter','Latex')
            legend_alphas = legend({'$\alpha_{v_{pot}}$', '$\alpha_{c_X}$','$\alpha_{c_r}$','$\alpha_{\chi_{model}}$','$\alpha_{S_v}$','$\alpha_{\chi_{c_{d}}}$','$\alpha_{\chi_R}$','$\alpha_{f_y}$','$\alpha_{f_c}$'}, 'Location','southeastoutside','Interpreter','Latex');
            set(legend_alphas,'FontSize',16);
            set(gca,'YTick',0:0.1:1)
            set(gca,'XTick',ax2)
            set(gca,'XTickLabel',x)
            axis([0.5 4.5 0 1])
            set(gca,'fontsize',fontsize);
            fix_xticklabels(gca,0.1,{'FontSize',fontsize,'Interpreter','Latex'});
            figurename = sprintf('%s%s%s%s','figures/alphas_paper');
            print([figurename],'-dpng','-r600')
            savefig([figurename])

        figure
        hold all
        bar_betas = bar(ax2,beta_values(:,1:2),0.4,'FaceColor','flat');
        for k = 1:size(beta_values,2)
            bar_betas(k).CData = k;
        end
%         line([0.5 6.5],[3.8 3.8],'LineWidth', 1, 'Color', 'red', 'HandleVisibility','off');
            legend({'without $S_v$','with $S_v$'}, 'Location','northwest','Interpreter','Latex');
            ylabel('Reliability index $\beta$','Interpreter','Latex')
            %legend({'$\alpha_V$', '$\alpha_C$','$\alpha_{cr}$','$\alpha_X$','$\alpha_{Sv}$','$\alpha_{c_{d}}$','$\alpha_R$'}, 'Location','southwest','Interpreter','Latex');
            x1 = 0.8;
            y1 = 4.0;
%             txt1 = '$\beta_{\mathrm{target}}$';
%             text(x1,y1,txt1,'FontSize',fontsize,'Interpreter','Latex')
            set(gca,'YTick',0:0.5:5.5)
            set(gca,'XTick',ax2)
            set(gca,'XTickLabel',x)
            grid minor
            axis([0.5 4.5 0 5.5])
            set(gca,'fontsize',fontsize);
            fix_xticklabels(gca,0.1,{'FontSize',fontsize,'Interpreter','Latex'});
            figurename = sprintf('%s%s%s%s','figures/betas_paper');
            print([figurename],'-dpng','-r600')
            savefig([figurename])
    else
    end
    else
end