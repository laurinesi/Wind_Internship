        
if bootstrapping == 2

    dim1 = Nreps+1 ;
else
    dim1 = 1 ;
end

X_t_analysis_mom = cell(1,3600);
X_t_analysis_for = cell(1,3600);
X_t_analysis_stat = cell(1,3600);
id = cell(1,3600);
id_stat = cell(1,3600);
ymin_transformed = cell(1,plotFs);
ymax_transformed = cell(1,plotFs);

if extremes == 1
F_N_Tfs_mom = cell(1,3600);
x_mom = cell(1,3600);
F_N_T3600_mom = cell(1,3600);
parGEVmom = cell(1,dim1);
Fx_GEV_MLE_Tfs_mom = cell(1,3600);
Fx_G_MLE_Tfs_mom = cell(1,3600);
Fx_GEV_MLE_T3600_mom = cell(1,3600);
Fx_G_MLE_T3600_mom = cell(1,3600);
x_fractiles_GEV_MLE_mom_Tfs = cell(1,3600);
x_fractiles_G_MLE_mom_Tfs = cell(1,3600);
x_fractiles_GEV_MLE_mom = cell(1,3600);
x_fractiles_G_MLE_mom = cell(1,3600);

[Fx_Tfs_mom{1:plotFs-1}] = deal([]);
Fx_Tfs_mom{plotFs} = deal(cell(1,b));

[Fx_T3600_mom{1:plotFs-1}] = deal([]);
Fx_T3600_mom{plotFs} = deal(cell(1,b));

x_fractiles_F_N_mom = cell(1,plotFs);
x_cook_F_N_mom = zeros(1,plotFs);
F_N_Tfs_mom_transformed = cell(1,plotFs);
F_N_T3600_mom_transformed = cell(1,plotFs);


elseif extremes == 2

F_N_Tfs_for = cell(1,3600);
x_for = cell(1,3600);
F_N_T3600_for = cell(1,3600);
parGEVfor = cell(1,dim1);
Fx_GEV_MLE_Tfs_for = cell(1,3600);
Fx_G_MLE_Tfs_for = cell(1,3600);
Fx_GEV_MLE_T3600_for = cell(1,3600);
Fx_G_MLE_T3600_for = cell(1,3600);
x_fractiles_GEV_MLE_for_Tfs = cell(1,3600);
x_fractiles_G_MLE_for_Tfs = cell(1,3600);
x_fractiles_GEV_MLE_for = cell(1,3600);
x_fractiles_G_MLE_for = cell(1,3600);

[Fx_Tfs_for{1:plotFs-1}] = deal([]);
Fx_Tfs_for{plotFs} = deal(cell(1,b));

[Fx_T3600_for{1:plotFs-1}] = deal([]);
Fx_T3600_for{plotFs} = deal(cell(1,b));

x_fractiles_F_N_for = cell(1,plotFs);
x_cook_F_N_for = zeros(1,plotFs);
F_N_Tfs_for_transformed = cell(1,plotFs);
F_N_T3600_for_transformed = cell(1,plotFs);

elseif extremes == 3

F_N_Tfs_stat = cell(1,3600);
x_stat = cell(1,3600);
F_N_T3600_stat = cell(1,3600);
parGEVstat = cell(1,dim1);
Fx_GEV_MLE_Tfs_stat = cell(1,3600);
Fx_G_MLE_Tfs_stat = cell(1,3600);
Fx_GEV_MLE_T3600_stat = cell(1,3600);
Fx_G_MLE_T3600_stat = cell(1,3600);
x_fractiles_GEV_MLE_stat_Tfs = cell(1,3600);
x_fractiles_G_MLE_stat_Tfs = cell(1,3600);
x_fractiles_GEV_MLE_stat = cell(1,3600);
x_fractiles_G_MLE_stat = cell(1,3600);


end

xmin = cell(1,plotFs);
xmax = cell(1,plotFs);

for i = T_fs
    X_t_analysis_mom{i} = zeros(dim1,size(Cm_max{i},2));
    X_t_analysis_for{i} = zeros(dim1,size(Cf_max{i},2));
    X_t_analysis_stat{i} = zeros(dim1,size(Cm_stat{i},2));
    id{i} = zeros(Nreps,size(Cm_max{i},2));
    id_stat{i} = zeros(Nreps,size(Cm_stat{i},2));
    
    if extremes == 1
        
        mean_X_mom = zeros(dim1,3600);
        std_X_mom = zeros(dim1,3600);
        a_X_mom = zeros(dim1,3600);
        F_N_Tfs_mom{i} = zeros(dim1,size(Cm_max{i},2));
        F_N_T3600_mom{i} = zeros(dim1,size(Cm_max{i},2));
        alpha1_MLE_mom = zeros(dim1,3600);
        u1_MLE_mom = zeros(dim1,3600);
        nlogL_G_MLE_mom = zeros(1,3600);
        
    elseif extremes == 2

        mean_X_for = zeros(dim1,3600);
        std_X_for = zeros(dim1,3600);
        a_X_for = zeros(dim1,3600);
        F_N_Tfs_for{i} = zeros(dim1,size(Cm_max{i},2));
        F_N_T3600_for{i} = zeros(dim1,size(Cm_max{i},2));
        alpha1_MLE_for = zeros(dim1,3600);
        u1_MLE_for = zeros(dim1,3600);
        nlogL_G_MLE_for = zeros(1,3600);

    elseif extremes == 3 
        
        mean_X_stat = zeros(dim1,3600);
        std_X_stat = zeros(dim1,3600);
        a_X_stat = zeros(dim1,3600);
        F_N_Tfs_stat{i} = zeros(dim1,size(Cm_max{i},2));
        F_N_T3600_stat{i} = zeros(dim1,size(Cm_max{i},2));
        alpha1_MLE_stat = zeros(dim1,3600);
        u1_MLE_stat = zeros(dim1,3600);
        nlogL_G_MLE_stat = zeros(1,3600);
        
    end
    
end