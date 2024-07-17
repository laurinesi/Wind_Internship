function linplot(settings, parameters, parametersOut)
%
% parameters
itmax = settings.itmax;
% make dummy grid
x1_grid = linspace(0,40,itmax);
x2_grid = linspace(0,40,itmax);
% step width
delta_x1 = range(parametersOut(1,1).x)/(itmax-1);
delta_x2 = range(parametersOut(2,1).x)/(itmax-1);

% generate a 2D meshgrid
[x1_mesh,x2_mesh] = meshgrid(x1_grid,x2_grid);
% evaluate pdf
f_x1 = normpdf(x1_mesh,parameters(1,1).gemX,parameters(1,1).sigX);
f_x2 = normpdf(x2_mesh, parameters(2,1).gemX,parameters(2,1).sigX);

f_x1x2 = f_x1 .* f_x2;
% identify failure domain where R<S
idx = x1_mesh <= x2_mesh;
Pf_joint = sum(f_x1x2(idx).*delta_x1.*delta_x2);

% plot the joint pdf
figure
%surf(x1_mesh, x2_mesh, f_x1x2)
plot(f_x1,f_x2)
xlabel('x_1')
ylabel('x_2')
zlabel('pdf')
title('joint pdf X_1 and X_2')

