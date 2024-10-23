% Stability of tail index over all heights
% sample fraction p = 10%

load("..\..\..\Alex's code\out\tail_values.mat")

% Extract list
p = tail.p;
val = tail.val;
lb = tail.lb;
ub = tail.ub;
groups = unique(tail.group);

% Initialization
tail_index_at_p_0_1 = [];


% Loop to find the tail index corresponding to p = 0.1
for i = 1:length(groups)
    group = groups(i);
    
    % Filter data for the current group
    group_indices = strcmp(tail.group, group);
    group_data = val(group_indices);
    group_proba = p(group_indices);
    group_lb = lb(group_indices);
    group_ub = ub(group_indices);

%     disp(['Processing group: ', group{1,1}])
    
    % Find the index corresponding to p = 0.1
    k = find(abs(group_proba-0.1) < 0.001);

    % Find the tail index value at that index
    tail_indices_at_p_0_1(i,1) = group_data(k);

    % Find the confidence interval
    lb_0_1(i) = group_lb(k);
    ub_0_1(i) = group_ub(k);

    disp(['Tail index at p = 0.1 for group ', group{1,1}, ': ', num2str(tail_indices_at_p_0_1(i))]);
end

% Put in order (10, 50, 100, 150, 200, 250, 300)
for j = 1:length(groups)-1
    % Fix tail index
    tail_index_at_p_0_1(1) = tail_indices_at_p_0_1(7);
    tail_index_at_p_0_1(j+1,1) = tail_indices_at_p_0_1(j);
    % Fix lower bound
    fixed_lb_0_1(1) = lb_0_1(7);
    fixed_lb_0_1(j+1,1) = lb_0_1(j);
    % Fix upper bound
    fixed_ub_0_1(1) = ub_0_1(7);
    fixed_ub_0_1(j+1,1) = ub_0_1(j);
end

% Calculate the mean of the tail indices for p = 0.1 across all groups
mean_tail_index = mean(tail_index_at_p_0_1);
std_tail_index = std(tail_index_at_p_0_1);
cov_tail_index = std_tail_index / mean_tail_index;

% Display the result
disp('Mean tail index at p = 0.1 across all groups:');
disp(mean_tail_index);
disp('Std tail index at p = 0.1 across all groups:');
disp(std_tail_index);
disp('COV tail index at p = 0.1 across all groups:');
disp(cov_tail_index);
%%
h = [10; 50; 100; 150; 200; 250; 300];
ci = [];

ci(:,1) = h;
ci(:,2) = fixed_lb_0_1;
ci(:,3) = fixed_ub_0_1;

% Create an interpolated set of values
numPoints = 100; % Number of interpolation points
xi = linspace(min(h), max(h), numPoints); % Interpolated x values
yi1 = interp1(h, fixed_lb_0_1, xi, 'linear'); % Interpolated y1 values
yi2 = interp1(h, fixed_ub_0_1, xi, 'linear'); % Interpolated y2 values

figure
plot(h, tail_index_at_p_0_1, 'k',LineStyle="-.", Marker="+", DisplayName='iHilli estimation',MarkerSize=7)
hold on
plot(h, fixed_lb_0_1,'b', LineWidth=1, HandleVisibility='off')
hold on
plot(h, fixed_ub_0_1,'b', LineWidth=1, HandleVisibility='off')
hold on
fill([xi fliplr(xi)], [yi1 fliplr(yi2)], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none', DisplayName='Confidence Interval');
% hold on
% yline(0.415,'w-', LineWidth=1, HandleVisibility='off')
xlim([10 300])
xlabel('Height (m)');
ylabel('Tail index');
title('Tail index over different heights, th=0.3, p=0.1');
legend;
grid minor;
