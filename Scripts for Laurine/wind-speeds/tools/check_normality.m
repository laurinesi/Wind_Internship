% Plot histogram of bootstrap estimates and compare/fit to a normal
% distribution

function [Qt, Qs, Ql] = check_normality(GEVparameters, parmhat, dist, binctrs)

% Qt            : quantiles of the tail for the cumulative 
%                 probabilities 0.025 and 0.975
% Qs            : quantiles of the scale for the cumulative 
%                 probabilities 0.025 and 0.975
% Ql            : quantiles of the location for the cumulative 
%                 probabilities 0.025 and 0.975
% GEVparameters : estimates of GEV dist using MLE and bootstrap method
% dist          : if a normal distribution is plotted on top of histograms or not
% binctrs       : if the bin centers are circled (easier to see the normality)



% Check if the tail is fixed
num_cols = size(GEVparameters, 2);


if num_cols == 3

    % Calling of ML estimates
    tail = parmhat(1);
    scale_model = parmhat(2);
    location_model = parmhat(3);

    % Compute the quantiles 
    Qt = quantile(GEVparameters.tail,[0.025 0.975]);
    Qs = quantile(GEVparameters.scale,[0.025 0.975]);
    Ql = quantile(GEVparameters.location,[0.025 0.975]);

    % Tail
    figure
    h = histogram(GEVparameters.tail,"LineWidth",1,"BinMethod","scott");
    set(h,'FaceColor',[.98 .98 .98],'EdgeColor', "#7600bc");
    % Compute mean and standard deviation of the bootstrap data
    mu = mean(GEVparameters.tail);
    sigma = std(GEVparameters.tail);
    % Generate x values for plotting the normal distribution
    x = linspace(mu - 4*sigma, mu + 4*sigma, 100);
%    x = linspace(min(GEVparameters.tail), max(GEVparameters.tail), 100);
    % Generate the normal distribution values (pdf) for the x range
    pdf_normal = normpdf(x, mu, sigma);
    % Scale the normal distribution to match the histogram counts
    % (The area under the histogram and the normal distribution should match)
    pdf_normal_scaled = pdf_normal * (h.BinWidth * sum(h.BinCounts));
    if binctrs == 1
        set(h,'FaceColor',[.98 .98 .98],'EdgeColor',[.94 .94 .94]);
        hold on
        plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts, "Marker","o", "LineStyle","none", "LineWidth",1)
%         % draw 95% confidence interval
%         hold on
%         % Draw a vertical line at q1 for 95% confidence interval
%         plot([Qt(1), Qt(1)], [0, mean(pdf_normal_scaled)], 'b', 'LineWidth', 1);
%         hold on
%         % Draw a vertical line at q2 for 95% confidence interval
%         plot([Qt(2), Qt(2)], [0, mean(pdf_normal_scaled)], 'b', 'LineWidth', 1);
    end
    if dist == 1
        hold on
        % Plot the normal distribution on top of the histogram
        plot(x, pdf_normal_scaled, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Normal PDF')
%         hold on
%         % Draw a vertical line at the mean of bootstrapped estimate
%         plot([mu, mu], [0, max(pdf_normal_scaled)], 'k-', 'LineWidth', 0.6);
%         hold on
%         % Draw a vertical line at the ML estimate
%         plot([tail, tail], [0, max(pdf_normal_scaled)], 'b-', 'LineWidth', 0.6);
    end
    xlim([min(x), max(x)])
    title('PDF of Shape Parameter');
    xlabel('Bootstrap tail estimates'); 
    ylabel('Frequency');
    
    
    % Scale
    figure
    h = histogram(GEVparameters.scale,"LineWidth",1,"BinMethod","scott");
    set(h,'FaceColor',[.98 .98 .98],'EdgeColor',"#41a317");
    % Compute mean and standard deviation of the bootstrap data
    mu = mean(GEVparameters.scale);
    sigma = std(GEVparameters.scale);
    % Generate x values for plotting the normal distribution
    x = linspace(mu - 4*sigma, mu + 4*sigma, 100);
%    x = linspace(min(GEVparameters.scale), max(GEVparameters.scale), 100);
    % Generate the normal distribution values (pdf) for the x range
    pdf_normal = normpdf(x, mu, sigma);
    % Scale the normal distribution to match the histogram counts
    % (The area under the histogram and the normal distribution should match)
    pdf_normal_scaled = pdf_normal * (h.BinWidth * sum(h.BinCounts));
    if binctrs == 1
        set(h,'FaceColor',[.98 .98 .98],'EdgeColor',[.94 .94 .94]);
        hold on
        plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts, "Marker","o", "LineStyle","none", "LineWidth",1)
    %     % draw 95% confidence interval
    %     hold on
    %     % Draw a vertical line at q1 for 95% confidence interval
    %     plot([Qs(1), Qs(1)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    %     hold on
    %     % Draw a vertical line at q2 for 95% confidence interval
    %     plot([Qs(2), Qs(2)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    end
    if dist == 1
        hold on
        % Plot the normal distribution on top of the histogram
        plot(x, pdf_normal_scaled, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Normal PDF')
%         hold on
%         % Draw a vertical line at the mean
%         plot([mu, mu], [0, max(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    end
    xlim([min(x), max(x)])
    title('PDF of Scale Parameter');
    xlabel('Bootstrap scale estimates'); 
    ylabel('Frequency');
    
    
    % Location
    figure
    h = histogram(GEVparameters.location,"LineWidth",1,"BinMethod","scott");
    set(h,'FaceColor',[.98 .98 .98],'EdgeColor',"#357ec7");
    % Compute mean and standard deviation of the bootstrap data
    mu = mean(GEVparameters.location);
    sigma = std(GEVparameters.location);
    % Generate x values for plotting the normal distribution
    x = linspace(mu - 4*sigma, mu + 4*sigma, 100);
%    x = linspace(min(GEVparameters.location), max(GEVparameters.location), 100);
    % Generate the normal distribution values (pdf) for the x range
    pdf_normal = normpdf(x, mu, sigma);
    % Scale the normal distribution to match the histogram counts
    % (The area under the histogram and the normal distribution should match)
    pdf_normal_scaled = pdf_normal * (h.BinWidth * sum(h.BinCounts));
    if binctrs == 1
        set(h,'FaceColor',[.98 .98 .98],'EdgeColor',[.94 .94 .94]);
        hold on
        plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts, "Marker","o", "LineStyle","none", "LineWidth",1)
    %     % draw 95% confidence interval
    %     hold on
    %     % Draw a vertical line at q1 for 95% confidence interval
    %     plot([Ql(1), Ql(1)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    %     hold on
    %     % Draw a vertical line at q2 for 95% confidence interval
    %     plot([Ql(2), Ql(2)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    end
    if dist == 1
        hold on
        % Plot the normal distribution on top of the histogram
        plot(x, pdf_normal_scaled, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Normal PDF')
    end
    xlim([min(x), max(x)])
    title('PDF of Location Parameter');
    xlabel('Bootstrap location estimates'); 
    ylabel('Frequency');
    

    % normal probability plot
    figure
    normplot(GEVparameters.tail)
    title('Normal Probability Plot (shape)');
    figure
    normplot(GEVparameters.scale)
    title('Normal Probability Plot (scale)');
    figure
    normplot(GEVparameters.location)
    title('Normal Probability Plot (location)');

elseif num_cols == 2

    % Compute the quantiles 
    Qt = 1;
    Qs = quantile(GEVparameters.scale,[0.025 0.975]);
    Ql = quantile(GEVparameters.location,[0.025 0.975]);
    

    % Scale
    figure
    h = histogram(GEVparameters.scale,"LineWidth",1,"BinMethod","scott");
    set(h,'FaceColor',[.98 .98 .98],'EdgeColor',"#41a317");
    % Compute mean and standard deviation of the bootstrap data
    mu = mean(GEVparameters.scale);
    sigma = std(GEVparameters.scale);
    % Generate x values for plotting the normal distribution
    x = linspace(mu - 4*sigma, mu + 4*sigma, 100);
%     x = linspace(min(GEVparameters.scale), max(GEVparameters.scale), 100);
    % Generate the normal distribution values (pdf) for the x range
    pdf_normal = normpdf(x, mu, sigma);
    % Scale the normal distribution to match the histogram counts
    % (The area under the histogram and the normal distribution should match)
    pdf_normal_scaled = pdf_normal * (h.BinWidth * sum(h.BinCounts));
    if binctrs == 1
        set(h,'FaceColor',[.98 .98 .98],'EdgeColor',[.94 .94 .94]);
        hold on
        plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts, "Marker","o", "LineStyle","none", "LineWidth",1)
    %     % draw 95% confidence interval
    %     hold on
    %     % Draw a vertical line at q1 for 95% confidence interval
    %     plot([Qs(1), Qs(1)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    %     hold on
    %     % Draw a vertical line at q2 for 95% confidence interval
    %     plot([Qs(2), Qs(2)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    end
    if dist == 1
        hold on
        % Plot the normal distribution on top of the histogram
        plot(x, pdf_normal_scaled, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Normal PDF')
    end
    xlim([min(x), max(x)])
    title('PDF of Scale Parameter');
    xlabel('Bootstrap scale estimates'); 
    ylabel('Frequency');
    
    
    % Location
    figure
    h = histogram(GEVparameters.location,"LineWidth",1,"BinMethod","scott");
    set(h,'FaceColor',[.98 .98 .98],'EdgeColor',"#357ec7");
    % Compute mean and standard deviation of the bootstrap data
    mu = mean(GEVparameters.location);
    sigma = std(GEVparameters.location);
    % Generate x values for plotting the normal distribution
    x = linspace(mu - 4*sigma, mu + 4*sigma, 100);
%     x = linspace(min(GEVparameters.location), max(GEVparameters.location), 100);
    % Generate the normal distribution values (pdf) for the x range
    pdf_normal = normpdf(x, mu, sigma);
    % Scale the normal distribution to match the histogram counts
    % (The area under the histogram and the normal distribution should match)
    pdf_normal_scaled = pdf_normal * (h.BinWidth * sum(h.BinCounts));
    if binctrs == 1
        set(h,'FaceColor',[.98 .98 .98],'EdgeColor',[.94 .94 .94]);
        hold on
        plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts, "Marker","o", "LineStyle","none", "LineWidth",1)
    %     % draw 95% confidence interval
    %     hold on
    %     % Draw a vertical line at q1 for 95% confidence interval
    %     plot([Ql(1), Ql(1)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    %     hold on
    %     % Draw a vertical line at q2 for 95% confidence interval
    %     plot([Ql(2), Ql(2)], [0, mean(pdf_normal_scaled)], 'k--', 'LineWidth', 1);
    end
    if dist == 1
        hold on
        % Plot the normal distribution on top of the histogram
        plot(x, pdf_normal_scaled, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Normal PDF')
    end
    xlim([min(x), max(x)])
    title('PDF of Location Parameter');
    xlabel('Bootstrap location estimates'); 
    ylabel('Frequency');
    

    % normal probability plot
    figure
    normplot(GEVparameters.scale)
    title('Normal Probability Plot (scale)');
    figure
    normplot(GEVparameters.location)
    title('Normal Probability Plot (location)');
end




