% Bootstrapping 

function [GEVparameters] = bootstrap(population, n)

% Initialization
samples = zeros(length(population), n);
parmhat = zeros(n, 3);  % [k, sigma, mu] parameters of GEV

% Generate bootstrap samples and fit GEV distribution to each
for i = 1:n
    % Generate a bootstrap sample with replacement
    samples(:, i) = datasample(population, length(population), 'Replace', true);
    
    % Fit GEV distribution to the bootstrap sample
    % gevfit(X) returns maximum likelihood estimates
    parmhat(i, :) = gevfit(samples(:, i));
end

GEVparameters = array2table(parmhat, "VariableNames",{'shape', 'scale', 'location'});

% Display the GEV parameters for each bootstrap sample
disp('GEV Parameters for each bootstrap sample :');
disp(['Result for ', num2str(n), ' estimations'])
disp(['ndraw = ', num2str(length(population)), 'years'])
disp(GEVparameters);