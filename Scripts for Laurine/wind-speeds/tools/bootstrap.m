% Bootstrapping function
%
% population : original dataset
% n : number of bootstrap samples
% tail : fixed tail index from weather model dataset

function [GEVparameters] = bootstrap(population, n, tail)

% Initialization
samples = zeros(length(population), n);
parmhat = zeros(n, 3);  % [k, sigma, mu] parameters of GEV

% if the tail is not in the input
if nargin < 3 || isempty(tail)
    % Generate bootstrap samples and fit GEV distribution to each
    for i = 1:n
        % Generate a bootstrap sample with replacement
        samples(:, i) = datasample(population, length(population), 'Replace', true);
        
        % Fit GEV distribution to the bootstrap sample
        % gevfit(X) returns maximum likelihood estimates
        parmhat(i, :) = gevfit(samples(:, i));
    end
    
    GEVparameters = array2table(parmhat, "VariableNames",{'tail', 'scale', 'location'});
    
    % Display the GEV parameters for each bootstrap sample
    disp('GEV Parameters for each bootstrap sample :');
    disp(['Result for ', num2str(n), ' estimations'])
    disp(['ndraw = ', num2str(length(population)), 'years'])
    disp(GEVparameters);

else

% Generate bootstrap samples and fit GEV distribution to each
for i = 1:n
    % Generate a bootstrap sample with replacement
    samples(:, i) = datasample(population, length(population), 'Replace', true);
    
    % Fit GEV distribution to the bootstrap sample
    % gevfit(X) returns maximum likelihood estimates
    parmhat(i, :) = gevfit_fixedtail(samples(:, i),tail);
end

GEVparameters = array2table(parmhat, "VariableNames",{'scale', 'location'});

% Display the GEV parameters for each bootstrap sample
disp('GEV Parameters for each bootstrap sample :');
disp(['Result for ', num2str(n), ' estimations'])
disp(['ndraw = ', num2str(length(population)), 'years'])
disp(GEVparameters);

end

