function [bootStats, ci, pVal] = f_BootstrapTrialStats(trialsMatrixA, numBoots, statFunc, trialsMatrixB)
    % f_BootstrapTrialStats
    % Performs bootstrap resampling on trial-based electrophysiological data
    % to estimate the variability of a statistic, and optionally compare two conditions.
    %
    % Inputs:
    %   trialsMatrixA : [time x trials] matrix of trials for condition A
    %   numBoots      : number of bootstrap iterations
    %   statFunc      : function handle to compute statistic (e.g., @mean, @(x) mean(abs(x)))
    %   trialsMatrixB : (optional) [time x trials] matrix of trials for condition B
    %
    % Outputs:
    %   bootStats     : [numBoots x 1] bootstrapped statistics or differences
    %   ci            : [1 x 2] 95% confidence interval (2.5th, 97.5th percentiles)
    %   pVal          : two-sided p-value (if B is provided); NaN otherwise
    %
    % Example:
    %   [stats, ci, p] = f_BootstrapTrialStats(trialsA, 1000, @mean, trialsB);

    % Defaults
    if nargin < 3 || isempty(statFunc)
        statFunc = @(x) mean(x, 'all'); % Mean across time and trials
    end

    if nargin < 4
        trialsMatrixB = [];
    end

    nA = size(trialsMatrixA, 2);
    bootStats = zeros(numBoots, 1);

    if isempty(trialsMatrixB)
        % --- Single group: estimate confidence interval ---
        for b = 1:numBoots
            idxA = randi(nA, [1, nA]); % resample with replacement
            bootSample = trialsMatrixA(:, idxA);
            bootStats(b) = statFunc(bootSample);
        end

        ci = prctile(bootStats, [2.5, 97.5]);
        pVal = NaN;

    else
        % --- Two groups: compare statistic difference ---
        nB = size(trialsMatrixB, 2);
        bootStats = zeros(numBoots, 1);

        for b = 1:numBoots
            idxA = randi(nA, [1, nA]);
            idxB = randi(nB, [1, nB]);
            bootA = trialsMatrixA(:, idxA);
            bootB = trialsMatrixB(:, idxB);

            statA = statFunc(bootA);
            statB = statFunc(bootB);
            bootStats(b) = statA - statB;
        end

        ci = prctile(bootStats, [2.5, 97.5]);
        pVal = 2 * min(mean(bootStats > 0), mean(bootStats < 0)); % two-sided
    end

end
