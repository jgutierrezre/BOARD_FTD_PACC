function f_tf_histogram(selectedTrialsMatrix, minFreqHz, maxFreqHz, sampleRateReduction, actualSampleRate, stDevCycles, normalizeFlag, headerInfo, selectedChannel)
    % f_tf_histogram
    % Computes time-frequency decompositions for a set of trials using Morse wavelets,
    % aggregates frequency-band-specific magnitudes, and plots histograms per EEG band.
    %
    % Inputs:
    %   selectedTrialsMatrix   : [time x trials] matrix of extracted trials
    %   minFreqHz              : minimum frequency for decomposition (Hz)
    %   maxFreqHz              : maximum frequency for decomposition (Hz)
    %   sampleRateReduction    : downsampling factor (to reduce time resolution)
    %   actualSampleRate       : original sampling rate (Hz)
    %   stDevCycles            : wavelet width in standard deviations (controls time-frequency tradeoff)
    %   normalizeFlag          : 1 to apply log-transform and z-score across trials, 0 to use raw magnitudes
    %   headerInfo             : structure containing channel names (for figure titles)
    %   selectedChannel        : index of the channel being analyzed
    %
    % Notes:
    %   - Applies the Morse analytic wavelet transform across trials.
    %   - Magnitude values are averaged within canonical EEG bands:
    %     Delta (0.5–4 Hz), Theta (4–8 Hz), Alpha (8–13 Hz), Beta (13–30 Hz),
    %     Gamma (30–120 Hz), and High-Frequency Oscillations (120+ Hz).
    %   - Histogram of average band-specific magnitude is plotted per band.
    %
    % Example:
    %   f_tf_histogram(trialsMatrix, 2, 120, 2, 1000, 4, 1, headerInfo, 1);

    % Adjust sampling rate and select single channel
    actualSampleRate = actualSampleRate / sampleRateReduction;
    selectedTrialsMatrix = downsample(selectedTrialsMatrix, sampleRateReduction);

    %% Configurable Parameters
    EEG_BANDS = struct( ...
        'Delta', [0.5, 4], ...
        'Theta', [4, 8], ...
        'Alpha', [8, 13], ...
        'Beta', [13, 30], ...
        'Gamma', [30, 120], ...
        'HFO', [120, maxFreqHz] ...
    );

    bandNames = fieldnames(EEG_BANDS);
    numBands = length(bandNames);

    % Frequency segmentation width (configurable scaling factor)
    SEGMENT_WIDTH_SCALE = 2;
    freqSegmentation = SEGMENT_WIDTH_SCALE * (maxFreqHz - minFreqHz);

    % Compute Morse wavelet transform for selected channel
    [tfData, ~, freqVector] = f_MorseAWTransformMatlab(selectedTrialsMatrix, actualSampleRate, ...
        minFreqHz, maxFreqHz, freqSegmentation, stDevCycles, ...
        1, 0, 0, 0, []); % Only compute magnitudes (no power, no phase)

    % Normalize if required
    if normalizeFlag == 1
        % Step 1: Log-transform to compress large variations
        tfData = 10 * log10(abs(tfData) + eps);

        % Step 2: Z-score normalization across trials
        tfData = f_Mat_to_zscore(tfData);
    else
        % Use raw magnitudes without transformation
        tfData = abs(tfData);
    end

    % Compute average magnitude per EEG band
    tfAverage = zeros(numBands, size(tfData, 2)); % Time dimension

    for j = 1:numBands
        bandLimits = EEG_BANDS.(bandNames{j});
        bandIndices = find(freqVector >= bandLimits(1) & freqVector <= bandLimits(2));

        if ~isempty(bandIndices)
            tfAverage(j, :) = mean(tfData(bandIndices, :), 1); % Average across frequencies
        end

    end

    % Plot histograms for each EEG band
    figure;

    for j = 1:numBands
        subplot(numBands, 1, j);
        histogram(tfAverage(j, :), 'Normalization', 'probability');
        % histogram(log10(tfAverage(j, :)), 'Normalization', 'probability');
        xlabel('Magnitude');
        ylabel('Probability');
        title([bandNames{j}, ' Band']);
    end

    sgtitle(headerInfo.recChNames(selectedChannel));
end
