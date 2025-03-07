function f_tf_histogram(selectedTrialsMatrix, minFreqHz, maxFreqHz, sampleRateReduction, actualSampleRate, stDevCycles, normalizeFlag, headerInfo, selectedChannel)

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