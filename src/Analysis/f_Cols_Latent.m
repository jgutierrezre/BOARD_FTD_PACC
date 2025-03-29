function [trialsMatrix, trialTimeVector] = f_Cols_Latent(dataMatrix, selectedChannel, threshold, trialWindowLength, exclusionWindowLength, numChannels, actualSampleRate, numTrials)
    % f_Cols_Latent
    % Extracts clean, non-stimulus (latent) trials from continuous electrophysiological data,
    % avoiding regions near stimulus events, and plots excluded zones and selected trials.
    %
    % Inputs:
    %   dataMatrix           : [time x channels] continuous data matrix
    %   selectedChannel      : index of the data channel to extract latent activity from
    %   threshold            : threshold to detect stimulus events (based on numChannels)
    %   trialWindowLength    : duration of each trial in seconds
    %   exclusionWindowLength: time to exclude before and after each stimulus (in seconds)
    %   numChannels          : index of the stimulus-detection channel
    %   actualSampleRate     : sampling rate of the data (Hz)
    %   numTrials            : number of random, non-overlapping latent trials to extract
    %
    % Outputs:
    %   trialsMatrix         : [samples x trials] matrix of extracted latent trials
    %   trialTimeVector      : [samples x 1] time vector (0 to trialWindowLength)
    %
    % Example:
    %   [trials, tVec] = f_Cols_Latent(data, 1, 0.6, 1.0, 0.2, 16, 1000, 50);

    %% --- Parameters & Setup ---
    threshold = abs(threshold);
    trialLength = round(trialWindowLength * actualSampleRate);
    exclusionLength = round(exclusionWindowLength * actualSampleRate);
    totalSamples = size(dataMatrix, 1);
    signal = dataMatrix(:, selectedChannel);
    timeVec = (0:totalSamples - 1) / actualSampleRate;

    %% --- Step 1: Detect stimulus onsets (abs threshold, both polarities) ---
    stimMask = abs(dataMatrix(:, numChannels)) > threshold;
    stimDiff = diff([0; stimMask]);
    onsetIndices = find(stimDiff == 1); % entering the threshold zone

    %% --- Step 2: Build exclusion mask around each stimulus ---
    excludeMask = false(totalSamples, 1);
    excludeRegions = zeros(length(onsetIndices), 2);

    for i = 1:length(onsetIndices)
        i1 = max(onsetIndices(i) - exclusionLength, 1);
        i2 = min(onsetIndices(i) + exclusionLength, totalSamples);
        excludeMask(i1:i2) = true;
        excludeRegions(i, :) = [i1, i2];
    end

    % Step 3: Build valid trial start vector
    valid = ~excludeMask;
    valid(end - trialLength + 1:end) = false; % prevent overflow at end

    % Find all valid trial start indices where the entire window is clean
    candidateStarts = find(movsum(valid, trialLength, 'Endpoints', 'discard') == trialLength);

    % Enforce non-overlapping selection
    trialStarts = zeros(1, min(numTrials, numel(candidateStarts)));
    taken = false(size(valid));
    count = 0;
    shuffledOrder = candidateStarts(randperm(numel(candidateStarts)));
    cursor = 1;

    while count < numTrials && cursor <= numel(shuffledOrder)
        ix = shuffledOrder(cursor);

        if ~any(taken(ix:ix + trialLength - 1))
            count = count + 1;
            trialStarts(count) = ix;
            taken(ix:ix + trialLength - 1) = true; % mark window as used
        end

        cursor = cursor + 1;
    end

    trialStarts = trialStarts(1:count);

    if count < numTrials
        warning('Only %d non-overlapping trials could be extracted.', count);
        numTrials = count;
    end

    %% --- Step 4: Extract trial windows ---
    trialsMatrix = zeros(trialLength, numTrials);

    for p = 1:numTrials
        ix = trialStarts(p);
        trialsMatrix(:, p) = signal(ix:ix + trialLength - 1);
    end

    %% --- Step 5: Create trial-local time vector ---
    trialTimeVector = linspace(0, trialWindowLength, trialLength)';

    %% --- Step 6: Plotting ---
    latentSignal = signal(~excludeMask);
    latentTime = linspace(0, length(latentSignal) / actualSampleRate, length(latentSignal));

    figure;

    % --- Top subplot: raw signal with shaded exclusions AND trials
    subplot(2, 1, 1); hold on;
    plot(timeVec, signal, 'Color', [0.4 0.4 1], 'DisplayName', 'Raw signal');

    yLimits = ylim;

    % 🔴 Shade excluded regions (light red)
    for k = 1:size(excludeRegions, 1)
        t1 = timeVec(excludeRegions(k, 1));
        t2 = timeVec(excludeRegions(k, 2));
        fill([t1 t2 t2 t1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            [1 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4, 'DisplayName', 'Excluded');
    end

    % 🟩 Shade extracted trials (light green)
    for p = 1:numTrials
        ix = trialStarts(p);
        t1 = timeVec(ix);
        t2 = timeVec(ix + trialLength - 1);
        fill([t1 t2 t2 t1], [yLimits(1) yLimits(1) yLimits(2) yLimits(2)], ...
            [0.7 1 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'Trial');
    end

    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Raw signal with excluded regions (red) and selected trials (green)');
    legend({'Raw signal'}, 'Location', 'best');
    xlim([timeVec(1), timeVec(end)]);

    % --- Bottom subplot: latent signal (cleaned)
    subplot(2, 1, 2);
    plot(latentTime, latentSignal, 'k');
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Latent signal (excluded regions removed and joined)');
    xlim([latentTime(1), latentTime(end)]);
end
