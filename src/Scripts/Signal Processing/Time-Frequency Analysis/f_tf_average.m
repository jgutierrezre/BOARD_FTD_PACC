function f_tf_average(v_Data, ps_MinFreqHz, ps_MaxFreqHz, srt, srate, v_Timew, ps_StDevCycles, norm, log, h, Ch)

    srate = srate / srt;
    v_Data = downsample(v_Data, srt);
    v_Timew = downsample(v_Timew, srt);

    %%wavellet

    %%Wavellet parameter setup
    ps_FreqSeg = 2 * (ps_MaxFreqHz - ps_MinFreqHz);
    ps_Magnitudes = 1;
    ps_SquaredMag = 0;
    ps_MakeBandAve = 0;
    ps_Phases = 0;
    ps_TimeStep = [];

    %frequency windows
    fstep = (ps_FreqSeg - 1) / ps_MaxFreqHz;

    %transforme

    [tf_Data, ti, fr] = f_MorseAWTransformMatlab(v_Data(:, 1), srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles, ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases, ps_TimeStep);
    tf_av = zeros(size(tf_Data));

    for i = 1:size(v_Data, 2)
        [tf_Data, ti, fr] = f_MorseAWTransformMatlab(v_Data(:, i), srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles, ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases, ps_TimeStep);
        tf_av = tf_av + tf_Data;
    end

    tf_2 = tf_av ./ size(v_Data, 2);
    %normalization
    if norm == 1
        [tf_fin] = f_Mat_to_zscore(tf_2);
    else
        tf_fin = tf_2;
    end

    if log == 1
        figure;
        v_FreqAxisTemp = fr;
        s_InvertImage = 0;
        s_NonEquAxis = 1;
        v_FreqAxisTemp = log10(v_FreqAxisTemp);
        ax(1) = subplot(1, 1, 1);
        f_ImageMatrix(tf_fin, v_Timew, v_FreqAxisTemp, [], 'jet', 256, 0, s_InvertImage, 1, s_NonEquAxis);
        v_YTick = get(ax, 'YTick');
        set(ax(1), 'YTick', v_YTick, 'YTickLabel', round(10 .^ v_YTick));
        ylabel('Frequency (Hz) - Log');
        xlabel('Time(s)');
        title(h.recChNames(Ch));
        colorbar
    else
        figure
        f_ImageMatrix(tf_fin, v_Timew, fr, [], 'jet', 256);
        ylabel('Frequency (Hz)');
        xlabel('Time(s)');
        title(h.recChNames(Ch));

    end

end
