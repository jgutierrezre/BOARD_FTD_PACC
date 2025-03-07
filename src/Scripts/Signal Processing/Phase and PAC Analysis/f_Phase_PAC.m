function f_Phase_PAC(P1, P2, A1, A2, v_Data, srate, srt, Tr, GF, vtime)

    if GF == 1
        Tr = 1;
    end

    srate = srate / srt;
    v_Data = downsample(v_Data(:, Tr), srt);
    vtime = downsample(vtime, srt);

    %lim del tf
    f_ini = A1;
    f_fin = A2;

    %band mod delta
    s_Liml = P1;
    s_Limh = P2;

    tflim = [];

    %% filt

    st_Filt1 = f_GetIIRFilter(srate, [s_Liml s_Limh]);
    v_Data_Filt = f_IIRBiFilter(v_Data, st_Filt1);

    %% tf
    ps_MinFreqHz = f_ini;
    ps_MaxFreqHz = f_fin;
    ps_FreqSeg = 4 * (ps_MaxFreqHz - ps_MinFreqHz);
    ps_StDevCycles = 1;
    ps_Magnitudes = 1;
    ps_SquaredMag = 0;
    ps_MakeBandAve = 0;
    ps_Phases = 0;
    ps_TimeStep = [];
    %frequency windows
    fstep = (ps_FreqSeg - 1) / ps_MaxFreqHz;
    %time windows
    tstep = (ps_FreqSeg - 1) / ps_MaxFreqHz;

    %% tranf mod
    [m_tf_Data, ti, fr] = f_MorseAWTransformMatlab(v_Data, srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles, ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases, ps_TimeStep);
    %phase matrix
    nbins = 100;
    l = 2 * pi / nbins;
    index = -pi:l:pi;

    m_phase = zeros(size(m_tf_Data, 1), size(index, 2) - 1);
    v_phase_num = zeros(1, size(index, 2) - 1);
    v_Data_phase = angle(hilbert(v_Data_Filt));

    %fill the matrix
    for countp = 1:nbins
        N = find(v_Data_phase < index(countp + 1) & v_Data_phase >= index(countp));
        m_phase(:, countp) = m_phase(:, countp) + sum(m_tf_Data(:, N), 2);
        v_phase_num(countp) = v_phase_num(countp) + size(N, 1);
    end

    %average
    m_phase_av = m_phase ./ (repmat(v_phase_num, size(m_phase, 1), 1));
    %histograms
    hist_pre_1 = sum(m_phase_av, 1);

    %index
    index_pi = index(1:nbins);
    %normalizacion
    m_phase_av_mean = mean(m_phase_av, 2);
    m_phase_av_std = std(m_phase_av, [], 2);
    m_phase_av_aux = m_phase_av - repmat(m_phase_av_mean, 1, size(m_phase_av, 2));
    m_phase_av_aux = m_phase_av_aux ./ repmat(m_phase_av_std, 1, size(m_phase_av, 2));
    m_phase_av_aux_pre = m_phase_av_aux;

    %% graficas
    % figure;
    %         f_ImageMatrix(m_phase_av,index_pi,fr,tflim, 'jet', 256);
    %        ylabel('Frequency (Hz) - Log');
    %         xlabel('Phase (rad)');
    %       title(['Trial ',num2str(Tr)]);

    figure;
    f_ImageMatrix(m_phase_av_aux, index_pi, fr, [], 'jet', 256);
    ylabel('Frequency (Hz) - Log');
    xlabel('Phase (rad)');
    title(['NORM Trial ', num2str(Tr)]);

    minh1 = min(hist_pre_1);
    maxh1 = max(hist_pre_1);

    dis1 = maxh1 - minh1;

    figure('name', 'histograms', 'position', [200, 50, 1000, 700]);
    bar(index_pi, hist_pre_1);
    set(gca, 'YLim', [(minh1 - dis1 * 0.1) (maxh1 + dis1 * 0.1)]);

    v_pha_av_pre1 = mean(hist_pre_1 .* exp(1i .* index_pi));

    figure('name', ' phase modulation average');
    max_lim = max([abs(real(v_pha_av_pre1)) abs(imag(v_pha_av_pre1))]);
    x_fake = [0 max_lim 0 -max_lim];
    y_fake = [max_lim 0 -max_lim 0];
    h_fake = compass(x_fake, y_fake);
    hold on;
    compass(real(v_pha_av_pre1), imag(v_pha_av_pre1), 'b')
    set(h_fake, 'Visible', 'off');

    angulo = imag(v_pha_av_pre1);
    amplitud = real(v_pha_av_pre1);

    v_Fc_P = [P1, P2];
    FiltroP = fdesign.bandpass('n,f3db1,f3db2', 8, v_Fc_P(1), v_Fc_P(2), srate);
    h_FiltP = design(FiltroP, 'butter');

    v_h_Filt_P = filter(h_FiltP, v_Data);
    v_h_Filt_P = filter(h_FiltP, fliplr(v_h_Filt_P));
    v_h_Filt_P = fliplr(v_h_Filt_P);

    v_Fc_A = [A1, A2];
    FiltroA = fdesign.bandpass('n,f3db1,f3db2', 8, v_Fc_A(1), v_Fc_A(2), srate);
    h_FiltA = design(FiltroA, 'butter');

    v_h_Filt_A = filter(h_FiltA, v_Data);
    v_h_Filt_A = filter(h_FiltA, fliplr(v_h_Filt_A));
    v_h_Filt_A = fliplr(v_h_Filt_A);

    phase_sig = angle(hilbert(v_h_Filt_P));

    figure;
    ax(1) = subplot(4, 1, 1);
    plot(vtime, v_Data);
    ylabel('Original');
    %xlim([1 2])
    ax(2) = subplot(4, 1, 2);
    plot(vtime, v_h_Filt_P);
    ylabel('Filtered PHASE');
    %xlim([1 2])
    ax(3) = subplot(4, 1, 3);
    plot(vtime, phase_sig);
    ylabel('Filtered PHASE');
    %xlim([1 2])
    ax(4) = subplot(4, 1, 4);
    plot(vtime, v_h_Filt_A);
    ylabel('Filtered AMPLITUDE');
    xlabel('Time(s)');
    %xlim([1 2]);
    linkaxes(ax, 'x');

end
