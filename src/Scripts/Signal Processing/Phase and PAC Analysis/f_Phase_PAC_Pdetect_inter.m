function f_Phase_PAC_Pdetect_inter(P1, P2, A1, A2, v_Data1, v_Data2, srate, srt)

    srate = srate / srt;
    v_Data1 = downsample(v_Data1, srt);
    v_Data2 = downsample(v_Data2, srt);
    A = 0;
    B = 0;

    for c = 1:size(v_Data1, 2)
        A = horzcat(A, v_Data1(:, c)');
        B = horzcat(B, v_Data2(:, c)');
    end

    v_Data1 = A;
    v_Data2 = B;

    v_Time = (0:length(v_Data1) - 1) / srate;
    v_Time = v_Time(:);

    f_ini = A1;
    f_fin = A2;

    s_Liml = P1;
    s_Limh = P2;

    tflim = [];

    nbins = 100;
    l = 2 * pi / nbins;
    index = -pi:l:pi;
    %% filtros
    %Phase
    s_Liml = P1;
    s_Limh = P2;
    st_Filt1 = f_GetIIRFilter(srate, [s_Liml s_Limh]);
    v_Data_Filt = f_IIRBiFilter(v_Data1, st_Filt1);

    % figure
    % ax(1) = subplot(2,1,1);
    %     plot(v_Time,v_Data);
    %     ylabel('Original');
    % ax(2) = subplot(2,1,2);
    %     plot(v_Time,v_Data_Filt);
    %     ylabel('Filtered');
    %     linkaxes(ax,'x');
    %% detec

    %Theta
    s_minT = 0.5 / s_Limh;
    s_maxT = 0.5 / s_Liml;

    %peak detection

    y = v_Data_Filt;
    x = (0:size(y, 1) - 1) / srate;
    % peaks

    [v_pk, v_ppk] = findpeaks(y); %,'MinPeakHeight', Min_volt); %,'MinPeakDistance', s_minT);
    v_tlk_t = v_ppk / srate;
    %figure;hist(v_pk, 20);
    %Max_volt = mean(v_pk);
    Max_volt = mean(v_pk) - std(v_pk);
    % figure
    % plot(x,y,v_tlk_t,v_pk,'ro')

    %lows
    [v_lw, v_plw] = findpeaks(-y);
    v_lw = -v_lw;
    v_tlw_t = v_plw / srate;

    % figure
    % plot(x,y,v_tlw_t,v_lw,'ro')

    % figure;hist(v_lw, 20);
    %Min_volt = mean(v_lw)+ 0.5*(std(v_lw));
    %Min_volt = mean(v_lw);
    Min_volt = mean(v_lw) + (std(v_lw));

    if length(v_lw) > length(v_pk)
        v_lw = v_lw(1:length(v_pk));
    elseif length(v_lw) < length(v_pk)
        v_pk = v_pk(1:length(v_lw));
    end

    %altura minima pp
    %  figure;hist(v_pk-v_lw, 20);
    Min_Vpp = mean(v_pk - v_lw) %- 1* std(v_lw-v_pk);

    [v_lw, v_plw] = findpeaks(-y); %,'MinPeakHeight', Min_volt); %,'MinPeakDistance', s_minT);
    v_lw = -v_lw;
    v_tlw_t = v_plw / srate;

    %zeros
    % Signal
    zci = @(v) find(v(:) .* circshift(v(:), [-1 0]) <= 0); % Returns Zero-Crossing Indices Of Argument Vector
    v_zs = zci(y); % Approximate Zero-Crossing Indices

    % figure
    % plot(x,y,v_tlk_t,v_pk,'ro')
    % hold on
    % plot(x,y,v_tlw_t,v_lw,'go')
    % hold on
    % plot(x(v_zs), y(v_zs),'ko')

    %% vect

    count_lw = 2;
    count_zs = 1;
    count_pk = 1;

    v_sel = ones(size(v_plw));
    v_sel(1) = 0;

    for p = 1:size(v_plw)

        if v_plw(count_lw) < v_zs(count_zs)
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        else
            break
        end

    end

    while count_lw < size(v_plw, 1) - 1
        idx_lw_actual = v_plw(count_lw);
        idx_lw_ant = v_plw(count_lw - 1);
        idx_lw_post = v_plw(count_lw + 1);

        a = find(v_zs < v_plw(count_lw));
        idx_zs_pre = v_zs(a(end));

        b = find(v_zs > v_plw(count_lw));
        idx_zs_post = v_zs(b(1));

        c = find(v_ppk < v_plw(count_lw));
        idx_pk_pre = v_ppk(c(end));

        d = find(v_ppk > v_plw(count_lw));
        idx_pk_post = v_ppk(d(1));

        if idx_lw_post < idx_zs_post
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif idx_lw_post < idx_pk_post
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif idx_pk_post < idx_zs_post
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif idx_lw_ant > idx_zs_pre
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif idx_lw_ant > idx_pk_pre
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif idx_pk_pre > idx_zs_pre
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif v_Data_Filt(idx_lw_actual) > Min_volt %|| abs(v_Data_Filt(idx_lw_actual)-v_Data_Filt(idx_pk_pre))> abs(Min_Vpp)
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif v_Data_Filt(idx_pk_pre) < Max_volt
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif v_Data_Filt(idx_pk_post) < Max_volt
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        elseif v_Time(idx_lw_actual) - v_Time(idx_pk_pre) < s_minT || ...
                v_Time(idx_lw_actual) - v_Time(idx_pk_pre) > s_maxT || ...
                v_Time(idx_pk_post) - v_Time(idx_lw_actual) < s_minT || ...
                v_Time(idx_pk_post) - v_Time(idx_lw_actual) > s_maxT
            v_sel(count_lw) = 0;
            count_lw = count_lw + 1;
        else
            v_sel(count_lw) = 1;
            count_lw = count_lw + 1;

        end

    end

    v_sel(end - 1) = 0;
    v_sel(end) = 0;

    v_sel_lws = v_sel .* v_plw;

    v_sel_lws(v_sel_lws == 0) = [];
    v_tlw_sel = v_sel_lws / srate;
    v_lw_sel = v_Data_Filt(v_sel_lws);

    s_tps_w = floor(2 * s_maxT * srate);

    m_lows_full = zeros(2 * s_tps_w + 1, size(v_sel_lws, 1) - 2);
    m_lows_theta = zeros(2 * s_tps_w + 1, size(v_sel_lws, 1) - 2);

    for count = 3:size(v_sel_lws, 1)
        m_lows_full(:, count - 2) = v_Data2((v_sel_lws(count) - s_tps_w):(v_sel_lws(count) + s_tps_w));
        m_lows_theta(:, count - 2) = v_Data_Filt((v_sel_lws(count) - s_tps_w):(v_sel_lws(count) + s_tps_w));
    end

    av = mean(m_lows_theta, 2);
    av_full = mean(m_lows_full, 2);

    %% tf

    ps_MinFreqHz = A1;
    ps_MaxFreqHz = A2;
    ps_FreqSeg = 2 * (ps_MaxFreqHz - ps_MinFreqHz);
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

    %% tranforme
    [ga, ~, fr] = f_GaborAWTransformMatlab(v_Data2((v_sel_lws(count) - s_tps_w):(v_sel_lws(count) + s_tps_w)), srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles, ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases, ps_TimeStep);

    [p_av, l_av] = findpeaks(-av);

    if length(l_av) == 3
        ct1 = l_av(1) - 1;
        ct2 = l_av(3) + 1;
    elseif length(l_av) == 1
        ct1 = 1;
        ct2 = length(av);
    elseif length(l_av) == 2
        ct1 = 1;
        ct2 = length(av);
    elseif length(l_av) > 3
        ct1 = l_av(2);
        ct2 = l_av(4);
    end

    av = av(ct1:ct2);
    av_full = av_full(ct1:ct2);
    m_lows_full = m_lows_full(ct1:ct2, :);
    m_lows_theta = m_lows_theta(ct1:ct2, :);

    ga = ga(:, ct1:ct2);

    ga_av = zeros(size(ga, 1), size(ga, 2));

    v_phase_num = zeros(1, size(index, 2) - 1);
    v_num_temp = zeros(1, size(index, 2) - 1);

    m_phase = zeros(size(ga, 1), size(index, 2) - 1);
    m_phase_temp = zeros(size(ga, 1), size(index, 2) - 1);

    s_tps_w = floor((size(m_lows_theta, 1) - 1) / 2);

    [p2_av, l2_av] = findpeaks(av);
    index_pi = index(1:nbins);
    %figure
    for count = 1:size(m_lows_theta, 2)
        [ga, ti, fr] = f_GaborAWTransformMatlab(m_lows_full(:, count), srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles, ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases, ps_TimeStep);

        % f_ImageMatrix(ga,ti,fr, [], 'jet', 256);
        if size(ga, 2) > size(ga_av, 2)
            ga_temp = ga(:, 1:size(ga_av, 2));
        elseif size(ga, 2) < size(ga_av, 2)
            ga_temp = zeros(size(ga_av));
            ga_temp(:, 1:size(ga, 2)) = ga;

        else
            ga_temp = ga;
        end

        ga_av = ga_av + ga_temp;

        ga = ga(:, l2_av(1):l2_av(2));
        v_Data_phase = angle(hilbert(m_lows_theta(l2_av(1):l2_av(2), count)));

        v_num_temp = zeros(1, size(index, 2) - 1);

        for countp = 1:nbins
            N = find(v_Data_phase < index(countp + 1) & v_Data_phase >= index(countp));
            m_phase(:, countp) = m_phase(:, countp) + sum(ga(:, N), 2);
            v_phase_num(countp) = v_phase_num(countp) + size(N, 1);
            v_num_temp(countp) = v_num_temp(countp) + size(N, 1);
        end

        m_phase_temp = m_phase ./ (repmat(v_num_temp, size(m_phase, 1), 1));
        hist_temp = sum(m_phase_temp, 1);
        %              figure;
        %              bar(index_pi,hist_temp);

        hist_temp(isinf(hist_temp) | isnan(hist_temp)) = 0;
        v_angles(count, 1) = mean(hist_temp .* exp(1i .* index_pi));

    end

    save('angulos.mat', 'v_angles');
    m_phase_av = m_phase ./ (repmat(v_phase_num, size(m_phase, 1), 1));
    %histograms
    hist_pre_1 = sum(m_phase_av, 1);

    m_phase_av_mean = mean(m_phase_av, 2);
    m_phase_av_std = std(m_phase_av, [], 2);
    m_phase_av_aux = m_phase_av - repmat(m_phase_av_mean, 1, size(m_phase_av, 2));
    m_phase_av_aux = m_phase_av_aux ./ repmat(m_phase_av_std, 1, size(m_phase_av, 2));
    m_phase_av_aux_pre = m_phase_av_aux;
    %

    l1 = -2;
    l2 = 2;

    %plot log
    %close all
    figure;
    ax(1) = subplot(1, 1, 1); %('position', [0.1 0.31 0.8 0.5]);
    v_FreqAxisTemp = fr;
    s_InvertImage = 0;
    s_NonEquAxis = 1;
    v_FreqAxisTemp = log10(v_FreqAxisTemp);
    f_ImageMatrix(m_phase_av_aux, index_pi, v_FreqAxisTemp, [l1 l2], 'jet', 256, 0, s_InvertImage, 1, s_NonEquAxis);
    v_YTick = [1.5, 2, 2.5];
    set(ax(1), 'YTick', v_YTick, 'YTickLabel', round(10 .^ v_YTick));
    ylabel('Frequency (Hz) - Log');
    colorbar

    save('m_PAC.mat', 'm_phase_av_aux');

    for c = 1:950
        b = m_phase_av_aux(c, :);
        v_mod(c) = max(b);
    end

    figure
    plot(fr, v_mod);

    v_mod2 = sum(m_phase_av_aux, 1);
    figure
    plot(index_pi, v_mod2);

    ga_av = ga_av ./ size(m_lows_theta, 2);

    ga_av = ga_av(:, l2_av(1):l2_av(2));
    av = av(l2_av(1):l2_av(2));

    t_tf = (0:size(ga_av, 2) - 1) / srate;

    [tf_av_aux] = f_Mat_to_zscore(ga_av);

    tf_lim1 = [];
    %
    figure;
    ax(1) = subplot(2, 1, 1); %('position', [0.1 0.31 0.8 0.5]);
    f_ImageMatrix(tf_av_aux, t_tf, fr, tf_lim1, 'jet', 256);
    ax(2) = subplot(2, 1, 2); %('position', [0.1 0.1 0.8 0.5]);
    plot(t_tf, av);
    %linkaxes(ax,'x');

    minh1 = min(hist_pre_1);
    maxh1 = max(hist_pre_1);

    dis1 = maxh1 - minh1;

    v_pha_av_pre1 = mean(hist_pre_1 .* exp(1i .* index_pi));

    figure('name', ' phase modulation average');
    max_lim = max(max(abs(real(v_pha_av_pre1))), max(abs(imag(v_pha_av_pre1))));
    %max_lim = 160; %Th cx
    %max_lim = 60; % cx Th
    x_fake = [0 max_lim 0 -max_lim];
    y_fake = [max_lim 0 -max_lim 0];
    h_fake = compass(x_fake, y_fake);
    hold on;
    compass(real(v_pha_av_pre1), imag(v_pha_av_pre1), 'b')
    set(h_fake, 'Visible', 'off');
    hold off;
    save('v_phase.mat', 'v_pha_av_pre1');

    i_mod = imag(v_pha_av_pre1);
    r_mod = real(v_pha_av_pre1);

    [angulo, amplitud] = cart2pol(real(v_pha_av_pre1), imag(v_pha_av_pre1));

    angulo = angle(v_pha_av_pre1)
    amplitud = abs(v_pha_av_pre1)
    num_cicles = size(m_lows_theta, 2)
    %close all
    Vpp = max(av) - min(av);

end
