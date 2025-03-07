function f_Filter_order_manual(F1, F2, s_Order, v_sig, srate, vtime, m_Sig, v_Timew, Tr, h, Ch)

    v_Fc = [F1, F2];
    Filtro1 = fdesign.bandpass('n,f3db1,f3db2', s_Order, v_Fc(1), v_Fc(2), srate);
    h_Filt = design(Filtro1, 'butter');

    fvtool(h_Filt)

    v_h_Filt = filter(h_Filt, v_sig);
    v_h_Filt = filter(h_Filt, fliplr(v_h_Filt));
    v_h_Filt = fliplr(v_h_Filt);

    % filter matrix

    figure;
    ax(1) = subplot(2, 1, 1);
    plot(vtime, v_sig);
    ylabel(h.recChUnits(Ch));
    title('Original');
    ax(2) = subplot(2, 1, 2);
    plot(vtime, v_h_Filt);
    ylabel(h.recChUnits(Ch));
    title('Filtered');
    xlabel('Time(s)');
    linkaxes(ax, 'x');

    m_Sig_Filt = filter(h_Filt, m_Sig); % filtrado
    m_Sig_Filt = filter(h_Filt, fliplr(m_Sig_Filt)); % se invierte la se�al y se filtra nuevamente para corregir el corrimiento de fase
    m_Sig_Filt = fliplr(m_Sig_Filt); % se invierte para tener la se�al filtrada en el orden correto sin corrimiento

    %plor filtered trial

    m_ph = angle(hilbert(m_Sig_Filt(:, Tr)));

    figure;
    ax(1) = subplot(3, 1, 1);
    plot(v_Timew, m_Sig(:, Tr));
    ylabel(h.recChUnits(Ch));
    title(['Original Trial ', num2str(Tr)]);
    ax(2) = subplot(3, 1, 2);
    plot(v_Timew, m_Sig_Filt(:, Tr));
    ylabel(h.recChUnits(Ch));
    title(['Filter Trial ', num2str(Tr)]);
    xlabel('Time(s)');
    ax(3) = subplot(3, 1, 3);
    plot(v_Timew, m_ph);
    ylabel(h.recChUnits(Ch));
    title(['Filter Trial ', num2str(Tr)]);
    xlabel('Time(s)');
    linkaxes(ax, 'x');

    %plot filtered
    resp = linspace(1, size(m_Sig, 2), size(m_Sig, 2));
    figure;
    f_ImageMatrix(m_Sig_Filt', v_Timew, resp, [], 'jet', 256);
    xlabel('Time(s)');
    ylabel('Trials');
    title('Filtered');

    % envelope

    m_H = hilbert(m_Sig_Filt);
    m_envelope = abs(m_H);
    figure;
    f_ImageMatrix(m_envelope', v_Timew, resp, [], 'jet', 256);
    xlabel('Time(s)');
    ylabel('Trials');
    title('Envelope');

    %phase
    m_phase = angle(m_H);
    figure;
    f_ImageMatrix(m_phase', v_Timew, resp, [], 'jet', 256);
    xlabel('Time(s)');
    ylabel('Trials');
    title('Phase');
end
