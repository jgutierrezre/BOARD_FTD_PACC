function f_Filter_auto(F1, F2, v_Sig, srate, vtime, m_Sig, v_Timew, Tr, h, Ch)

    v_Fc = [F1 F2];

    h_Filt = f_GetIIRFilter(srate, v_Fc);
    fvtool(h_Filt)
    v_h_Filt = f_IIRBiFilter(v_Sig, h_Filt);

    % filter matrix

    figure;
    ax(1) = subplot(2, 1, 1);
    plot(vtime, v_Sig);
    ylabel(h.recChUnits(Ch));
    title('Original');
    ax(2) = subplot(2, 1, 2);
    plot(vtime, v_h_Filt);
    ylabel(h.recChUnits(Ch));
    title('Filtered');
    xlabel('Time(s)');
    linkaxes(ax, 'x');

    m_Sig_Filt = f_IIRBiFilter(m_Sig, h_Filt);

    %plor filtered trial

    figure;
    ax(1) = subplot(2, 1, 1);
    plot(v_Timew, m_Sig(:, Tr));
    ylabel(h.recChUnits(Ch));
    title(['Original Trial ', num2str(Tr)]);
    ax(2) = subplot(2, 1, 2);
    plot(v_Timew, m_Sig_Filt(:, Tr));
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
