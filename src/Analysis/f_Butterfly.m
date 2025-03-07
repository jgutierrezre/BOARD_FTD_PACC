function f_Butterfly(v_Timew, m_Trialswt, h, Ch)

    av_Data = mean(m_Trialswt, 2);
    resp = linspace(1, size(m_Trialswt, 2), size(m_Trialswt, 2));
    l1 = 0.1 * min(min(m_Trialswt));
    l2 = 0.5 * max(max(m_Trialswt));

    figure;
    ax(1) = subplot(2, 1, 1);
    plot(v_Timew, m_Trialswt);
    ylabel('Butterfly');
    title(h.recChNames(Ch));

    ax(2) = subplot(2, 1, 2);
    f_ImageMatrix(m_Trialswt', v_Timew, resp, [l1 l2], 'jet', 256);
    ylabel('Trials');

    linkaxes(ax, 'x');

end
