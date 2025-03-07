function [m_TrialsAVw] = f_measure_plot(m_Trials_wt, num_stim, v_Timew, option, GF, x, h, Ch)

    if GF == 1
        num_stim = 0;
    end

    if num_stim == 0
        figure
        m_TrialsAVw = m_Trials_wt;
        plot(v_Timew, m_Trials_wt);
        ylabel(h.recChUnits(Ch));
        xlabel('Time(s)');
        title(h.recChNames(Ch));
    elseif option == 1
        v_Average = mean(m_Trials_wt, 2);
        ct = 1;
        cti = 1;

        while ct < size(m_Trials_wt, 2) - (x - 1)
            A = m_Trials_wt(:, ct:ct + (x - 1));
            m_TrialsAVw(:, cti) = mean(A, 2);
            ct = ct + x;
            cti = cti + 1;
        end

        figure
        plot(v_Timew, m_TrialsAVw)
        ylabel(h.recChUnits(Ch));
        xlabel('Time(s)');
        title(h.recChNames(Ch));
        hold on
        plot(v_Timew, v_Average, 'k', 'LineWidth', 2);
        hold off
    else
        v_Average = mean(m_Trials_wt, 2);
        m_TrialsAVw = v_Average;
        figure
        plot(v_Timew, v_Average);
        ylabel(h.recChUnits(Ch));
        xlabel('Time(s)');
        title(h.recChNames(Ch));
    end

end
