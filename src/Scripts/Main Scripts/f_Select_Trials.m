function [m_Trials_selected, m_Trials_freq_sel] = f_Select_Trials(m_Trials_wt, v_trials, v_Timew, all_trials, m_Trials_freq, GF, h, Ch)

    if GF == 1
        m_Trials_freq_sel = m_Trials_freq;
        m_Trials_selected = m_Trials_wt;
    else

        if all_trials == 1
            m_Trials_selected = m_Trials_wt;
            m_Trials_freq_sel = m_Trials_freq;
        else
            v_trials = v_trials(2:end - 1);
            m_Trials_selected = m_Trials_wt(:, v_trials);
            m_Trials_freq_sel = m_Trials_freq(:, v_trials);
        end

    end

    figure
    plot(v_Timew, m_Trials_selected(:, 1));
    ylabel('Amplitude');
    title(h.recChNames(Ch));

end
