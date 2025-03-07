function [m_Trials_w, m_Trials_wt, v_Timew, m_Trials_freq] = f_Load_GF(m_Data, v_Time, ti, tf, srate, h, Ch)
    %UNTITLED3 Summary of this function goes here

    if tf == inf
        tf = v_Time(end);
    end

    if tf >= v_Time(end)
        tf = v_Time(end);
    end

    if ti <= 0
        t1 = 1;
        ti = 0;
    else
        t1 = floor(ti * srate);
    end

    t2 = floor(tf * srate);

    if t2 <= t1
        t2 = t1 + 1;
    end

    m_Trials = m_Data(:, Ch);
    m_Trials_selected = m_Trials(t1:t2);
    m_Trials_w = m_Trials(t1:t2);
    m_Trials_wt = m_Trials(t1:t2);
    m_Trials_freq = m_Trials(t1:t2);
    v_Timew = linspace(ti, tf, size(m_Trials_selected, 1))';

    figure
    plot(v_Timew, m_Trials_selected);
    ylabel(h.recChUnits(Ch));
    xlabel('Time(s)');
    title(h.recChNames(Ch));
end
