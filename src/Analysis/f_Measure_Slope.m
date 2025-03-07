function [S1, S2] = f_Measure_Slope(m_TrialsAVw, r1t1, r1t2, r2t1, r2t2, srate, v_Timew, Stim_num)

    if Stim_num == 2

        a = find(v_Timew >= r1t1);
        s_11 = a(1) - 1;
        b = find(v_Timew >= r1t2);
        s_12 = b(1) - 1;
        c = find(v_Timew >= r2t1);
        s_21 = c(1) - 1;
        d = find(v_Timew >= r2t2);
        s_22 = d(1) - 1;

        m_r1 = m_TrialsAVw(s_11:s_12, :);
        m_r2 = m_TrialsAVw(s_21:s_22, :);

        for cd = 1:size(m_r1, 2)
            m_r1diff (:, cd) = diff(m_r1(:, cd));
            m_r2diff (:, cd) = diff(m_r2(:, cd));
        end

        S1 = max(abs(m_r1diff)) * srate;
        S2 = max(abs(m_r2diff)) * srate;

        S1'
        S2'
    elseif Stim_num == 1

        a = find(v_Timew >= r1t1);
        s_11 = a(1) - 1;
        b = find(v_Timew >= r1t2);
        s_12 = b(1) - 1;

        m_r1 = m_TrialsAVw(s_11:s_12, :);

        for cd = 1:size(m_r1, 2)
            m_r1diff (:, cd) = diff(m_r1(:, cd));
        end

        S1 = max(abs(m_r1diff)) * srate;
        S1'
        S2 = 0;
    else
        S1 = 0;
        S2 = 0;
    end

end
