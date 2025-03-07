function f_RMS_zscore(m_Trials_selected, srate, v_Timew, GF, h, Ch)

    m_TrialsRMS = zeros(size(m_Trials_selected));

    for cr = 1:size(m_Trials_selected, 2)
        m_TrialsRMS(:, cr) = rms(m_Trials_selected(:, cr), 2);
    end

    v_BL_mean = mean(mean(m_TrialsRMS(1:round(.15 * srate), :), 2));

    if GF == 0

        if size(m_Trials_selected, 2) == 1
            v_std_mean = std(m_TrialsRMS(1:round(.15 * srate), :));
            m_TrialsZ = (m_TrialsRMS - v_BL_mean) ./ v_std_mean;
            figure
            plot(v_Timew, m_TrialsZ);
            ylabel('RMS zscore');
            xlabel('Time(s)');
            title(h.recChNames(Ch));

        else

            v_std_mean = std(std(m_TrialsRMS(1:round(.15 * srate), :)));
            m_TrialsZ = (m_TrialsRMS - v_BL_mean) ./ v_std_mean;
            v_t = 1:1:size(m_Trials_selected, 2);
            [m_Time, m_trial] = meshgrid(v_Timew, v_t);

            figure;
            surf(m_Time, m_trial, m_TrialsZ');

            colormap(jet);
            shading interp;
            view([4.4200 71.6000]);
            zlabel('RMS zscore');
            xlabel('Time(s)');
            ylabel('Trials');
            title(h.recChNames(Ch));
            colorbar
            zlim([-5, 60]);

        end

    else
        v_std_mean = std(m_TrialsRMS(1:round(.15 * srate), :));
        m_TrialsZ = (m_TrialsRMS - v_BL_mean) ./ v_std_mean;
        figure
        plot(v_Timew, m_TrialsZ);
        ylabel('RMS zscore');
        xlabel('Time(s)');
        title(h.recChNames(Ch));

    end

end
