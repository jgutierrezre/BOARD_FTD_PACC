function [m_Data, srate, h, v_Time, num_Ch] = f_LoadSignal(str_FullFileName)
    %F_LOADSIGNAL generates a matrix with the data in columns
    %   str_FullFileName name of the file
    %   m_Data data in columns for each channel
    %   srate sample rate Hz
    %   h structure of the file header
    %   v_Time time vector
    %   num_Ch number of channels

    if contains(str_FullFileName, 'a') == 1

        [m_Data, si, h] = abfload(str_FullFileName);
        srate = 1 / (si * 10 ^ -6);
        v_Time = (0:size(m_Data, 1) - 1) * si * 10 ^ -6;
        v_Time = v_Time(:);
        num_Ch = size(m_Data, 2);
        Filtro1 = fdesign.bandstop('n,f3db1,f3db2', 100, 59, 61, srate);
        h_Filt = design(Filtro1, 'butter');

        for i = 1:(num_Ch - 1)
            m = mean(m_Data(:, i));
            m_Data(:, i) = m_Data(:, i) - m;
        end

    elseif contains(str_FullFileName, 'e') == 1

        [h, m_Data] = edfread('0001.edf');
        m_Data = m_Data';
        srate = h.frequency(1);
        v_Time = (0:size(m_Data, 1) - 1) / srate;
        v_Time = v_Time(:);
        num_Ch = size(m_Data, 2);
        Filtro1 = fdesign.bandstop('n,f3db1,f3db2', 100, 59, 61, srate);
        h_Filt = design(Filtro1, 'butter');

        for i = 1:(num_Ch - 1)
            m = mean(m_Data(:, i));
            m_Data(:, i) = m_Data(:, i) - m;
        end

        h.recChNames = h.label';
        h.recChUnits = h.units';
        si = (1 / srate) / 10 ^ -6;
    else

    end

end
