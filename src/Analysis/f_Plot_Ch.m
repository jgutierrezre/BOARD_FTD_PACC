function f_Plot_Ch(v_Time, m_Data, h, num_Ch)
    %generates a plot of the raw signal, according to the channel number
    %v_Time time vector
    % Data matrix in columns
    %structure of the header
    %number of channels
    figure

    for c = 1:num_Ch
        ax(c) = subplot(num_Ch, 1, c);
        plot(v_Time, m_Data(:, c));
        ylabel(h.recChUnits(c));
        xlabel('Time(s)');
        title(h.recChNames(c));
        linkaxes(ax, 'x');
    end

end
