function f_view_art(m_Trials, h, Ch)
    a = mean(m_Trials, 2);

    figure;
    plot(a);
    ylabel(h.recChUnits(Ch));
    xlabel('Position');
    title(h.recChNames(Ch));

end
