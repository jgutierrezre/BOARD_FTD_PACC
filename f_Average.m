function    f_Average(m_Trials_selected,v_Timew,h,Ch)

av_Data = mean(m_Trials_selected,2);
figure
plot(v_Timew,av_Data);
    ylabel(h.recChUnits(Ch));

    xlabel('Time(s)');
    title(h.recChNames(Ch)); 

end

