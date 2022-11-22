function   f_plot_single(m_Data,h,Tr,sel_Ch)

figure;
plot(m_Data(:,Tr))
    ylabel(h.recChUnits(sel_Ch));
    xlabel('Time(s)');
    title(['Trial ',num2str(Tr)]);

end

