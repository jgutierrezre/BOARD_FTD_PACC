function  f_PSD(v_Data,srt,srate,Tr,h,Ch)


srate =srate/srt;
v_Data= downsample(v_Data,srt);

[pxx,f] = pwelch(v_Data,[],[],[],srate);

   
    figure;
    plot(f,pxx)
    xlabel('Frequency (Hz)')
    ylabel('PSD');
    title(h.recChNames(Ch));
    
     figure;
    plot(f,pxx(:,Tr))
    xlabel('Frequency (Hz)')
    ylabel('PSD');
    title(['Trial ',num2str(Tr)]);
    

end

