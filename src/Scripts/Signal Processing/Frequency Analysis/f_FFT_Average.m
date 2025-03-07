function  f_FFT_Average(v_Data,srate,srt,GF,h,Ch)

srate =srate/srt;
v_Data= downsample(v_Data,srt);


m_Fourier = fft(v_Data,size(v_Data,1));
m_Fourier = m_Fourier(1:floor(size(m_Fourier,1)/2)+1,:);
v_Freq = (srate/2).*linspace(0,1,size(m_Fourier,1));
v_Freq = v_Freq(:);
m_FFT = abs(m_Fourier);

if GF==1
    v_Prom =m_FFT;
else
    v_Prom = mean(m_FFT');
end

figure
plot(v_Freq ,v_Prom); 
xlabel('Frequency (Hz)')
ylabel('FFT Average');
title(h.recChNames(Ch));  

end

