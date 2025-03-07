function  f_FFT(v_Data,srate,srt,Tr,h,Ch)

srate =srate/srt;
v_Data= downsample(v_Data,srt);



m_Fourier = fft(v_Data,size(v_Data,1));
m_Fourier = m_Fourier(1:floor(size(m_Fourier,1)/2)+1,:);
v_Freq = (srate/2).*linspace(0,1,size(m_Fourier,1));
v_Freq = v_Freq(:);
m_FFT = abs(m_Fourier);


figure
plot(v_Freq ,m_FFT); 
xlabel('Frequency (Hz)')
ylabel('FFT');
title(h.recChNames(Ch));  

figure
plot(v_Freq ,m_FFT(:,Tr)); 
xlabel('Frequency (Hz)')
ylabel('FFT');
 title(['Trial ',num2str(Tr)]); 
end

