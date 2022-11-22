function [MInorm fm fp]=f_PAC_sing( v_Data, srate, p_ini, p_fin, pstep, a_ini, a_fin, astep, srt,Tr,GF)

%  pstep phase step
% v_Data data colum vector
% p_ini initial phase frequency
% p_fin final phase frequency
% astep amplitude step
% a_ini inital amplitude frequency
% a_fin final amplitude frequency
% srt = resample rate

if GF==1
    Tr=1;
end


%resample
srate =srate/srt;
v_Data= downsample(v_Data(:,Tr),srt);


data_length = size(v_Data,1);

%% initial filter
% Phase
    phase_range = [p_ini:pstep:p_fin];
    m_pha = zeros(size(v_Data,1),size(phase_range,2));
    for countp = 1:size(phase_range,2) 
        st_Filt1 = f_GetIIRFilter(srate, [phase_range(countp) phase_range(countp)+1]);
       
        PhaseFil = f_IIRBiFilter(v_Data, st_Filt1);
       
        m_phaFil(:,countp) = PhaseFil;
        m_pha(:,countp) = angle(hilbert(PhaseFil));
    end

% Amplitude
    amp_range = [a_ini:astep:a_fin];
   m_amp = zeros(size(v_Data,1),size(amp_range,2));
   %for counta = 1:(size(amp_range,2) 
        for counta = 1:178
        st_Filt1 = f_GetIIRFilter(srate, [amp_range(counta) amp_range(counta)+10]);
        AmpFil = f_IIRBiFilter(v_Data, st_Filt1);
       
        m_ampFil(:,counta) = AmpFil;
       m_amp(:,counta) = abs(hilbert(AmpFil));
       %counta
    end
    

    
            
%% Phase & amplitude vectors


for countp = 1:size(phase_range,2)   
    for counta = 1:size(amp_range,2)
        
       
            a =m_amp(:,counta).*exp(1i.*m_pha(:,countp));
            b = mean(a);  
            
            MI(counta,countp)= abs(b);
            d = sum(a);
            e = sum(m_amp(:,counta));
            f = d/e;
            MInorm(counta,countp)=abs(f);  
            
    end
end
fp = phase_range(1):pstep:phase_range(end);
fm = amp_range(1):astep:amp_range(end);


% figure
% f_ImageMatrix(MI,fp,fm, [], 'jet', 256);
%     ylabel('Freq Amplitude');
%     xlabel('Freq Phase');
%     title(['MI for trial ',num2str(Tr)]);
 figure
f_ImageMatrix(MInorm,fp,fm, [], 'jet', 256);
    ylabel('Freq Amplitude');
    xlabel('Freq Phase');
    title(['MI for trial ',num2str(Tr)]);
    %saveas(gcf,'CGR.jpg');
  
%save('MI.mat','MInorm');

end

