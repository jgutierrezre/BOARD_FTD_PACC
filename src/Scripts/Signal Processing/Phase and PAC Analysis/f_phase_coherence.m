function f_phase_coherence( m_Ch1, m_Ch2,srate, s_Step, ps_MinFreqHz, ps_MaxFreqHz, wd ,v_Timew)
%F_COHERENCY Summary of this function goes here
% m_Ch1 signal 1 
% m_Ch2 signal 2
% srate 
% s_Step resample rate 
% ps_MinFreqHz min frequnecy
% ps_MaxFreqHz max frequency
% wd  tamaño de ventana 0.05; % depende de las frecuencias que me importan (100ms para 100Hz, 500Hz 20ms diferentes ventas para diferentes frecuencias 



srate=srate/s_Step;
m_Ch1 = downsample(m_Ch1,s_Step);
m_Ch2 = downsample(m_Ch2,s_Step);
v_Time = downsample(v_Timew,s_Step);

%% TF

ps_FreqSeg = round((ps_MaxFreqHz-ps_MinFreqHz));
ps_StDevCycles = 1; 
ps_Magnitudes = 0;
ps_SquaredMag = 0; 
ps_MakeBandAve = 0; 
ps_Phases = 1;
ps_TimeStep = []; 
        

SxyGlb=0;
RGlb=0;
p=1;


Coh_Glb = zeros(ps_FreqSeg,floor(size(v_Time,1)/floor(wd*srate))); 
count=1;
for ct=1:floor(wd*srate):(size(m_Ch1,1)-wd*srate)
    for trail=1:size(m_Ch1,2)
        [waveCh2, v_TimeAxis, v_FreqAxis] = f_GaborAWTransformMatlab(m_Ch2(ct:floor(ct+floor(wd*srate)),trail), srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles,ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases,ps_TimeStep);
        [waveCh1, v_TimeAxis, v_FreqAxis] = f_GaborAWTransformMatlab(m_Ch1(ct:floor(ct+floor(wd*srate)),trail), srate, ps_MinFreqHz, ps_MaxFreqHz, ps_FreqSeg, ps_StDevCycles,ps_Magnitudes, ps_SquaredMag, ps_MakeBandAve, ps_Phases,ps_TimeStep); 
        m_resta = waveCh1-waveCh2;
        v_Coh = sum(exp(1i.*m_resta),2);
        Coh_Glb(:,count) =Coh_Glb(:,count)+v_Coh;
    end
    
    Coh_Glb(:,count) =abs(Coh_Glb(:,count)/(size(v_TimeAxis,2)*size(m_Ch1,2))); %mirar la fase
    count=count+1;
end
 t = linspace(v_Time(1),v_Time(end),size(Coh_Glb,2));
%     figure;
%         f_ImageMatrix(Coh_Glb,t,v_FreqAxis,[0 0.5], 'jet', 256);
%        
 figure;
    v_FreqAxisTemp=v_FreqAxis;
    s_InvertImage = 0;
    s_NonEquAxis = 1;
    v_FreqAxisTemp = log10(v_FreqAxisTemp);
    ax(1)=subplot(1,1,1);
   f_ImageMatrix(Coh_Glb,t,v_FreqAxisTemp,[0 0.5], 'jet', 256,0, s_InvertImage, 1, s_NonEquAxis);
    v_YTick = get(ax, 'YTick');
    set(ax(1), 'YTick', v_YTick, 'YTickLabel', round(10.^v_YTick));
    ylabel('Frequency (Hz) - Log');
    xlabel('Time(s)');
    title('Phase Coherence');
    colorbar
end


