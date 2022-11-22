function [m_Trials_selected,m_Trials_freq_sel] = f_Resample_mat(m_Trials_selected,m_Trials_freq_sel,v_Timew,f_target,srate)

srt =round(srate/f_target);

srate = srate/srt;
m_Trials_selected = downsample(m_Trials_selected,srt);
m_Trials_freq_sel = downsample(m_Trials_freq_sel,srt);
v_Timew = downsample(v_Timew,srt);
end

