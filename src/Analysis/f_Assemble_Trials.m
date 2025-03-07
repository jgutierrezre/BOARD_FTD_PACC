function [m_Trialsw,t_end,m_Trials_freq] = f_Assemble_Trials(m_Trials,Stim_num,tis,srate)
%asembles a matrix of the selected trials in the determines time window

if Stim_num==0
    m_Trialsw = m_Trials;
    m_Trials_freq = m_Trials;
elseif Stim_num==1
    m_Trialsw = m_Trials;
    m_Trials_freq = m_Trials;
elseif Stim_num==2;    
    ct=1;
    for c=1:2:size(m_Trials,2)-1
        m_Trialsw(:,ct)=vertcat(m_Trials(1:floor(tis*srate),c),m_Trials(:,c+1));
        m_Trials_freq(:,ct) = m_Trials(:,c+1);
        ct=ct+1;
    end
end
  v_Time_Temp = (0:size(m_Trialsw,1)-1)./srate;

     t_end =v_Time_Temp(end);
end

