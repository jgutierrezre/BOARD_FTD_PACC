function [m_Trials,v_ini] = f_Cols_Trials(m_Data,Ch,Th,Stim_num,ti,num_Ch,srate)
%generates a matrix of trials in columns, started by each stimulus)
% m_Trials trial matrix in columns
% m_Data data matrix
% Ch channel for the analisis
% Th threshold for the stimuly
% ti=initial time

Th=abs(Th);

    if Stim_num == 0

       L= length(m_Data);
        
        c=ceil(srate*(-ti));
        
        v_random=randi([c+1 (size(m_Data,1)-(L+1))],[1 60]);
        v_ini(1:c)=m_Data(v_random(1):c+v_random(1)-1,Ch);
        for p=1:60
            m_Trials(:,p) = m_Data(v_random(p):v_random(p)+L-1,Ch);    
        end
        
   
        
    else
        count=0;
        p = 1;
        q = 1;
        ci=1;
        for i=1:size(m_Data, 1)
            if abs(m_Data(i,num_Ch))<Th && count==0
                 v_ini(ci)=m_Data(i,Ch);
                ci=ci+1;
            elseif abs(m_Data(i,num_Ch))>Th && count==0
                m_Trials(p,q) = m_Data(i,Ch);
                p=p+1;
                count=1;
            elseif abs( m_Data(i,num_Ch))>Th&& count==1
                m_Trials(p,q) = m_Data(i,Ch);
                p=p+1;
                count=1;
            elseif abs(m_Data(i,num_Ch))<Th && count==1
                m_Trials(p,q) = m_Data(i,Ch);
                p=p+1;
                count=2;
            elseif abs(m_Data(i,num_Ch))<Th && count==2
                m_Trials(p,q) = m_Data(i,Ch);
                p=p+1;
            elseif abs(m_Data(i,num_Ch))>Th && count==2
                count=1;
                q=q+1;
                p=1;
                m_Trials(p,q) = m_Data(i,Ch);
            end
        end    
    end
     v_Time_Temp = (0:size(m_Trials,1)-1)./srate;



end

