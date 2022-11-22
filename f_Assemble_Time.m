function [m_Trialswt,v_Timew] = f_Assemble_Time(m_Trialsw,v_ini,ti,tf,t_end,srate)

if tf==inf
    tf=t_end;
end
if tf>=t_end;
        tf=t_end;
end

t2=tf;



col1 = zeros(size(m_Trialsw,1),1);

if size(v_ini,2)<=size(m_Trialsw,1)
    num_ini = size(v_ini,2);
else num_ini = size(m_Trialsw,1);
end
for p=1:num_ini
    col1(size(m_Trialsw,1)-p+1,1) = v_ini(num_ini-p+1);
end
m_Trialsw = horzcat(col1,m_Trialsw);



if ti<0
    t1=-ti;

    %time adjustment
    m_Trialswt = vertcat(m_Trialsw(size(m_Trialsw,1)-floor(t1*srate):end,1:end-1),m_Trialsw(1:floor(t2*srate),2:end));
    v_Timew = linspace(-t1,t2,size(m_Trialswt,1))';
elseif ti==0
    t1=0;
    m_Trialswt = m_Trialsw(1:end,1:end-1);
    v_Timew = linspace(t1,t2,size(m_Trialswt,1))';
else
    t1=ti;
    %time adjustment
    m_Trialswt = m_Trialsw(floor(t1*srate):floor(t2*srate),1:end-1);
    v_Timew = linspace(t1,t2,size(m_Trialswt,1))';
end


end

