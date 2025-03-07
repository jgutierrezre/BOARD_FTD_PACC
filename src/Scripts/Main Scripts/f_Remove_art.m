function  m_Tmod = f_Remove_art(m_Trials,option,L1,L2,p)
if option==0
    m_Tmod = m_Trials;
elseif option ==1
    m_Tmod = m_Trials(L1+1:end,:);
else
    v1=m_Trials(L1+1:p-1,:);
    v2=m_Trials(L2+1:end,:);
    
    m_Tmod=vertcat(v1,v2);
end
end

