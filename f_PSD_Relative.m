function [pow,powr,pr,p] = f_PSD_Relative(v_Data,srt,srate,F1,F2,GF)

srate =srate/srt;
v_Data= downsample(v_Data,srt);

[pxx,f] = pwelch(v_Data,[],[],[],srate);

if GF==1
    pxx_av = pxx;
else
    pxx_av = mean(pxx');
end

ad1=find(f>=F1);
ad2=find(f>=F2);

id1=ad1(1);
id2=ad2(1);



l1=find(f>=0.5);
lim1=l1(1);
l2=find(f>=500);
lim2=l2(1);
Ptot= sum(pxx_av(lim1:lim2));
powr = sum(pxx_av(id1:id2))/Ptot; 
pow = sum(pxx_av(id1:id2));
 


d1=0.5;
d2=3.9;
t1=4;
t2=7.9;
a1=8;
a2=11.9;
b1=12;
b2=24.9;
g1=25;
g2=59;
g3=61;
g4=120;
HFO11=121;
HFO12=250;
HFO21=251;
HFO22=500;

ad1=find(f>=d1);
ad2=find(f>=d2);
aa1=find(f>=a1);
aa2=find(f>=a2);
ab1=find(f>=b1);
ab2=find(f>=b2);
at1=find(f>=t1);
at2=find(f>=t2);
ag1=find(f>=g1);
ag2=find(f>=g2);
ag3=find(f>=g3);
ag4=find(f>=g4);
aHFO11=find(f>=HFO11);
aHFO12=find(f>=HFO12);
aHFO21=find(f>=HFO21);
aHFO22=find(f>=HFO22);

id1=ad1(1);
id2=ad2(1);
ia1=aa1(1);
ia2=aa2(1);
ib1=ab1(1);
ib2=ab2(1);
it1=at1(1);
it2=at2(1);
ig1=ag1(1);
ig2=ag2(1);
ig3=ag3(1);
ig4=ag4(1);
iHFO11=aHFO11(1);
iHFO12=aHFO12(1);
iHFO21=aHFO21(1);
iHFO22=aHFO22(1);




pr(1) = round(sum(pxx(id1:id2))/Ptot,4);
pr(2) = round(sum(pxx(it1:it2))/Ptot,4);
prg1 = round(sum(pxx(ig1:ig2))/Ptot,4);
prg2 = round(sum(pxx(ig3:ig4))/Ptot,4);
pr(3) = round(sum(pxx(ia1:ia2))/Ptot,4);
pr(4) = round(sum(pxx(ib1:ib2))/Ptot,4);
pr(5)=prg1+prg2;
pr(6) = round(sum(pxx(iHFO11:iHFO12))/Ptot,4);
pr(7) = round(sum(pxx(iHFO21:iHFO22))/Ptot,4);

p(1) = round(sum(pxx(id1:id2)),4);
p(2) = round(sum(pxx(it1:it2)),4);
pg1 = round(sum(pxx(ig1:ig2)),4);
pg2 = round(sum(pxx(ig3:ig4)),4);
p(3) = round(sum(pxx(ia1:ia2)),4);
p(4) = round(sum(pxx(ib1:ib2)),4);
p(5)=pg1+pg2; 
p(6) = round(sum(pxx(iHFO11:iHFO12)),4);
p(7) = round(sum(pxx(iHFO21:iHFO22)),4);
end


