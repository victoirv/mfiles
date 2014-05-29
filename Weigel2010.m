close all;clear all;
load('OMNI_OMNI2_merged')
%getyears=Year>=2000;
%getyears=getyears+(Year<=2011);
%getyears=(getyears==2);
%KP=Kp_index(getyears);
VBS=1/2*Plasma_bulk_speed.*(abs(Bz_GSM)-Bz_GSM);
%VBS=VBS(getyears);
%DST=Dst_index(getyears);

N=50;
Na=0;
Nb=36;

if(Na>0)
cas=zeros(N,Na);
ca2s=cas;
cbs=zeros(N,Nb);
cb2s=cbs;
for i=1:N
   rstart=floor(rand(1)*length(x));
   rend=rstart+floor(rand(1)*(length(x)-rstart));
   [xnew,xnew2,corr,corr2,cas(i,:),ca2s(i,:),cbs(i,:),cb2s(i,:)]=IRsort(Dst_index(rstart:rend),VBS(rstart:rend),Ion_density(rstart:rend),Na,Nb);
end
else
cbs=zeros(N,Nb);
cb2s=cbs;
for i=1:N
   rstart=floor(rand(1)*length(VBS));
   rend=rstart+floor(rand(1)*(length(VBS)-rstart));
   [xnew,xnew2,corr,corr2,ca,ca2,cbs(i,:),cb2s(i,:)]=IRsort(Dst_index(rstart:rend),VBS(rstart:rend),Ion_density(rstart:rend),Na,Nb);
   if(mod(i/N*100,10)==0)
       fprintf('%d%% ... ',i/N*100)
   end
end
    
end

cb=mean(cbs);
cb2=mean(cb2s);



figure; plot(flipud(cb)); hold on; plot(flipud(cb2),'r'); legend('Low Density','High Density','Location','SouthEast')
plot(cb+2*std(cbs),'b:');plot(cb-2*std(cbs),'b:');
plot(cb2+2*std(cb2s),'r:');plot(cb2-2*std(cb2s),'r:');
ylabel('DST response to VBs')
xlabel('Time Lags')
title(sprintf('N:%d Nx:%d Nf:%d',N,Na,Nb))
print -dpng dstvbs.png