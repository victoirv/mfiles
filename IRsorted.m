function IRsorted
numcoef=10;
advance=0; %How much future to use to predict the present

IN=dlmread('/media/D/Denton/massdensitypruned.txt');
fid=fopen('/media/D/Denton/headers.txt');
headers=textscan(fid,'%s');
fclose(fid);

OMNI=dlmread('omni2_test.dat');



%IN=sortrows(IN,6);
%IN(IN(:,1)~=6,:)=[]; %Only GOES 6

%Sort out unique values
Dt=(IN(:,2)-1980).*(24*60*365)+IN(:,3).*(24*60)+IN(:,4).*60+(IN(:,5)-5);

[b,m,n]=unique(Dt,'rows');
IN=IN(m,:);
Dt=(IN(:,2)-1980).*(24*60*365)+IN(:,3).*(24*60)+IN(:,4).*60+(IN(:,5)-5);

%Denton notable column values
%10:93, 49 is f1
%13 is dst,
%32 is Vsw, 31 is Bz_sw
%88 is mass density


OMNI(OMNI==9999)=NaN;
OMNI(OMNI==999.9)=NaN;


Ot=(OMNI(:,1)-1980).*(24*365)+OMNI(:,2).*24+OMNI(:,3);
[b,m,n]=unique(Ot,'rows');
OMNI=OMNI(m,:);
Ot=(OMNI(:,1)-1980).*(24*365)+OMNI(:,2).*24+OMNI(:,3);

tnew=Ot;
t=tnew;

for sortx=10:93
t=tnew;
run=sortx-9;
x=OMNI(:,41);
f=OMNI(:,25).*1/2.*(abs(OMNI(:,17))-OMNI(:,17));


%Sort by mass density 

%sorter=interp1(t./60,IN(:,88),tnew,'nearest')';
%sorter=OMNI(:,24);
sorter=interp1(Dt./60,IN(:,sortx),tnew,'nearest')';

x(sorter<median(sorter(~isnan(sorter))),:)=[];
f(sorter<median(sorter(~isnan(sorter))),:)=[];
t(sorter<median(sorter(~isnan(sorter))),:)=[];


[xnew, corr, ca, cb]=IR(t,x,f,numcoef);
errormatcoef(run)=sortx;
errormatval(run)=corr;

close all
figure
plot(0:length(ca)-advance-1,flipud(ca)/norm(flipud(ca)),'b');
hold on

x=OMNI(:,41);
f=OMNI(:,25).*1/2.*(abs(OMNI(:,17))-OMNI(:,17));
t=tnew;
x(sorter>median(sorter(~isnan(sorter))),:)=[];
f(sorter>median(sorter(~isnan(sorter))),:)=[];
t(sorter>median(sorter(~isnan(sorter))),:)=[];


[xnew2, corr2, ca2, cb2]=IR(t,x,f,numcoef);
%errormatcoef(run)=sortx;
%errormatval(run)=corr2;

diff=abs(corr2-corr);

plot(0:length(ca2)-advance-1,flipud(ca2)/norm(flipud(ca2)),'r');

legend(sprintf('upper corr %1.3f',corr),sprintf('lower corr %1.3f',corr2));
title(sprintf('%1.3f - %s',diff,cell2mat(headers{1}(sortx))));

saveas(gcf,sprintf('IRcoef/%1.3f-%s.png',diff,cell2mat(headers{1}(sortx))));
close all;



end

plot(errormatcoef,errormatval,'k+');


stuff=1; %For a breakpoint