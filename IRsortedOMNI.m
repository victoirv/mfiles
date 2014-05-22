function IRsortedOMNI

numcoef=10;
advance=0; %How much future to use to predict the present

OMNI=dlmread('omni2_2003.dat');



%{
[b,m,n]=unique(Ot,'rows');
OMNI=OMNI(m,:);
Ot=(OMNI(:,1)-1980).*(24*365)+OMNI(:,2).*24+OMNI(:,3);
%}

%tnew=t(1):60:t(end)-20;
%{
N=hist(t,tnew);
IN2(N==1,:)=IN(1:sum(N==1),:);
IN2(N<1,:)=9999; %Add missing data as bad so it's pruned later in algorithm
%}

%Denton column values
%10:93, 49 is f1
%13 is dst,
%32 is Vsw, 31 is Bz_sw
%88 is mass density

%x=interp1(t,IN(:,13),tnew)';

%{
vx=IN(:,25);
Bz=IN(:,17);
D=IN(:,41);
I=vx.*1/2.*(abs(Bz)-Bz);
%}
%f=interp1(t,IN(:,32).*(1/2.*(abs(IN(:,31))-IN(:,31))),tnew)';

OMNI(OMNI==9999)=NaN;
OMNI(OMNI==999.9)=NaN;
sorter=OMNI(:,24); %Proton density
f=OMNI(:,25).*1/2.*(abs(OMNI(:,17))-OMNI(:,17));
x=OMNI(:,41);
t=(OMNI(:,1)-1980).*(24*365)+OMNI(:,2).*24+OMNI(:,3);
%{
Unnecessary for omni data
tnew=t(1):1:t(end);
f=interp1(t,f,tnew);
x=interp1(t,x,tnew);
t=tnew;
%}


%sorter=interp1(t,IN(:,88),tnew)';
%{
x=IN(:,13);
f=IN(:,32).*(1/2.*(abs(IN(:,31))-IN(:,31)));
tnew=t;
%}

x(sorter>median(sorter(~isnan(sorter))),:)=[];
f(sorter>median(sorter(~isnan(sorter))),:)=[];
t(sorter>median(sorter(~isnan(sorter))),:)=[];

%IN=IN2;
%{
clear IN2;
clear N;
clear tnew;
%}
run=1;

len=floor(length(t)-numcoef*2);
for start=1:len:length(t)-len-numcoef
    
    A=zeros(len,numcoef*2+1);
    %A=zeros(length(x),numcoef*2+1);
    for i=1:numcoef
        A(1:len,i)=x(i+start-1:i+start+len-2)';
        A(1:len,i+numcoef)=f(i+start-1:i+start+len-2)';
        %        A(1:end-i+1,i)=x(i:end)';
        %A(1:end-i+1,i+numcoef)=f(i:end)';
    end
    
    A(:,end)=1;
    
    b=x(start+numcoef-advance:start-1+len-advance);
    %b=x(1+numcoef-advance:end-advance);
    
    A=[A(1:end-numcoef,:) b];
    
    for a=1:numcoef*2+1
        A(isnan(A(:,a)),:)=[];
    end
    %{
    if(size(A)<500)
        fprintf('Not enough coef: %d\n',size(A));
        errormatcoef(run)=t(start);
        errormatval(run)=0.0;
        run=run+1;
        continue
    end
    A=A(1:500,:); %Make each segment the same length
    %}
    b=A(:,end);
    A=A(:,1:end-1);
    
    
    coef=A(1:end,:)\b;
    %{
    m=(length(coef)-1);
    b=coef(1:end-1);
    a=b.*0;
    %}
    m=(length(coef)-1)/2;
    ca=coef(1:(end-1)/2);
    cb=coef((end-1)/2+1:end-1);
    cc=coef(end);
    
    xtemp=x;
    ftemp=f;
    xtemp(isnan(f))=NaN;
    ftemp(isnan(x))=NaN;
    xtemp=xtemp(~isnan(xtemp));
    ftemp=ftemp(~isnan(ftemp));
  
    xtest=zeros(length(xtemp),1);
    xtest(1:m)=xtemp(1:m);
    for i=m+1:length(xtemp)-advance-1
        xtest(i)=(xtest(i-m+advance:1:i-1+advance)'*ca)+(ftemp(i-m+advance:1:i-1+advance)'*cb)+cc;
    end
    
    corr=corrcoef(xtest(m:end),xtemp(m:end)); %Ignore first added bit
    corr=corr(1,2);
    
    fprintf('%f %d %d\n',corr,numcoef,start);
    %errormatcoef(run)=numcoef;
    errormatcoef(run)=t(start);
    errormatval(run)=corr;
    run=run+1;

end

plot(errormatcoef,errormatval,'k+');


stuff=1; %For a breakpoint