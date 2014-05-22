function IRold

numcoef=10; %each
advance=0; %How much future to use to predict the present

IN=dlmread('omni2_test.dat');
f=zeros(1,length(IN(:,1)));
vx=IN(:,25);
Bz=IN(:,17);
D=IN(:,41);
t=IN(:,1).*(24*365)+IN(:,2).*24+IN(:,3); %hours

%Prune bad values
D(Bz==999.9)=[];
t(Bz==999.9)=[];
vx(Bz==999.9)=[];
Bz(Bz==999.9)=[];

D(vx==9999)=[];
t(vx==9999)=[];
Bz(vx==9999)=[];
vx(vx==9999)=[];   

I=vx.*1/2.*(abs(Bz)-Bz);
f=I;
x=D;  

%{
%Use Denton Data
IN=dlmread('/media/D/Denton/massdensitypruned.txt');
%IN=sortrows(IN,2);
%IN(IN(:,1)~=6,:)=[]; %Only use satellite 2
SAT=IN(:,1);
x=IN(:,88);
f=IN(:,49);
t=IN(:,6);
x=IN(:,13);
f=IN(:,32).*(1/2.*(abs(IN(:,31))-IN(:,31)));
%}


usefx=1;

if usefx==1
    A=zeros(length(x),numcoef*2+1);
    for i=1:numcoef
        A(1:end-i+1,i)=x(i:end);
        A(1:end-i+1,i+numcoef)=f(i:end);
    end
    A(:,end)=1;
    
else
    A=zeros(length(x),numcoef+1);
    for i=1:numcoef
        A(1:end-i+1,i)=f(i:end)';
    end
    A(:,end)=1;
end

b=x(1+numcoef-advance:end-advance);
coef=A(1:end-numcoef,:)\b;

if usefx==1
    m=(length(coef)-1)/2;   
    a=coef(1:(end-1)/2);
    b=coef((end-1)/2+1:end-1);
    c=coef(end);
else
    m=(length(coef)-1);
    a=coef(1:end-1).*0;
    b=coef(1:end-1);
    c=coef(end);
end

a'
b'

xtest=zeros(length(x),1);
xtest(1:m)=x(1:m);
for i=m+1:length(x)-advance-1
    xtest(i)=(xtest(i-m+advance:1:i-1+advance)'*a)+(f(i-m+advance:1:i-1+advance)'*b)+c;
end

corr=corrcoef(xtest(m:end),x(m:end));
corr

one=plot(t,x);
hold on;
two=plot(t,xtest,'k');
legend([one, two],'Actual','Mine');

figure
plot(-advance:1:length(b)-advance-1,flipud(b)/norm(flipud(b)));
title('Impulse response');
xlabel('Hours before current');
ylabel('Normalized Coefficient');

stuff=1;