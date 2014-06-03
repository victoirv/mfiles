function [xnew, xnew2, corr, corr2, ca, ca2, cb, cb2, cc, cc2] = IRsort(x,f,sorter,numxcoef,numfcoef,lag)
%Usage: [xnew, xnew2, corr, corr2, ca, ca2, cb, cb2, cc, cc2] = IRsort(x,f,sorter,numxcoef,numfcoef,lag)
%Where ca are the x coefficients, cb the f coefficients
%Allows for a matrix of impulses
%***Important: Assumes more data points than impulses***%
 
if (nargin < 5) || (nargin > 6)
    disp('Usage: [xnew, xnew2, corr, corr2, ca, ca2, cb, cb2, cc, cc2] = IRsort(x,f,sorter,numxcoef,numfcoef,lag)');
    disp('Where ca are the x coefficients, cb the f coefficients');
    disp('***Important: Assumes more data points than impulses***');
    error('');
end
if nargin == 5
    lag=0;
end

%Make x and f row vectors for standardization purposes
if(length(x)~=size(x,2))
    x=x';
end
if(length(f)~=size(f,2))
    f=f';
end
if(length(sorter)~=size(sorter,2))
    sorter=sorter';
end

predstart=max(numxcoef,numfcoef)+1+lag;

xstart=predstart-numxcoef-lag;
fstart=predstart-numfcoef-lag;


len=floor(length(x)-predstart);

numimpulses=min(size(f));
    
A=zeros(len,numxcoef+numfcoef*numimpulses+1+1); %+1 for mean-subtraction, +1 for sorter

%Must add +i-1 to shift column start point
for i=1:numxcoef
   A(1:len,i)=x(xstart+i-1:xstart+i-1+len-1); 
end
for i=1:numfcoef
    for j=1:numimpulses
        A(1:len,i+numxcoef+(j-1)*numfcoef)=f(j,fstart+i-1:fstart+i-1+len-1);
    end
end

%Add sorter. A=[xs... fs... 1 sort]
for i=1:len
   A(i,end)=mean(sorter(xstart+i:xstart+numxcoef+i));
end

A(:,end-1)=1;

b=x(predstart:predstart+len-1);
A=[A(1:end,:) b'];

%A=[xs... fs... 1 sort b]

for a=1:(numxcoef+numfcoef+1+1+1) %+1 for mean-normalization coef (cc, column of 1s), +1 for sorter, and +1 for column of 'b'
    A(isnan(A(:,a)),:)=[];
end

%Now separate into two parts via sorter
sorter=A(:,end-1);
A2=A;
split=mean(sorter);
A(sorter>split,:)=[];
A2(sorter<=split,:)=[];


b=A(:,end);
A=A(:,1:end-2);
b2=A2(:,end);
A2=A2(:,1:end-2);


coef=A(1:end,:)\b;
coef2=A2(1:end,:)\b2;

ca=coef(1:numxcoef);
cb=coef(numxcoef+1:end-1);
cc=coef(end);
ca2=coef2(1:numxcoef);
cb2=coef2(numxcoef+1:end-1);
cc2=coef2(end);

xtemp=x;
ftemp=f;

xnew=zeros(1,length(x));
xnew(1:predstart)=xtemp(1:predstart);
xnew2=xnew;


%Anywhere f is nan, don't predict, just copy data
iter=1:(length(f));
iter=iter(iter>=predstart); %Don't use copied variables
iter=iter(iter<=length(f)-lag); %Allow space to predict 

for i=iter
    %xnew(i)=(xnew(i-numxcoef:1:i-1)'*ca)+(ftemp(i-numfcoef+1:1:i)'*cb)+cc;
    xnew(i+lag)=(xtemp(i-numxcoef:1:i-1)*ca)+(reshape(ftemp(:,i-numfcoef:1:i-1),1,[])*cb)+cc;
    xnew2(i+lag)=(xtemp(i-numxcoef:1:i-1)*ca2)+(reshape(ftemp(:,i-numfcoef:1:i-1),1,[])*cb2)+cc2;
    %xnew(i+lag)=(xnew(i-numxcoef:1:i-1)*ca)+(reshape(ftemp(:,i-numfcoef:1:i-1),1,[])*cb)+cc;
end

%xnew(isnan(f))=NaN;





%Calculate correlation here to save program from needing to strip NaNs
skip=(isnan(xnew) | isnan(xtemp));
skip(1:predstart+lag)=1;
corr=corrcoef(xnew(~skip),xtemp(~skip)); %Ignore first added bit
corr=corr(1,2);
corr2=corrcoef(xnew2(~skip),xtemp(~skip)); %Ignore first added bit
corr2=corr2(1,2);
