function [ca, cb, cc, xnew, corr, eff, xnew1, xnew2] = IRboot(x,f,numxcoef,numfcoef,lag,advance)
%Usage: [xnew, corr, ca, cb] = IRm(t,x,f,numxcoef,numfcoef)
%Where ca are the x coefficients, cb the f coefficients
%Allows for a matrix of impulses
%***Important: Assumes more data points than impulses***%

if (nargin < 4) || (nargin > 6)
    disp('Usage: [xnew, corr, ca, cb] = IRm(t,x,f,numxcoef,numfcoef)');
    disp('Where ca are the x coefficients, cb the f coefficients');
    disp('***Important: Assumes more data points than impulses***');
    error('');
end
if nargin == 4
    lag=0;
    advance=0;
end
if nargin == 5
    advance=0;
end

%Make x and f row vectors for standardization purposes
if(length(x)~=size(x,2))
    x=x';
end
if(length(f)~=size(f,2))
    f=f';
end

predstart=max(numxcoef,numfcoef)+1+lag;

xstart=predstart-numxcoef-lag;
fstart=predstart-numfcoef-lag;
numimpulses=min(size(f));
datalen=length(x);

xnewtotal=x;

[ca1,cb1,cc1]=IR(x(1:floor(end/2)-1),f(1:floor(end/2)-1),numxcoef,numfcoef,lag,advance);
[ca2,cb2,cc2]=IR(x(floor(end/2):end),f(floor(end/2):end),numxcoef,numfcoef,lag,advance);


%xnew=xnewtotal;
xnew=x;

ca=(ca2+ca1)./2;
cb=(cb2+cb1)./2;
cc=(cc2+cc1)./2;

    iter=1:(length(x));
    iter=iter(iter>=predstart); %Don't use copied variables
    iter=iter(iter<=length(x)-lag); %Allow space to predict
    
    for i=iter
        %xnew(i)=(xnew(i-numxcoef:1:i-1)'*ca)+(ftemp(i-numfcoef+1:1:i)'*cb)+cc;
        xnew(i+lag)=(x(i-numxcoef:1:i-1)*ca)+(reshape(f(:,i-numfcoef:1:i-1),1,[])*cb)+cc;
        xnew1(i+lag)=(x(i-numxcoef:1:i-1)*ca1)+(reshape(f(:,i-numfcoef:1:i-1),1,[])*cb1)+cc1;
        xnew2(i+lag)=(x(i-numxcoef:1:i-1)*ca2)+(reshape(f(:,i-numfcoef:1:i-1),1,[])*cb2)+cc2;
        
    end



%Calculate correlation here to save program from needing to strip NaNs
skip=(isnan(xnew) | isnan(x));
skip(1:predstart+lag)=1;
corr=corrcoef(xnew(~skip),x(~skip)); %Ignore first added bit
corr=corr(1,2);

%HSS?
thresh=quantile(xnew,0.9);
hit=0;
miss=0;
for test=2:length(xnew)
   if(x(test)>=thresh && x(test-1)<thresh)
       if(xnew(test)>=thresh)
           hit=hit+1;
       else
           miss=miss+1;
       end
   end
end

eff=hit/miss;
xnew=xnew';

