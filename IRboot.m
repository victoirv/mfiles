function [xnew, corr, ca, cb] = IRboot(x,f,numxcoef,numfcoef,lag)
%Usage: [xnew, corr, ca, cb] = IRm(t,x,f,numxcoef,numfcoef)
%Where ca are the x coefficients, cb the f coefficients
%Allows for a matrix of impulses
%***Important: Assumes more data points than impulses***%

if (nargin < 4) || (nargin > 5)
    disp('Usage: [xnew, corr, ca, cb] = IRm(t,x,f,numxcoef,numfcoef)');
    disp('Where ca are the x coefficients, cb the f coefficients');
    disp('***Important: Assumes more data points than impulses***');
    error('');
end
if nargin == 4
    lag=0;
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


for half=1:2
    
    if half==1
    xtrain=x(1:floor(end/2)-1);
    ftrain=f(1:floor(end/2)-1);
    xtest=x(floor(end/2):end);
    ftest=f(floor(end/2):end);
    
    else
    xtrain=x(floor(end/2):end);
    ftrain=f(floor(end/2):end);
    xtest=x(1:floor(end/2)-1);
    ftest=f(1:floor(end/2)-1);
    end   
    
    len=floor(length(xtrain)-predstart);
    
    A=zeros(len,numxcoef+numfcoef*numimpulses+1);
    
    %Must add +i-1 to shift column start point
    for i=1:numxcoef
        A(1:len,i)=xtrain(xstart+i-1:xstart+i-1+len-1);
    end
    for i=1:numfcoef
        for j=1:numimpulses
            A(1:len,i+numxcoef+(j-1)*numfcoef)=ftrain(j,fstart+i-1:fstart+i-1+len-1);
        end
    end
    A(:,end)=1;
    
    b=xtrain(predstart:predstart+len-1);
    A=[A(1:end,:) b'];
    
    for a=1:(numxcoef+numfcoef+1+1) %+1 for mean-normalization coef (cc, column of 1s), and +1 for column of 'b'
        A(isnan(A(:,a)),:)=[];
    end
    b=A(:,end);
    A=A(:,1:end-1);
    
    %Randomly sample rows to draw coefficients from, then average all sets
    nsamp=100;
    coefs=zeros(nsamp,numxcoef+numfcoef+1);
    nrows=floor((rand(1)*size(A,1))/10);
    for i=1:nsamp
        rows=floor(rand(1,nrows)*size(A,1))+1;
        coefs(i,:)=A(rows,:)\b(rows);
    end
    
    coef=mean(coefs);
    ca=coef(1:numxcoef);
    cb=coef(numxcoef+1:end-1);
    cc=coef(end);
    
    
    xnew=zeros(1,length(xtest));
    xnew(1:predstart)=xtest(1:predstart);
    
    %Anywhere f is nan, don't predict, just copy data
    iter=1:(length(ftest));
    iter=iter(iter>=predstart); %Don't use copied variables
    iter=iter(iter<=length(ftest)-lag); %Allow space to predict
    
    for i=iter
        %xnew(i)=(xnew(i-numxcoef:1:i-1)'*ca)+(ftemp(i-numfcoef+1:1:i)'*cb)+cc;
        xnew(i+lag)=(xtest(i-numxcoef:1:i-1)*ca')+(reshape(ftest(:,i-numfcoef:1:i-1),1,[])*cb')+cc;
    end
    
    %xnew(isnan(f))=NaN;
    if half==1
    xnewtotal(floor(end/2):end)=xnew;
    else
    xnewtotal(1:floor(end/2)-1)=xnew;
    end
    
end

xnew=xnewtotal;

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

hit/miss

