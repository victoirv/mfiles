function [xnew, corr, ca, cb] = IRf(t,x,f,numfcoef)
%Usage: [xnew, corr, ca, cb] = IRf(t,x,f,numfcoef)
%Where cb are the f coefficients, ca=0

if nargin ~= 4 
    disp('Usage: [xnew, corr, ca, cb] = IR(t,x,f,numfcoef)');
    disp('Where ca are the x coefficients, cb the f coefficients');
    return
end

len=floor(length(t)-2*numfcoef);
start=1;
advance=0;
    
A=zeros(len,numfcoef+1);
%A=zeros(length(x),numfcoef*2+1);
for i=1:numfcoef
    A(1:len,i)=f(i+start-1:i+start+len-2)';
end

A(:,end)=1;

b=x(start+numfcoef-1-advance:start-2+len-advance);
%b=x(1+numfcoef-advance:end-advance);

A=[A(1:end-numfcoef,:) b];

for a=1:numfcoef+1
    A(isnan(A(:,a)),:)=[];
end

b=A(:,end);
A=A(:,1:end-1);

%Get coefficients
coef=A(1:end,:)\b;

m=(length(coef)-1);
cb=coef(1:(end-1));
ca=zeros(length(coef)-1,1);
cc=coef(end);

xtemp=x;
ftemp=f;
xtemp(isnan(f))=NaN;
ftemp(isnan(x))=NaN;
xtemp=xtemp(~isnan(xtemp));
ftemp=ftemp(~isnan(ftemp));

xnew=zeros(length(xtemp),1);
xnew(1:m)=xtemp(1:m);
for i=m+1:length(xtemp)-advance
    xnew(i)=(ftemp(i-m+advance+1:1:i+advance)'*cb)+cc;
end

%Calculate correlation here to save program from needing to strip NaNs
corr=corrcoef(xnew(~isnan(xnew)),xtemp(~isnan(xnew))); %Ignore first added bit
corr=corr(1,2);
