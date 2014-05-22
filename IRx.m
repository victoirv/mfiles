function [xnew, corr, ca, cb] = IRx(t,x,f,numxcoef,numfcoef)
%Usage: [xnew, corr, ca, cb] = IRx(t,x,f,numxcoef)
%Where ca are the x coefficients, cb the f coefficients

if nargin ~= 4 
    disp('Usage: [xnew, corr, ca, cb] = IRx(t,x,f,numxcoef)');
    disp('Where ca are the x coefficients, cb the f coefficients');
    return
end

len=floor(length(t)-2*numxcoef);
start=1;
advance=0;
    
A=zeros(len,numxcoef+1);
for i=1:numxcoef
    A(1:len,i)=x(i+start-1:i+start+len-2)';
end

A(:,end)=1;

b=x(start+numxcoef-advance:start-1+len-advance);
%b=x(1+numxcoef-advance:end-advance);

A=[A(1:end-numxcoef,:) b];

for a=1:numxcoef+1
    A(isnan(A(:,a)),:)=[];
end

b=A(:,end);
A=A(:,1:end-1);

%Get the coefficients
coef=A(1:end,:)\b;

m=(length(coef)-1);
ca=coef(1:(end-1));
cb=zeros(length(coef)-1,1); %For returning
cc=coef(end);


%Despite only using x, get rid of bad f values to match data samples in other tests
xtemp=x;
%xtemp(isnan(f))=NaN;
%xtemp=xtemp(~isnan(xtemp));
xtemp(isnan(f))=NaN;

xnew=zeros(length(x),1);
xnew(1:m)=xtemp(1:m);

iter=1:(length(f));
xnew(isnan(f))=x(isnan(f));
cutout=1:length(f);
cutout(isnan(f))=NaN;
for i=1:m
    cutout(circshift(isnan(f),[0 i]))=NaN; %Cut out all that would be used
end
iter=iter(~isnan(cutout)); %Don't predict using NaNs 
iter=iter(iter>=m+1); %Don't use copied variables
iter=iter(1:end-advance); %Only matters if advance>0


for i=iter
    xnew(i)=(xnew(i-m+advance:1:i-1+advance)'*ca)+cc;
end

%Calculate correlation here to save program from needing to strip NaNs
corr=corrcoef(xnew(~isnan(f)),xtemp(~isnan(f))); %Ignore first added bit
corr=corr(1,2);
