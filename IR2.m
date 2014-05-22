function IR2

numcoef=24; %each
advance=0; %How much future to use to predict the present 


IN=dlmread('/media/C/Downloads/massdensitypruned.txt');
fid=fopen('/media/C/Downloads/headers.txt');
headers=textscan(fid,'%s');
fclose(fid);

IN=sortrows(IN,2);
IN(IN(:,1)~=6,:)=[]; %Only use satellite 6

for a=1:93
    IN(IN(:,a)==9999,:)=[];
    IN(IN(:,a)==999.9,:)=[];
end

SAT=IN(:,1);
x=IN(:,88);
t=IN(:,6);

mincorr=0;

for numcoef=3:60
for input=49 %10:93, 49 is f1
    if(input==88)
        continue
    end
f=IN(:,input);

usefx=0;


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
    b=coef(1:end-1);
    a=b.*0;
    c=coef(end);
end

xtest=zeros(length(x),1);
xtest(1:m)=x(1:m);
for i=m+1:length(x)-advance-1
    xtest(i)=(xtest(i-m+advance:1:i-1+advance)'*a)+(f(i-m+advance:1:i-1+advance)'*b)+c;
end

corr=corrcoef(x(m:end),xtest(m:end)); %Ignore first added bit
corr=corr(1,2);

fprintf('%s %f %d\n',cell2mat(headers{1}(input)),corr,numcoef);
end
%{
if(corr>mincorr)
   mincorr=corr;
   %fprintf('%s %d\n',cell2mat(headers{1}(input)),corr);
   figure
   one=plot(t,x);
   hold on;
   two=plot(t,xtest,'k');
   legend([one, two],'Actual','Mine');
   title(cell2mat(headers{1}(input)));
   hold off;
   figure
   plot(-advance:1:length(b)-advance-1,flipud(b)/norm(flipud(b)));
   title(cell2mat(headers{1}(input)));
   xlabel('Hours before current');
   ylabel('Normalized Coefficient');
end
%}

end


stuff=1; %For a breakpoint