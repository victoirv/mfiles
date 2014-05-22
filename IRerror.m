function IRerror

numcoef=24; %each
advance=0; %How much future to use to predict the present


IN=dlmread('/media/D/Denton/massdensitypruned.txt');
fid=fopen('/media/D/Denton/headers.txt');
headers=textscan(fid,'%s');
fclose(fid);

%IN=sortrows(IN,6);
IN(IN(:,1)~=6,:)=[]; %Only GOES 6

%Sort out unique values
t=(IN(:,2)-1980).*(24*60*365)+IN(:,3).*24*60+IN(:,4).*60+IN(:,5);
[b,m,n]=unique(t,'rows');
IN=IN(m,:);
IN2=IN.*0;

tnew=t(1):10:t(end);
N=hist(t,tnew);
IN2(N==1,:)=IN(1:sum(N==1),:);
IN2(N~=1,:)=9999; %Add missing data as bad so it's pruned later in algorithm
IN=IN2;
t=tnew;
clear IN2;
clear N;
clear tnew;

run=1;

len=floor(length(t)/200);
x=IN(:,88);

for numcoef=10
    for input=49 %10:93, 49 is f1
        f=IN(:,input);
        for start=1:len:length(t)-len
            
            A=zeros(len,numcoef+1);
            for i=1:numcoef
                A(1:len,i)=f(i+start:i+start+len-1)';
            end
            A(:,end)=1;
            
            b=x(start+numcoef-advance:start-1+len-advance);
            
            A=[A(1:end-numcoef,:) b];
            for a=1:numcoef+2
              A(A(:,a)==9999,:)=[];
              A(A(:,a)==999.9,:)=[];
            end
            if(size(A)<100)
                fprintf('Not enough coef: %d\n',size(A));
                errormatcoef(run)=t(start);
                errormatval(run)=0.0;
                run=run+1;
                continue
            end
            A=A(1:100,:); %Make each segment the same length
            b=A(:,end);
            A=A(:,1:end-1);
            
            
            coef=A(1:end,:)\b;
            
            m=(length(coef)-1);
            b=coef(1:end-1);
            a=b.*0;
            c=coef(end);
            
            xtemp=x(x~=9999);
            xtest=zeros(length(xtemp),1);
            xtest(1:m)=xtemp(1:m);
            ftemp=f(f~=9999);
            for i=m+1:length(xtemp)-advance-1
                xtest(i)=(xtest(i-m+advance:1:i-1+advance)'*a)+(ftemp(i-m+advance:1:i-1+advance)'*b)+c;
            end
            
            corr=corrcoef(xtemp(m:end),xtest(m:end)); %Ignore first added bit
            corr=corr(1,2);
            
            fprintf('%s %f %d %d\n',cell2mat(headers{1}(input)),corr,numcoef,start);
            %errormatcoef(run)=numcoef;
            errormatcoef(run)=t(start);
            errormatval(run)=corr;
            run=run+1;
        end
    end    
end

plot(errormatcoef,errormatval,'k+');


stuff=1; %For a breakpoint