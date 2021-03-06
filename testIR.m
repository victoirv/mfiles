close all; clear all;

len=10000;
nx=5;
nf=5;

x=rand(1,len);
f=rand(1,len);

xco=[-0.5, 0.4, 0.1, -0.4, 0.4];
fco=[0.2, 0.3, -0.6, -0.4, 0.5];
disp('Real x coefficients:');
disp(xco);
disp('Real f coefficients:');
disp(fco);


%Test 1, only persist
disp('Test 1: Only persist');
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)';
end

%Correlation should be exact
[ca,cb,cc,xnew,corr]=IR(x,f,5,0);
disp(sprintf('Test 1 corr %2.3f and coef:',corr));
disp(ca');


%Add some randomness to make sure it's not falsely exact
[ca,cb,cc,xnew,corr]=IR(x+rand(1,length(x))*0.01,f,5,0);
disp(sprintf('Test 1+rand corr %2.3f and coef:',corr));
disp(ca');


%Test 2: Only persist, but with a lag
lag=1;
x=rand(1,len);
for i=nx+1:len-lag
    x(i+lag)=xco*x(i-nx:i-1)';
end
[ca,cb,cc,xnew,corr]=IR(x,f,5,0,lag);
disp(sprintf('Test 2 corr %2.3f and coef:',corr));
disp(ca');

[ca,cb,cc,xnew,corr]=IR(x+rand(1,length(x))*0.01,f,5,0,lag);
disp(sprintf('Test 2+rand corr %2.3f and coef:',corr));
disp(ca');


%Test 3: Only response
disp('Test 3: Only response');

x=rand(1,len);
for i=nx+1:len
    x(i)=fco*f(i-nx:i-1)';
end
[ca,cb,cc,xnew,corr]=IR(x,f,0,5);
disp(sprintf('Test 3 corr %2.3f and coef:',corr));
disp(cb');


%Test 4: Both, no lags
disp('Test 4: Both, no lags');
x=rand(1,len);
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[ca,cb,cc,xnew,corr]=IR(x,f,5,5);
disp(sprintf('Test 4 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');


%Test 5: Both, 1 lag
disp('Test 5: Both, 1 lag');
x=rand(1,len);
for i=nx+1:len-lag
    x(i+lag)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[ca,cb,cc,xnew,corr]=IR(x,f,5,5,1);
disp(sprintf('Test 5 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');


%------------------.=^=.=^=.\/\/\/\/\/\/ testing comment dividers
disp('Test 6: Both, 1 lag, +rand');
x=rand(1,len);
for i=nx+1:len-lag
    x(i+lag)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[ca,cb,cc,xnew,corr]=IR(x+rand(1,length(x))*0.01,f+rand(1,length(x))*0.01,5,5,1);
disp(sprintf('Test 6 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');

%------------------
disp('Test 7: Both, no lags, constant coef term');
x=rand(1,len);
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)'+0.5;
end
[ca,cb,cc,xnew,corr]=IR(x,f,5,5);
disp(sprintf('Test 7 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');
disp(cc');

%------------------
disp('Test 8: Both, no lags, guess wrong number of coeffs');
x=rand(1,len);
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)'+0.5;
end
[ca,cb,cc,xnew,corr]=IR(x,f,6,6);
disp(sprintf('Test 8 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');
disp(cc');



%------------------
disp('Test 9: Sorted IR, but with same coefs');
x=rand(1,len);
sorter=1:len;
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[xnew,xnew2,corr,corr2,ca,ca2,cb,cb2,cc]=IRsort(x,f,sorter,5,5);
disp(sprintf('Test 9 corr %2.3f and coefs:',corr));
disp(ca');
disp(ca2');
disp(cb');
disp(cb2');
disp(cc');

%------------------
disp('Test 10: Sorted IR, different x coef, overlap not accounted for');
x=rand(1,len);
sorter=1:len;
for i=nx+1:len/2
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
    %x2(i)=fliplr(xco)*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
xco2=[-0.3, 0.1, 0.1, -0.2, 0.3];
for i=(len/2+1):len
    x(i)=xco2*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
%x(sorter>mean(sorter))=x2(sorter>mean(sorter));
[xnew,xnew2,corr,corr2,ca,ca2,cb,cb2,cc]=IRsort(x,f,sorter,5,5);
disp(sprintf('Test 10 corr %2.3f and coefs:',corr));
disp(ca');
disp(ca2');
disp(cb');
disp(cb2');
disp(cc');