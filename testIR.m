close all; clear all;

len=10000;
nx=5;
nf=5;

x=rand(1,len);
f=rand(1,len);

xco=[-0.5, 0.2, 0.3, -0.4, 0.4];
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
[xnew,corr,ca,cb]=IR(x,f,5,0);
disp(sprintf('Test 1 corr %2.3f and coef:',corr));
disp(ca');


%Add some randomness to make sure it's not falsely exact
[xnew,corr,ca,cb]=IR(x+rand(1,length(x))*0.01,f,5,0);
disp(sprintf('Test 1+rand corr %2.3f and coef:',corr));
disp(ca');


%Test 2: Only persist, but with a lag
lag=1;
x=rand(1,len);
for i=nx+1:len-lag
    x(i+lag)=xco*x(i-nx:i-1)';
end
[xnew,corr,ca,cb]=IR(x,f,5,0,lag);
disp(sprintf('Test 2 corr %2.3f and coef:',corr));
disp(ca');

[xnew,corr,ca,cb]=IR(x+rand(1,length(x))*0.01,f,5,0,lag);
disp(sprintf('Test 2+rand corr %2.3f and coef:',corr));
disp(ca');


%Test 3: Only response
disp('Test 3: Only response');

x=rand(1,len);
%NOTE: Don't think it should be 6 coef, but currently have IR set to think
%we have current F but one-lag X
for i=nx+1:len
    x(i)=fco*f(i-nx:i-1)';
end
[xnew,corr,ca,cb]=IR(x,f,0,5);
disp(sprintf('Test 3 corr %2.3f and coef:',corr));
disp(cb');


%Test 4: Both, no lags
disp('Test 4: Both, no lags');
x=rand(1,len);
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[xnew,corr,ca,cb]=IR(x,f,5,5);
disp(sprintf('Test 4 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');


%Test 5: Both, 1 lag
disp('Test 5: Both, 1 lag');
x=rand(1,len);
for i=nx+1:len-lag
    x(i+lag)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[xnew,corr,ca,cb]=IR(x,f,5,5,1);
disp(sprintf('Test 5 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');


%Test 6: Both, 1 lag, rand
disp('Test 6: Both, 1 lag, +rand');
x=rand(1,len);
for i=nx+1:len-lag
    x(i+lag)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)';
end
[xnew,corr,ca,cb]=IR(x+rand(1,length(x))*0.01,f+rand(1,length(x))*0.01,5,5,1);
disp(sprintf('Test 6 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');

%Test 7: Both, no lags, constant term
disp('Test 7: Both, no lags');
x=rand(1,len);
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)'+0.5;
end
[xnew,corr,ca,cb,cc]=IR(x,f,5,5);
disp(sprintf('Test 7 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');
disp(cc');

%Test 7: Both, no lags, constant term, guess wrong number
disp('Test 8: Both, no lags, guess wrong number of coeffs');
x=rand(1,len);
for i=nx+1:len
    x(i)=xco*x(i-nx:i-1)'+fco*f(i-nx:i-1)'+0.5;
end
[xnew,corr,ca,cb,cc]=IR(x,f,6,6);
disp(sprintf('Test 8 corr %2.3f and coefs:',corr));
disp(ca');
disp(cb');
disp(cc');