function corrmat = CorrWithLags(x,f,numxcoef,numfcoef,lags)
%corrmat = CorrWithLags(x,f,numxcoef,numfcoef,lags)
%Where lags is the maximum number to check to


corrmat=zeros(1,lags+1);
for lag=0:lags
    [xnew,corr]=IR(x,f,numxcoef,numfcoef,lag);
    corrmat(lag+1)=corr;
end

figure
plot(0:lags,corrmat)
xlabel('Lags')
ylabel('Correlation')
title(sprintf('IR Correlation with lags, %d persist coef, %d impulse coef',numxcoef,numfcoef))