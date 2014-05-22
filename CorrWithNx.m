function corrmat = CorrWithNx(x,f,numxcoef,numfcoef,lag)
%function corrmat = CorrWithNx(x,f,numxcoef,numfcoef,lag)

corrmat=zeros(1,numxcoef);
camat=zeros(numxcoef);
for nx=0:numxcoef
    [xnew,corr,ca]=IR(x,f,nx,numfcoef,lag);
    corrmat(nx+1)=corr;
    
    camat(nx+1,:)=[ca' zeros(1,numxcoef-nx)];
    
end

camat
plot(0:numxcoef,corrmat)