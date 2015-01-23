function map=LagMap(x,f,nf,nlag)


map=zeros(nf,nlag+1);

for i=1:nf
    for j=0:nlag
        
        %[xnew,map(i,j+1)]=IR(x,f,0,i,j);
        [xnew,map(i,j+1)]=IR(x,f,0,i,j);
        if(mod((i*(nlag-1)+j)/(nf*nlag)*100,10)==0)
            fprintf('%d%% .. ',(i*(nlag-1)+j)/(nf*nlag)*100)
        end
    end
end
figure
imagesc(0:nlag,1:nf,map)
colorbar
