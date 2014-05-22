function map=CorrMap(x,f,nx,ny,lag)


map=zeros(nx+1,ny+1);

for i=0:nx
    for j=0:ny
        if(i~=0 && j~=0)
        [xnew,map(i+1,j+1)]=IR(x,f,i,j,lag);
        end
    end
end
figure
imagesc(0:nx,0:ny,map)
        