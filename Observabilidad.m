function [S1] = Observabilidad(q,W,nx,ny)
matrixdiag=zeros(nx,ny);
Wo=zeros(nx);
n=19;
for i=1:ny
        Wo=Wo+q(i)*W(:,:,i);
end
S1=rank(Wo);
end
