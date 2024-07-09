


% Colocar las coordenadas
XSourc=5000;
YSourc=5000;
%cargar el archivo meteorologico
load meteorologia.dat;
% encuentra las distancias
[X,Y]=Dist_XY(meteorologia,XSourc,YSourc);
%introduce la estabilidad y la zona
kst=2;
OpDis='rural';
% encuentra los sigmas
for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if X(j,jj,i)<=0;
              SZout1(j,jj,i)=0;
            else
            SZout1(j,jj,i)=SIGZ(X(j,jj,i),OpDis,kst);
            SYout1(j,jj,i)=SIGY(X(j,jj,i),OpDis,kst);
            end
        end
    end
end
