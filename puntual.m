

function puntual(M,I,N)
% PUNTO, CALCULA LA CONCENTRACION HORARIA PARA FUENTES PUNTUALES
%  SUBRUTINAS QUE SON EMPLEADAS
%  FPLUME --- ALTURA FINAL DE SOBREELEVACION DE LA PLUMA (BRIGGS)
%  XPLUME --- ALTURA EFECTIVA DE LEVACION A UNA DISTANCIS X (BRIGGS)
%  RCONCP --- DETERMINA LA CONCENMTRACION RELATIVA EN UN RECEPTOR DADA 
%             SU DISTANCIA DESDE UNA FUENTE
%  FUNCIONES LLAMADAS
%  XVY ---- CALCULA LA DISTNACIA NECESARIA 
%  
% XVZ
% QP            EMISION  (G/SEG)
% HPP           ALTURA FISICA DE LA CHIMENEA
% TSP           TEMPERATURA DE SALIDAD DEL GAS
% VSP           VELOCIDAD DE SALIDAD DEL GAS
% DP            DIAMETRO DE LA CHIMENEA
% VFP           VOLUMEN DEL FLUJO DEL GAS
% RQP           COORDENADAS ESTE DE LA FUENTE
% SQR           COORDENADAS NORTE DE LA CHIMENEA
% SYOP          SIGMA Y
% SZOP          SIGMA Z
vs=10;h=30;d=1.2; Tg=400;
puz=[0.15;0.15;0.20;0.25;0.40;0.60] ;
% for nm=1:lengeth(meteorologia(:,1))
angulo =meteorologia(:,5);
tr=angulo/57.2958;
sint=sin(tr);
cost=cos(tr);
U=meteorologia(:,6);
UZ=U;
kst=meteorologia(:,8);
DTHDZ=0;
hl=meteorologia(:,9);
Ta=meteorologia(:,7);
pwr=puz(kst(:,1));

sp=[200  100;100  200;50 50];
recp=[450  300;560  800;750 350;1050 350];
for i=1:length(sp(:,1))
    for n=1:length(recp(:,1))
    R(i,n)=(recp(n,1)-sp(i,1));
    S(i,n)=recp(n,2)-sp(i,2);
    end
end
for i=1:length(sp(:,1))
    for n=1:length(recp(:,1))
        for mi=1:length(angulo(:,1))
    X(i,n,mi)=S(i,n)*cost(mi,1)+R(i,n)*sint(mi,1);
    Y(i,n,mi)=S(i,n)*sint(mi,1)+R(i,n)*cost(mi,1);
        end
    end
end% end
X=abs(X);
Y=abs(Y);

%Cálculo del factor de boyancia (Fo)

for i=1:length(Ta(:,1))
fo(i,1)=3.12*0.785*vs*(d^2)*((Tg-Ta(i,1))/Tg);
if fo(i,1)>55,
    x0(i,1)=34*exp(0.4*log(fo(i,1)));
else
    x0(i,1)=14*exp(0.625*log(fo(i,1)));
end
end
for i=1:length(fo(:,1))
%Cálculo de delta h (dh)
dh(i,1)= (1.6*exp((log(fo(i,1)))/3)*exp((2*log(3.5*x0(i,1)))/3))/UZ(i,1);

%Cálculo de la altura efectiva (H)
H(i,1)=(h+dh(i,1));
%Cálculo del coeficiente de dispersion horizontal y vertical(G)
end

[Gy,Gz]=dispersion(kst,X);
% concentracion= (Q/(2*3.1416*v*Gz*Gy))*(exp((-0.5)*(y)/Gy)^2)*((exp((-0.5)*((z-H)/Gz)^2))+(exp((-0.5)*(z+H)/Gz)^2));


%________________________________
function [Gy,Gz]=dispersion(kst,X);
[n1,m1,s1]=size(X(:,:,:));

% Calculo los coeficientes de dispersion horizontal(Gy) y vertical (Gz)
for ii=1:length(kst(:,1))
    if kst(ii,1)==1,
        for i1=1:n1
            for i2=1:m1
                for i3=1:s1
                    Gy(i1,i2,i3)= exp(5.357+0.8828*(log(X(i1,i2,i3)/1000))-0.0076*(log(X(i1,i2,i3)/1000))^2);
                    Gz(i1,i2,i3)= exp(6.035+2.1097*(log(X(i1,i2,i3)/1000))+0.277*(log(X(i1,i2,i3)/1000))^2);
                end
            end
        end
    else if kst(ii,1)==2,
            for i1=1:n1
                for i2=1:m1
                    for i3=1:s1
                        Gy(i1,i2,i3)= exp(5.058+0.9024*(log(X(i1,i2,i3)/1000))-0.0096*(log(X(i1,i2,i3)/1000))^2);
                        Gz(i1,i2,i3)= exp(4.694+1.0629*(log(X(i1,i2,i3)/1000))+0.0136*(log(X(i1,i2,i3)/1000))^2);
                    end
                end
            end
        else if kst(ii,1)==3,
                for i1=1:n1
                    for i2=1:m1
                        for i3=1:s1
                            Gy(i1,i2,i3)= exp(4.651+0.9181*(log(X(i1,i2,i3)/1000))-0.0076*(log(X(i1,i2,i3)/1000))^2);
                            Gz(i1,i2,i3)= exp(4.11+0.9201*(log(X(i1,i2,i3)/1000))-0.002*(log(X(i1,i2,i3)/1000))^2);
                        end
                    end
                end
            else if kst(ii,1)==4,
                    for i1=1:n1
                        for i2=1:m1
                            for i3=1:s1
                                Gy(i1,i2,i3)= exp(4.23+0.9222*(log(X(i1,i2,i3)/1000))-0.0087*(log(X(i1,i2,i3)/1000))^2);
                                Gz(i1,i2,i3)= exp(3.414+0.7371*(log(X(i1,i2,i3)/1000))-0.0316*(log(X(i1,i2,i3)/1000))^2);
                            end
                        end
                    end
                else if kst(ii,1)==5,
                        for i1=1:n1
                            for i2=1:m1
                                for i3=1:s1 ,
                                    Gy(i1,i2,i3)= exp(3.922+0.9222*(log(X(i1,i2,i3)/1000))-0.0064*(log(X(i1,i2,i3)/1000))^2);
                                    Gz(i1,i2,i3)= exp(3.057+0.6794*(log(X(i1,i2,i3)/1000))-0.045*(log(X(i1,i2,i3)/1000))^2);
                                end
                            end
                        end
                    else if kst(ii,1)==6,
                            for i1=1:n1
                                for i2=1:m1
                                    for i3=1:s1 ,
                                        Gy(i1,i2,i3)= exp(3.533+0.9181*(log(X(i1,i2,i3)/1000))-0.007*(log(X(i1,i2,i3)/1000))^2);
                                        Gz(i1,i2,i3)= exp(2.621+0.6564*(log(X(i1,i2,i3)/1000))-0.054*(log(X(i1,i2,i3)/1000))^2);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

   