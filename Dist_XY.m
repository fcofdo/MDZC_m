

function [X,Y]=Dist_XY(angulo,XSourc,YSourc,Xmax,N)

angulo(angulo>360) = angulo(angulo>360) - 360;
angulo(angulo<0)   = angulo(angulo(angulo<0)) + 360;
tr=angulo*pi/180;
WDSIN=sin(tr);
WDCOS=cos(tr);

a=0:N:Xmax;
[YRec,XRec] = meshgrid(a);
XRec=10000-XRec;
DD=(XRec-XSourc);
FF=(YRec-YSourc);
% http://en.wikipedia.org/wiki/Rotation_matrix
for f=1:length(WDSIN(:,1));
for I=1:length(DD(:,1))
  for N=1:length(DD(1,:))
      DD1(I,N,f)=DD(I,N)*WDCOS(f)+FF(I,N)*WDSIN(f);
      FF1(I,N,f)=-DD(I,N)*WDSIN(f)+FF(I,N)*WDCOS(f);
     if DD1(I,N,f)<0;
         DD1(I,N,f)=0;
     end
%      if FF1(I,N,f)<0;
%          FF1(I,N,f)=0;
%      end
  end
end
end
 X=DD1;
Y=FF1;
ssa2=zeros(size(X));ssa3=zeros(size(X));
for i1=1:length(X(1,1,:))
    ssa2(:,:,i1)=isinf( X(:,:,i1));
    ssa3(:,:,i1)=isnan( X(:,:,i1));
end

for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if  ssa2(j,jj,i)==1;
                X(j,jj,i)=0;
            end
            if  ssa3(j,jj,i)==1;
                X(j,jj,i)=0;
            end
        end
    end
end



% X = -((XR-XS)*WDSIN + (YR-YS)*WDCOS)
%       Y =   (XR-XS)*WDCOS - (YR-YS)*WDSIN