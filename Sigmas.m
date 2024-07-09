
function [SZout1,SYout1,H]=Sigmas(X,Ta,Ts,D,Vs,Uadj,hs,kst,OpDis)

SZout1=zeros(size(X));
SYout1=zeros(size(X));
H=zeros(length(Uadj(:,1)),1);
for i1=1:length(Uadj);
 H(i1,1)=Fplume(Ta(i1),Ts,D,Vs,Uadj(i1),hs,kst(i1));
end
for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if X(j,jj,i)<=.1;
                SZout1(j,jj,i)=0;
                SYout1(j,jj,i)=0;
            else
                SZout1(j,jj,i)=SIGZ(X(j,jj,i),OpDis,kst(i));
                SYout1(j,jj,i)=SIGY(X(j,jj,i),OpDis,kst(i));                             
            end
        end
    end
end