
function [Conc]=ConcSourPoint(X,Q,Uadj,Y,SZout1,SYout1,z,H,kst,hl)
% encuentra los sigmas

Conc=zeros(size(X));
for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if X(j,jj,i)<=0;
                Conc(j,jj,i)=0;
            else
                if kst(i,1)>=5
                    Conc(j,jj,i)= (Q*1000000/(2*pi()*Uadj(i)*SZout1(j,jj,i)...
                        *SYout1(j,jj,i)))*(exp((-0.5)*(Y(j,jj,i)/...
                        SYout1(j,jj,i))^2))*(exp((-0.5)*((z(j,jj)-H(i))...
                        /SZout1(j,jj,i))^2)+exp((-0.5)*((z(j,jj)+H(i))/SZout1(j,jj,i))^2));
                end
                if kst(i,1)<=4
                    if SZout1(j,jj,i)> 1.6
                        if hl(i)<H(i)  || hl(i)<z(j,jj)
                            Conc(j,jj,i)=0;
                        end
                        Conc(j,jj,i)= (Q*1000000/(((2*pi())^0.5)*Uadj(i)*...
                            SYout1(j,jj,i)*hl(i)))*(exp((-0.5)*(Y(j,jj,i)/...
                            SYout1(j,jj,i))^2));
                    else
                        g(1,1)=0;
                        for ii1=-4:4
                                g31(1,1)=exp(-0.5*(((z(j,jj)-H(i)+2*ii1*hl(i))^2)/((SZout1(j,jj,i))^2)))+...
                                exp(-0.5*(((z(j,jj)+H(i)+2*ii1*hl(i))^2)/((SZout1(j,jj,i))^2)));
                            g(1,1)=g31+g(1,1);
                        end
                        Conc(j,jj,i)=(Q*1000000/(2*pi()*Uadj(i)*SZout1(j,jj,i)*...
                         SYout1(j,jj,i)))*(exp((-0.5)*(Y(j,jj,i)/SYout1(j,jj,i))^2))*g(1,1);
                    end
                end
            end
        end
    end
end


ssa=zeros(size(X));ssa1=zeros(size(X));
for i1=1:length(X(1,1,:))
    ssa(:,:,i1)=isinf( Conc(:,:,i1));
    ssa1(:,:,i1)=isnan( Conc(:,:,i1));
end

for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if  ssa(j,jj,i)==1;
                Conc(j,jj,i)=0;
            end
            if  ssa1(j,jj,i)==1;
                Conc(j,jj,i)=0;
            end
        end
    end
end

% for i=1:length(X(1,1,:))
%     for j=1:length(X(:,1,1))
%         for jj=1:length(X(1,:,1))
%             if  Conc(j,jj,i)<.1;
%                 Conc(j,jj,i)=NaN;
%             end
%         end
%     end
% end
