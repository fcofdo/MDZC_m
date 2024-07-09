
function [Conc]=Refina(X,Conc)
for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if  Conc(j,jj,i)<.1;
                Conc(j,jj,i)=NaN;
            end
        end
    end
end