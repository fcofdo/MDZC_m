
function [ejeX,ejeY]=ejes(Xmax,Ymax,N)
ejeY=0:N:Ymax;
ejeX=Xmax:-N:0;
ejeX=ejeX';
ejeY=ejeY';