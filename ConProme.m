
function ConProme(Conc,Xmax,Ymax,N,YSourc, XSourc)
 Conc1=zeros(size(Conc(:,:,1)));
 for i=1:length(Conc(1,1,:))
     for j=1:length(Conc(:,1,1))
         for jj=1:length(Conc(1,:,1))
             Conc1(j,jj)=Conc1(j,jj)+Conc(j,jj,i);
         end
     end
 end
 Concpro=zeros(size(Conc(:,:,1)));
 Concpro=Conc1./length(Conc(1,1,:));
     for j=1:length(Concpro(:,1))
        for jj=1:length(Concpro(1,:))
            if  Concpro(j,jj)<.01;
                Concpro(j,jj)=NaN;
            end
        end
    end
[ejeX,ejeY]=ejes(Xmax,Ymax,N);
strColorMap = 'jet_30Intervals.txt';
mCB=IMDC_Read_Matrix(strColorMap,'%f %f %f',' ');
%  h = pcolor(ejeY,ejeX,Concpro(:,:)); shading interp; hold on;
    % contour(ejeY,ejeX,Conc(:,:,n));
%     set(gca,'clim',[20,50]);%colorbar
clist = [5, 10, 20, 30, 40, 50, 60];
[ch] = contourf(ejeY,ejeX,Concpro(:,:), clist );
    grid on;
    axis([0 Xmax 0 Ymax]);
    % axis tight;
    % set(gca,'box','on');
    hold on
    plot(YSourc, XSourc, 'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r' )
    text(YSourc+250, XSourc+250, 'Fuente 1', 'FontSize', 10, 'FontWeight','bold' );
    hold off
    xlabel('East (m)');
    ylabel('North (m)');
    mColorMap=colormap(mCB);
    caxis([0 50]);
hax=gca;
hc=colorbar;
title('Concentraci�n Promedio [microgram/m3]');
set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperUnits', 'centimeters', 'PaperPosition', [0.25 0.25 28 19]);
strOut = ['concentracionmedia','.bmp'];
print( '-dtiff', '-r200', '-opengl', strOut);
close(gcf);
 