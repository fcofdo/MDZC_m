
function  grabamovie(X,Conc,Xmax,Ymax,N,YSourc, XSourc)
for i=1:length(X(1,1,:))
    for j=1:length(X(:,1,1))
        for jj=1:length(X(1,:,1))
            if  Conc(j,jj,i)<10;
                Conc(j,jj,i)=NaN;
            end
        end
    end
end
[ejeX,ejeY]=ejes(Xmax,Ymax,N);
nr_fr =length(X(1,1,:));
% Initialize matrix using 'moviein'
frames = moviein(nr_fr);

strColorMap = 'jet_30Intervals.txt';
mCB=IMDC_Read_Matrix(strColorMap,'%f %f %f',' ');
% strColorMap = 'jet';
% Generate frames with any plotting function.
% We use a cosine with variable frequency.
figure(1)
filename = 'gauss.gif';
for n = 1: nr_fr
    mesh(ejeY,ejeX,Conc(:,:,n));
%     h = pcolor(ejeY,ejeX,Conc(:,:,n)); shading interp; hold on;
    % contour(ejeY,ejeX,Conc(:,:,n));
%     set(gca,'clim',[20,50]);%colorbar
    grid on;
    axis([0 Xmax 0 Ymax]);
    % axis tight;
    % set(gca,'box','on');
    hold on
    plot(YSourc, XSourc, 'ro', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r' )
    text(YSourc, XSourc, 'Fuente 1', 'FontSize', 10, 'FontWeight','bold' );
    hold off
    xlabel('East (m)');
    ylabel('North (m)');
    mColorMap=colormap(mCB);
caxis([10 50]);
hax=gca;
hc=colorbar;

    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if n == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end