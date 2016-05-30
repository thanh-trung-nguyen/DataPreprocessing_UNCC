function draw_cross_section_two_targets(inpfile1,inpfile2,figname,colormaptype)


%draw the cross-section of a data set with two targets:

eval(['load ' inpfile1]);
u1 = u; X1 = X; Y1 = Y; 
eval(['load ' inpfile2]);


hold off;
u = u + u1; 
u(u==0) = 63; u(u< 0) = 30; 

image(u); truesize(size(u)*5);set(gca,'fontsize',20);
if nargin < 4
    colormap(gray);
else
    colormap(colormaptype);
end
    

%imagesc(u+u1); colorbar; 
% plot the true cross-section:
hold on; plot(X,Y,'-k','linewidth',2); hold off;
hold on; plot(X1,Y1,'-k','linewidth',2); hold off;
set(gca,'xtick',Tick_x,'xticklabel',x(Tick_x));
set(gca,'ytick',Tick_y,'yticklabel',y(Tick_y));    
xlabel('x (m)'); ylabel('y (m)');
if strcmp(figname(end-3:end),'.jpg')
    print('-djpeg90',figname);      
elseif strcmp(figname(end-3:end),'.eps')
    print('-depsc2',figname);
end

