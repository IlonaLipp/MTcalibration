function plotschematic
    %close all

    %figure('Position',[69 1 1852 961])
    %figure

    colors=get(gca,'colororder');

    plotalpha=0.4;

    scale=1;

    % choose max angle where f0 almost zero
    angles=linspace(0,sqrt(-1/scale*log(0.001)));

    % angle where magnitude of curvature maximised
    maxAbsCurvature=sqrt(3/2/scale);

    % first angle where linear Taylor coefficient equals quadratic
    equal12TaylorA=(-1+sqrt(1+2*scale))/2/scale;
    equal12TaylorB=(+1+sqrt(1+2*scale))/2/scale;

    %% Taylor series of saturation around arbitrary point
    f=@(s,a) 1-exp(-a.^2*s);
    m=[0.5;1;1.5];

    offset=+0.3;

    hold on
    shadedrectangle([0,equal12TaylorA],[0,max(m)+offset],[colors(1,:),plotalpha],[1,1,1,plotalpha])
    text(equal12TaylorA/2,max(m)+offset,'quad.','Rotation',0,'HorizontalAlignment','center','VerticalAlignment',"top",'FontSize',11)

    shadedrectangle([equal12TaylorA,max(maxAbsCurvature-equal12TaylorA,0)],[0,max(m)+offset],[1,1,1,plotalpha],[1,1,1,plotalpha])
    text((equal12TaylorA+maxAbsCurvature)/2,max(m)+offset,'transition','Rotation',0,'HorizontalAlignment','center','VerticalAlignment',"top",'FontSize',11)

    shadedrectangle([maxAbsCurvature,max(angles(end)-maxAbsCurvature,0)],[0,max(m)+offset],[1,1,1,plotalpha],[colors(2,:),plotalpha])
    text((maxAbsCurvature+angles(end))/2,max(m)+offset,'saturation','Rotation',0,'HorizontalAlignment','center','VerticalAlignment',"top",'FontSize',11)

    plot(angles,m*f(scale,angles),'k','LineWidth',1)
    %title('Taylor series decomposition of bound pool saturation at each \beta')

    xlabel('\beta_{loc}')
    ylabel('\delta')
    xlim([0,max(angles)])
    ylim([0,max(m)+offset])
    xticks('');
    yticks('');%[0,1]);

    toffset=0.02;
    text(angles(end),m(1)-toffset,'low myelin','HorizontalAlignment','right','VerticalAlignment','top')
    text(angles(end),m(2)-toffset,'medium myelin','HorizontalAlignment','right','VerticalAlignment','top')
    text(angles(end),m(3)-toffset,'high myelin','HorizontalAlignment','right','VerticalAlignment','top')

%exportgraphics(gcf,'schematic.png')
end

function shadedrectangle(x,y,c1,c2)

    c=[1,1,0,0].*reshape(c1(1:3),1,1,3)+[0,0,1,1].*reshape(c2(1:3),1,1,3);

    patch(x(1)+[0,0,x(2),x(2)],y(1)+[0,y(2),y(2),0],c,'LineStyle','none','FaceAlpha',0.6)

end
