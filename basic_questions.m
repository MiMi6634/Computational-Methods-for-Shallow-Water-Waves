%%%%%%%%%%%%%TSUNAMI SIMULATION
clc;
clear all;
load('bathymetry_MATLAB.mat')
load('IC_MATLAB.mat');
y=(X-X(1))*111000;%unit m
x=(Y-Y(1))*90000;%unit m
delta_x=x(2)-x(1);
delta_y=y(2)-y(1);
p=1:1:1200;
yi=1:1:1200;
xi=1:1:1200;
y1=interp1(p,y,yi,'linear');%%%%interpt
x1=interp1(p,x,xi,'linear');
[x2 y2]= meshgrid(x,y);
x1i=x(1):1*delta_x:x(end);
y1i=y(1):1*delta_y:y(end);
[x2i y2i]= meshgrid(x1i,y1i);
eta01=interp2(x2,y2,eta0,x2i,y2i,'cubic');%%%interpt the initial condition 
bathy1=-interp2(x2,y2,bathy,x2i,y2i,'cubic');%%%interpt the bathy
g=9.81;
%%%%%%%%create zero matrix and add the initail condition
alpha=zeros(1,length(x1)+4);
alpha1=zeros(length(y1)+4,1);
alpha3=zeros(length(y1)+4,1);
alpha2=zeros(1,length(x1)+4);
U1=zeros(length(y1)+4,length(x1)+4);
V1=zeros(length(x1)+4,length(x1)+4);
ita2=zeros(length(x1)+4,length(y1)+4);
ita2(3:length(x1)+2,3:length(y1)+2)=eta01;
h=zeros(length(x1)+4,length(x1)+4);
h(3:length(x1)+2,3:length(y1)+2)=bathy1;
y=y1;
x=x1;
h(1,:)=h(5,:);
h(2,:)=h(4,:);
h(length(y)+4,:)=h(length(y),:);
h(length(y)+3,:)=h(length(y)+1,:);
h(:,1)=h(:,5);
h(:,2)=h(:,4);
h(:,length(x)+4)=h(:,length(x));
h(:,length(x)+3)=h(:,length(x)+1);
%%%%%%%%sponge layer
    for i=1:length(h)
        for j=1:length(h)
            if h(j,i)<0
                h(j,i)=0;
            end
        end
    end
    delta_x=x(2)-x(1);
    delta_y=y(2)-y(1);
    delta_t=min(delta_x,delta_y)/sqrt(g*9000)*0.9;
    t=0:delta_t:3600;
    C={};
    ls=5*delta_x;
    ls1=5*delta_y
    n=1;
for i=1:length(x1)  %%%left boundary
    if x1(i)<=ls;
        alpha(1,i+2)=0.5-0.5*cos((x1(i))*pi/ls);
    else
        alpha(1,i+2)=1;
    end
end
for j=1:length(y1)
    if y1(j)<=ls1 ; %%% up boundary
        alpha1(j+2,1)=0.5+0.5*cos((y1(j)-ls1)*pi/ls1);
    else
        alpha1(j+2,1)=1;
    end
end
for i=1:length(x1)%%%right boundary
    if x1(i)>=x1(end)-ls ;
        alpha2(1,i+2)=0.5-0.5*cos((x1(end)-x1(i))*pi/ls);
    else
        alpha2(1,i+2)=1;
    end
end
for j=1:length(y1)%%%down boundary
    if y1(j)>=y1(end)-ls1;
        alpha3(j+2,1)=0.5-0.5*cos((y1(end)-y1(j))*pi/ls1);
    else
        alpha3(j+2,1)=1;
    end
end
%%%%%%%%%%simulation begin
for k=1:length(t)
    %%%%%%%% mirror condition
    ita2(1,:)=ita2(5,:);
    ita2(2,:)=ita2(4,:);
    ita2(length(y)+4,:)=ita2(length(y),:);
    ita2(length(y)+3,:)=ita2(length(y)+1,:);
    ita2(:,1)=ita2(:,5);
    ita2(:,2)=ita2(:,4);
    ita2(:,length(x)+4)=ita2(:,length(x));
    ita2(:,length(x)+3)=ita2(:,length(x)+1);
    U1(1,:)=U1(5,:);
    U1(2,:)=U1(4,:);
    U1(length(y)+4,:)=U1(length(y),:);
    U1(length(y)+3,:)=U1(length(y)+1,:);
    U1(:,1)=-U1(:,5);
    U1(:,2)=-U1(:,4);
    U1(:,length(x)+4)=-U1(:,length(x));
    U1(:,length(x)+3)=-U1(:,length(x)+1);
    V1(1,:)=-V1(5,:);
    V1(2,:)=-V1(4,:);
    V1(length(y)+4,:)=-V1(length(y),:);
    V1(length(y)+3,:)=-V1(length(y)+1,:);
    V1(:,1)=V1(:,5);
    V1(:,2)=V1(:,4);
    V1(:,length(x)+4)=V1(:,length(x));
    V1(:,length(x)+3)=V1(:,length(x)+1);
    %%%%%%%%%%R-K method first round
           for i=3:length(x)+2
        for j=3:length(y)+2
            ita_star(j,i)=ita2(j,i)-(delta_t/(12*delta_x))*(-U1(j,i+2)*h(j,i+2)+8*U1(j,i+1)*h(j,i+1)-8*U1(j,i-1)*h(j,i-1)+U1(j,i-2)*h(j,i-2))-(delta_t/(12*delta_y))*(-V1(j+2,i)*h(j+2,i)+8*V1(j+1,i)*h(j+1,i)-8*V1(j-1,i)*h(j-1,i)+V1(j-2,i)*h(j-2,i));
            U1_star(j,i)=U1(j,i)-(delta_t/(12*delta_x))*g*(-ita2(j,i+2)+8*ita2(j,i+1)-8*ita2(j,i-1)+ita2(j,i-2));
            V1_star(j,i)=V1(j,i)-(delta_t/(12*delta_y))*g*(-ita2(j+2,i)+8*ita2(j+1,i)-8*ita2(j-1,i)+ita2(j-2,i));
        end
           end
        ita_star(1,:)=ita_star(5,:);
    ita_star(2,:)=ita_star(4,:);
    ita_star(length(y)+4,:)=ita_star(length(y),:);
    ita_star(length(y)+3,:)=ita_star(length(y)+1,:);
    ita_star(:,1)=ita_star(:,5);
    ita_star(:,2)=ita_star(:,4);
    ita_star(:,length(x)+4)=ita_star(:,length(x));
    ita_star(:,length(x)+3)=ita_star(:,length(x)+1);
    U1_star(1,:)=U1_star(5,:);
    U1_star(2,:)=U1_star(4,:);
    U1_star(length(y)+4,:)=U1_star(length(y),:);
    U1_star(length(y)+3,:)=U1_star(length(y)+1,:);
    U1_star(:,1)=-U1_star(:,5);
    U1_star(:,2)=-U1_star(:,4);
    U1_star(:,length(x)+4)=-U1_star(:,length(x));
    U1_star(:,length(x)+3)=-U1_star(:,length(x)+1);
    V1_star(1,:)=-V1_star(5,:);
    V1_star(2,:)=-V1_star(4,:);
    V1_star(length(y)+4,:)=-V1_star(length(y),:);
    V1_star(length(y)+3,:)=-V1_star(length(y)+1,:);
    V1_star(:,1)=V1_star(:,5);
    V1_star(:,2)=V1_star(:,4);
    V1_star(:,length(x)+4)=V1_star(:,length(x));
    V1_star(:,length(x)+3)=V1_star(:,length(x)+1);
%%%%%%%%%%R-K method first round
%%%%%%%%%%R-K method second round
           for i=3:length(x)+2
        for j=3:length(y)+2
            ita_double_star(j,i)=3/4*ita2(j,i)+1/4*ita_star(j,i)-1/4*(delta_t/(12*delta_x))*(-U1_star(j,i+2)*h(j,i+2)+8*U1_star(j,i+1)*h(j,i+1)-8*U1_star(j,i-1)*h(j,i-1)+U1_star(j,i-2)*h(j,i-2))-1/4*(delta_t/(12*delta_y))*(-V1_star(j+2,i)*h(j+2,i)+8*V1_star(j+1,i)*h(j+1,i)-8*V1_star(j-1,i)*h(j-1,i)+V1_star(j-2,i)*h(j-2,i));
            V1_double_star(j,i)=3/4*V1(j,i)+1/4*V1_star(j,i)-1/4*(delta_t/(12*delta_x))*g*(-ita_star(j+2,i)+8*ita_star(j+1,i)-8*ita_star(j-1,i)+ita_star(j-2,i));
            U1_double_star(j,i)=3/4*U1(j,i)+1/4*U1_star(j,i)-1/4*(delta_t/(12*delta_y))*g*(-ita_star(j,i+2)+8*ita_star(j,i+1)-8*ita_star(j,i-1)+ita_star(j,i-2));
        end
           end
        ita_double_star(1,:)=ita_double_star(5,:);
    ita_double_star(2,:)=ita_double_star(4,:);
    ita_double_star(length(y)+4,:)=ita_double_star(length(y),:);
    ita_double_star(length(y)+3,:)=ita_double_star(length(y)+1,:);
    ita_double_star(:,1)=ita_double_star(:,5);
    ita_double_star(:,2)=ita_double_star(:,4);
    ita_double_star(:,length(x)+4)=ita_double_star(:,length(x));
    ita_double_star(:,length(x)+3)=ita_double_star(:,length(x)+1);
    U1_double_star(1,:)=U1_double_star(5,:);
    U1_double_star(2,:)=U1_double_star(4,:);
    U1_double_star(length(y)+4,:)=U1_double_star(length(y),:);
    U1_double_star(length(y)+3,:)=U1_double_star(length(y)+1,:);
    U1_double_star(:,1)=-U1_double_star(:,5);
    U1_double_star(:,2)=-U1_double_star(:,4);
    U1_double_star(:,length(x)+4)=-U1_double_star(:,length(x));
    U1_double_star(:,length(x)+3)=-U1_double_star(:,length(x)+1);
    V1_double_star(1,:)=-V1_double_star(5,:);
    V1_double_star(2,:)=-V1_double_star(4,:);
    V1_double_star(length(y)+4,:)=-V1_double_star(length(y),:);
    V1_double_star(length(y)+3,:)=-V1_double_star(length(y)+1,:);
    V1_double_star(:,1)=V1_double_star(:,5);
    V1_double_star(:,2)=V1_double_star(:,4);
    V1_double_star(:,length(x)+4)=V1_double_star(:,length(x));
    V1_double_star(:,length(x)+3)=V1_double_star(:,length(x)+1);
    %%%%%%%%%%R-K method second round
    %%%%%%%%%%R-K method third round
               for i=3:length(x)+2
        for j=3:length(y)+2
            ita2(j,i)=1/3*ita2(j,i)+2/3*ita_double_star(j,i)-2/3*(delta_t/(12*delta_x))*(-U1_double_star(j,i+2)*h(j,i+2)+8*U1_double_star(j,i+1)*h(j,i+1)-8*U1_double_star(j,i-1)*h(j,i-1)+U1_double_star(j,i-2)*h(j,i-2))-2/3*(delta_t/(12*delta_y))*(-V1_double_star(j+2,i)*h(j+2,i)+8*V1_double_star(j+1,i)*h(j+1,i)-8*V1_double_star(j-1,i)*h(j-1,i)+V1_double_star(j-2,i)*h(j-2,i));
            V1(j,i)=1/3*V1(j,i)+2/3*V1_double_star(j,i)-2/3*(delta_t/(12*delta_x))*g*(-ita_double_star(j+2,i)+8*ita_double_star(j+1,i)-8*ita_double_star(j-1,i)+ita_double_star(j-2,i));
            U1(j,i)=1/3*U1(j,i)+2/3*U1_double_star(j,i)-2/3*(delta_t/(12*delta_x))*g*(-ita_double_star(j,i+2)+8*ita_double_star(j,i+1)-8*ita_double_star(j,i-1)+ita_double_star(j,i-2));
        end
               end
  %%%%%%%%%%R-K method third round
  %%%%%%%%%%disappear wave by sponge layer
   for i=3:length(x)+2
                   for j=3:length(y)+2
                       ita2(j,i)=ita2(j,i)*alpha(1,i)*alpha1(j,1)*alpha2(1,i)*alpha3(j,1);
                       V1(j,i)=V1(j,i)*alpha1(j,1)*alpha3(j,1);
                       U1(j,i)=U1(j,i)*alpha(1,i)*alpha2(1,i);
                   end
   end
    %%%%%%%%%%disappear wave by sponge layer
     %%%%%%%%%%select data because the data is too bid that saving data has
     %%%%%%%%%%problem and increase the computing Efficiency
        ita3=ita2;
        if k==20*n
        for i=3:length(x)+2
        for j=3:length(y)+2
            if h(j,i)==0
                ita3(j,i)=NaN;
            end
        end 
        end
              C(n)= {ita3};
              n=n+1;
              t1(n)=t(k);
        end
               end
%%%%%%%%%plot tsunami process
for n=1:1:length(C)
               ita4=C{1,n};
               surf(x,y,ita4(3:length(x)+2,3:length(y)+2));
               colormap jet;
               colorbar ;
               shading interp;
               c.Label.String = '\eta(m)'
               xlabel('x(m)')
               ylabel('y(m)')
               zlabel('\eta(m)')
               title(['t=',num2str(t1(n)),'s'])
               axis([x(1) x(end) y(1) y(end) -8.6 8.6])
               view(0,90)
               caxis([-8.6 8.6])
               pause(0.0001)
               
end
%%save as a gif
clc 
clear all
load('1200s.mat');
load('IC_MATLAB.mat');
filename='homework.gif'
for i=1:length(y)
    if Y(i)>=40
        k=Y(i)
        break
    end
end
for j=1:length(X)
    if X(j)>=144.5
        k=X(i)
        break
    end
end
for n=1:1:length(C); %%%save as a gif
               ita4=C{1,n};
               surf(X,Y,ita4(3:length(x)+2,3:length(y)+2));
               box on
               grid off
               colormap jet ;
               colorbar
               c=colorbar
               c.Location ='southoutside'
               c.Label.String = '\eta(m)'
               caxis([-8.6 8.6]);
               shading interp;
               view(0,90)
               hold on
               contour3(X,Y,h(3:length(x)+2,3:length(x)+2),1000:1000:8000,'w');
               xlabel('Longitude(¢XE)')
               ylabel('Latitude(¢XN)')
               title(['t=2400s'])
               axis([X(1) X(j) Y(1) Y(i)]);
               hold off
               pause(0.01)    
               F=getframe(gcf);
               im=frame2im(F);
               [I,map]=rgb2ind(im,256);
               k=n-0;
              if k==1;
              imwrite(I,map,filename,'gif','Loopcount',inf, 'DelayTime',0.02);
              else
              imwrite(I,map,filename,'gif','WriteMode','append', 'DelayTime',0.02);
              end
end
convergence test
load('bathymetry_MATLAB.mat');
load('savecell_dxdy1127')
load('savecell_dxdy751')
load('savecell_dxdy375')
load('savecell_dxdy281')
 
x1 = linspace(X(1),X(1200),400);
y1 = linspace(Y(1),Y(1200),493);
x2 = linspace(X(1),X(1200),600);
y2 = linspace(Y(1),Y(1200),739);
x3 = linspace(X(1),X(1200),1200);
y3 = linspace(Y(1),Y(1200),1480);
x4 = linspace(X(1),X(1200),1600);
y4 = linspace(Y(1),Y(1200),1973);
deltax1 = (x1(2)-x1(1))*90*1000;
deltay1 = (y1(2)-y1(1))*111*1000;
deltax2 =(x2(2)-x2(1))*90*1000;
deltay2 = (y2(2)-y2(1))*111*1000;
deltax3 = (x3(2)-x3(1))*90*1000;
deltay3 = (y3(2)-y3(1))*111*1000;
deltax4 = (x4(2)-x4(1))*90*1000;
deltay4 = (y4(2)-y4(1))*111*1000;
 
figure(2)
plot(x1,savecell_dxdy1127{1,length(savecell_dxdy1127(1,:))}(247,:),'Linewidth',1.5);%247
hold on
plot(x2,savecell_dxdy751{1,length(savecell_dxdy751(1,:))}(370,:),'-g','Linewidth',1);%370
hold on
plot(x3,savecell_dxdy375{1,length(savecell_dxdy375(1,:))}(740,:),'k','Linewidth',1.5);%740
hold on
plot(x4,savecell_dxdy281{1,2}(987,:),'r','Linewidth',1);%987
hold off
xlim([X(220) X(350)])
title('At t=3600s','Fontsize',14)
h = legend('\Deltax=\Deltay=1127m','\Deltax=\Deltay=751m','\Deltax=\Deltay=375m','\Deltax=\Deltay=281m')
h.FontSize=14;
xlabel('Longitude (\circE)','Fontsize',14)
ylabel('\eta (m)','Fontsize',14)
%%%%%%%%%%compare real data and simulate data
clc
clear all
load 1200.mat
load('IC_MATLAB.mat');
%%%%%%%locking the range
for i=1:length(y)
    if Y(i)>=40
        k=Y(i)
        break
    end
end
for j=1:length(X)
    if X(j)>=144.5
        k=X(i)
        break
    end
end
%%%%%%%add poistion of station
xp=[38.17873 39.21195 38.80402 39.59104 36.92341 38.44083 38.58083];
yp=[141.69138 142.10587 141.90346 142.19832 141.19092 142.51006 142.59586] ;
for n=1:length(xp)
[~,i1]=min(abs(Y(:)-xp(n)));
[~,j1]=min(abs(X(:)-yp(n)));
a(n)=i1
b(n)=j1
end
for n=1;
               ita4=C{1,n};
               surf(X,Y,ita4(3:length(x)+2,3:length(y)+2));
               hold on
               plot3(X(b(1)),Y(a(1)),5,'r*')
               plot3(X(b(2)),Y(a(2)),5,'r*')
               plot3(X(b(3)),Y(a(3)),5,'r*')
               plot3(X(b(4)),Y(a(4)),5,'r*')
               plot3(X(b(5)),Y(a(5)),5,'r*')
               plot3(X(b(6)),Y(a(6)),5,'r*')
               plot3(X(b(7)),Y(a(7)),5,'r*')
               box on
               grid off
               colormap jet ;
               surf(X,Y,ita4(3:length(x)+2,3:length(y)+2));
               caxis([-8.6 8.6]);
               shading interp;
               view(0,90)
               contour3(X,Y,h(3:length(x)+2,3:length(x)+2),1000:1000:8000,'w');
               xlabel('Longitude(¢XE)')
               ylabel('Latitude(¢XN)')
               axis([X(1) X(j) Y(1) Y(i)]);
               hold off
end
