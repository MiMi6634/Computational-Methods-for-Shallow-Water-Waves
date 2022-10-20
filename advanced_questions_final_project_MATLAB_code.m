%%%%%%%%%%%%%maxheight map%%%%%%%%%%
clc
clear
load('savecell_dxdy375')
load('bathymetry_MATLAB.mat');
for i = 1:length(savecell_dxdy375(1,:))
    for D = 1:length(savecell_dxdy375{1,i}(1,:))
        for W = 1:length(savecell_dxdy375{1,i}(:,1))
            etacell{W,D}(i) = savecell_dxdy375{1,i}(W,D); 
        end
    end
end
for D = 1:length(savecell_dxdy375{1,i}(1,:))
        for W = 1:length(savecell_dxdy375{1,i}(:,1))
            maxeta(W,D) = max(etacell{W,D}(:)); 
        end
end
x3 = linspace(X(1),X(1200),1200);
y3 = linspace(Y(1),Y(1200),1480);
 
h = interp2(X,Y',bathy,x3,y3');
 
for T = 1 : length(h(:,1))
    for F = 1:length(h(1,:))
        if h(T,F)>=0
           maxeta(T,F) = NaN;
        end
    end
end
 
geo = zeros(length(h(:,1)),length((h(1,:))));
%{
for T = 1 : length(h(:,1))
    for F = 1:length(h(1,:))
        if h(T,F)>0
           geo(T,F) = NaN;
       
        else 
           geo(T,F) = abs(h(T,F));
        end
    end
end
%} 
load('bathymetry_MATLAB.mat');
figure(1)
[xx, yy] = meshgrid(x3, y3);
P = surf(xx, yy,maxeta );
ylim([Y(1) 40.2])
xlim([140 144.5])
xlabel('Longitude (\circE)','Fontsize',14)
ylabel('Latitude (\circN)','Fontsize',14)
set(gca, 'FontName', 'Helvetica')
shading interp
view(0,90)
c=colorbar('Location','southoutside');
c.Ticks = [-8:1:8];
caxis([-8 8])
c.Label.String = '\eta(m)';
c.Label.FontSize = 14;
colormap jet
box on
grid off
hold on
contour3(x3,y3(1:1332),abs(h(1:1332,:)),'Color','w');
ylim([Y(1) 40])
xlim([X(1) 144.5])
zlim([0 10000])
hold off
%%%%%%expand simulation
clc
clear all
load('bathymetry_MATLAB.mat')
load('IC_MATLAB.mat');
lon = ncread('gebco_2020_n50.0_s30.0_w130.0_e155.0.nc','lon');
lat = ncread('gebco_2020_n50.0_s30.0_w130.0_e155.0.nc','lat');
info = ncread('gebco_2020_n50.0_s30.0_w130.0_e155.0.nc','elevation');
bathy=fliplr(rot90(double(info),3));%%%%adujst the direction of the map
y=(lon-lon(1))*111000;%unit m
x=(lat-lat(1))*90000;%unit m
delta_x=x(2)-x(1);
delta_y=delta_x;
p1=1:1:6000;
p2=1:1:4800;
yi=1:2:length(y);%%%%%%%selected the range
xi=1:2:length(x);%%%%%%%selected the range
y1=interp1(p1,y,yi,'linear');
x1=interp1(p2,x,xi,'linear');
delta_x=x1(2)-x1(1)
delta_y=delta_x
%%%%%%%%%%%%%%%%deal the initial poisition of tsunami
for i=1:length(x)
    if lon(i)>=140
        k=lon(i)
        break
    end
end
for i1=1:length(x)
    if lon(i1)>=145
        k1=lon(i1)
        break
    end
end
for j=1:length(y)
    if lat(j)>=35.5
        k2=lat(j)
        break
    end
end
for j1=1:length(y)
    if lat(j1)>=40.5
        k3=lat(j1)
        break
    end
end
eta_zero=zeros(length(x),length(y));
eta_zero(j:j1-1,i:i1-1)=eta0;
eta0=eta_zero;
[x2 y2]= meshgrid(y,x);
[x2i y2i]= meshgrid(y1,x1);
eta01=interp2(x2,y2,eta0,x2i,y2i,'cubic');%%%interpt the initial condition
bathy1=-interp2(x2,y2,bathy,x2i,y2i,'cubic');%%%interpt the bathy
g=9.81;
%%%%%%%%create zero matrix and add the initail condition
alpha=zeros(1,length(y1)+4);
alpha1=zeros(length(x1)+4,1);
alpha3=zeros(length(x1)+4,1);
alpha2=zeros(1,length(y1)+4);
U1=zeros(length(x1)+4,length(y1)+4);
V1=zeros(length(x1)+4,length(y1)+4);
ita2=zeros(length(x1)+4,length(y1)+4);
ita2(3:length(x1)+2,3:length(y1)+2)=eta01;
h=zeros(length(x1)+4,length(y1)+4);
h(3:length(x1)+2,3:length(y1)+2)=bathy1;
y=y1;
x=x1;
h(1,:)=h(5,:);
h(2,:)=h(4,:);
h(length(x)+4,:)=h(length(x),:);
h(length(x)+3,:)=h(length(x)+1,:);
h(:,1)=h(:,5);
h(:,2)=h(:,4);
h(:,length(y)+4)=h(:,length(y));
h(:,length(y)+3)=h(:,length(y)+1);
    for i=1:length(y)
        for j=1:length(x)
            if h(j,i)<0
                h(j,i)=0;
            end
        end
    end
    delta_x=x(2)-x(1);
    delta_t=min(delta_x,delta_y)/sqrt(g*9000)*0.9;
    t=0:delta_t:3600;
    C={};
    ls=5*delta_x
    ls1=5*delta_y
     n=1;
     %%%%%%%%sponge layer
for i=1:length(y1)  %%%left boundary
    if y1(i)<=ls;
        alpha(1,i+2)=0.5-0.5*cos((y1(i))*pi/ls);
    else
        alpha(1,i+2)=1;
    end
end
for j=1:length(x1)
    if x1(j)<=ls1 ; %%% up boundary
        alpha1(j+2,1)=0.5+0.5*cos((x1(j)-ls1)*pi/ls1);
    else
        alpha1(j+2,1)=1;
    end
end
for i=1:length(y1)%%%right boundary
    if y1(i)>=y1(end)-ls ;
        alpha2(1,i+2)=0.5-0.5*cos((y1(end)-y1(i))*pi/ls);
    else
        alpha2(1,i+2)=1;
    end
end
for j=1:length(x1)%%%down boundary
    if x1(j)>=x1(end)-ls1;
        alpha3(j+2,1)=0.5-0.5*cos((x1(end)-x1(j))*pi/ls1);
    else
        alpha3(j+2,1)=1;
    end
end
x=y1;
y=x1;
%%%%%%%%%%simulation begin
for k=1:length(t)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%first around above%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    %%%%%%%%%%%%%%%%%%%%%%second round above%%%%%%%%%%%%%%%%%%%%%%%
               for i=3:length(x)+2
        for j=3:length(y)+2
            ita2(j,i)=1/3*ita2(j,i)+2/3*ita_double_star(j,i)-2/3*(delta_t/(12*delta_x))*(-U1_double_star(j,i+2)*h(j,i+2)+8*U1_double_star(j,i+1)*h(j,i+1)-8*U1_double_star(j,i-1)*h(j,i-1)+U1_double_star(j,i-2)*h(j,i-2))-2/3*(delta_t/(12*delta_y))*(-V1_double_star(j+2,i)*h(j+2,i)+8*V1_double_star(j+1,i)*h(j+1,i)-8*V1_double_star(j-1,i)*h(j-1,i)+V1_double_star(j-2,i)*h(j-2,i));
            V1(j,i)=1/3*V1(j,i)+2/3*V1_double_star(j,i)-2/3*(delta_t/(12*delta_x))*g*(-ita_double_star(j+2,i)+8*ita_double_star(j+1,i)-8*ita_double_star(j-1,i)+ita_double_star(j-2,i));
            U1(j,i)=1/3*U1(j,i)+2/3*U1_double_star(j,i)-2/3*(delta_t/(12*delta_x))*g*(-ita_double_star(j,i+2)+8*ita_double_star(j,i+1)-8*ita_double_star(j,i-1)+ita_double_star(j,i-2));
        end
               end
               %%%%%%%%%%disappear wave by sponge layer
   for i=3:length(x)+2
                   for j=3:length(y)+2
                       ita2(j,i)=ita2(j,i)*alpha(1,i)*alpha1(j,1)*alpha2(1,i)*alpha3(j,1);
                       V1(j,i)=V1(j,i)*alpha1(j,1)*alpha3(j,1);
                       U1(j,i)=U1(j,i)*alpha(1,i)*alpha2(1,i);
                   end
   end
        ita3=ita2;
%%%%%%%%%%select data because the data is too bid that saving data has
     %%%%%%%%%%problem and increase the computing Efficiency
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
               surf(x,y,ita4(3:length(y)+2,3:length(x)+2));
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
