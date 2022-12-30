%This ia a final project of Computational Methods for Shallow Water Waves
%We used 2D Linear Shallow water to simulate tsunami at Japan on March 11th in 2011
%The method we use includes finite difference methods, sponger layer and three orders of strong-stability preserving runge-kutta time-steppers

%%%%%%%%%%%%%TSUNAMI SIMULATION%%%%%%%%%%%%%%

clc
clear

%%load intial wave height and water depth
load('bathymetry_MATLAB.mat');
load('IC_MATLAB.mat');

%%set constant parameter
g = 9.81;
geox = X;
geoy = Y;
drawh = bathy;

%%depth of land is set NaN 
for j = 1 : length(geoy)
    for i = 1 : length(geox)
        if bathy(j,i) >= 0
           drawh(j,i) = NaN;
        end
    end
end
for j = 1 : length(geoy)
    for i = 1 : length(geox)
        if bathy(j,i) >= 0
           eta0(j,i) = NaN;
        end
    end
end
for j = 1 : length(geoy)
    for i = 1 : length(geox)
        if bathy(j,i) >= 0
           bathy(j,i) = 0;
        end
    end
end

%%check simulational depth
figure(1)
contour(geox,geoy,abs(bathy),'k')
hold on
surf(geox,geoy,abs(drawh));
view(2)
xlabel('Longitude (\circE)','Fontsize',14)
ylabel('Latitude (\circN)','Fontsize',14)
shading interp
c=colorbar('Location','southoutside');
c.Ticks = [0:1000:10000];c.Label.String = 'Depth (m)';
c.Label.FontSize = 14;
caxis([0 10000])
hold on
contour(geox,geoy,abs(bathy),'k')
hold off

%%check initial wave height
figure(2)
surf(geox,geoy,eta0);
view(2)
xlabel('Longitude (\circE)','Fontsize',14)
ylabel('Latitude (\circN)','Fontsize',14)
shading interp
colormap jet
c=colorbar('Location','southoutside');
c.Ticks = [-8:1:8];
c.Label.String = '\eta (m)';
c.Label.FontSize = 18;
caxis([-8 8])
hold on
contour3(geox,geoy,abs(bathy),'w')
zlim([0 10000])
grid off
hold off

%%K0 = 1/h0*sqrt(3*H/4/h0);
X = X.*90*1000;
Y = Y.*111*1000;

%%set x,y matrix & combine all space x,y to cell
x = linspace(X(1),X(1200),1200);
y = linspace(Y(1),Y(1200),1200);
deltax = x(2)-x(1);
deltay = y(2)-y(1);
xycell = {x;y};

%%set start time,end time,delta time
starttime = 0;
endtime = 3600;

%%define each eta,velocity zero matrix to save space & time's matrix
ng1 = zeros(length(x),1);
U1 = zeros(length(x),1);
V1 = zeros(length(y),1);

%%calculate initial wave height and velocity(1)
R = 0;
for B = 1 : 1
   for j = 1 : length(xycell{2,B})
    for i = 1 : length(xycell{1,B})
    wcell{2+R,1}(j,i) = 0;
    wcell{3+R,1}(j,i) = 0;
    end
   end
   wcell{1+R,1} = interp2(X,Y',eta0,x,y');
   R = R + 4;
end

%%check initial wave height which was interped
figure(3)
contour3(geox,geoy,abs(bathy),'w')
grid off
hold on
surf(geox,geoy,eta0);
view(2)
xlabel('Longitude (\circE)','Fontsize',14)
ylabel('Latitude (\circN)','Fontsize',14)
shading interp
colormap jet
c=colorbar('Location','southoutside');
c.Ticks = [-8:1:8];
c.Label.String = '\eta (m)';
c.Label.FontSize = 14;
caxis([-8 8])
hold on
contour(geox,geoy,abs(bathy),'k')
hold on
plot([x(1)+5*deltax x(1)+5*deltax]/90/1000,[y(1) y(1200)]/111/1000,'-.k','Linewidth',1.5)
hold on
plot([x(1200)-5*deltax x(1200)-5*deltax]/90/1000,[y(1) y(1200)]/111/1000,'-.k','Linewidth',1.5)
hold on
plot([x(1) x(1200)]/90/1000,[y(1)+5*deltax y(1)+5*deltax]/111/1000,'-.k','Linewidth',1.5)
hold on
plot([x(1) x(1200)]/90/1000,[y(1200)-5*deltax y(1200)-5*deltax]/111/1000,'-.k','Linewidth',1.5)
hold off

%%define water depth matrix
hcell = cell(1,length(deltax));
for p = 1 : 1
    hcell{1,p} = interp2(X,Y',bathy,x,y'); 
end
hmax = max(abs(hcell{1,1}),[],'all');

%%courant number limit
for i = 1 : 1
    tcell{1,i} = linspace(starttime,endtime,floor((endtime-starttime)/(0.9*deltax(1,i)/sqrt(g*hmax))+2));
end

%water depth must be positive. negative data need to tranfer to positive
for j = 1 : length(y)
    for i = 1 : length(x)
        if hcell{1,1}(j,i) < 0
           hcell{1,1}(j,i) = -hcell{1,1}(j,i);
        else
           hcell{1,1}(j,i) = 0;
        end
    end
end

%%check water depth which was interped
%{
figure(1)
surf(x,y,hcell{1,1});
shading interp

figure(2)
surf(x,y,wcell{1,1});
shading interp
%}

%%save time interval
savetime = tcell{1,1}(1,2);
%user want the data related time
savecell = wcell;
savecell{4,1} = 0;
LL = 0;

%%simulation begin
for L = 1 : length(deltax)%different dx or dy set
  %%ng,U,V,ng1st,U1st,V1st,ng2nd,U2nd,V2nd,ng3rd,U3rd,V3d, h are temporary variable in this
  %%function at every iteration. Data will save in wcell finally.
  if L ~= 1
      LL = LL+4;
  end
  U = zeros(length(x)+4,length(y)+4);
  V = zeros(length(x)+4,length(y)+4);
  ng = zeros(length(x)+4,length(y)+4);
  h = zeros(length(x)+4,length(y)+4);
  U(3:length(x)+2,3:length(y)+2) = cell2mat(savecell(2+LL,1));
  V(3:length(x)+2,3:length(y)+2) = cell2mat(savecell(3+LL,1));
  ng(3:length(x)+2,3:length(y)+2) = cell2mat(savecell(1+LL,1));
  tt = cell2mat(tcell(1,L));
  h(3:length(x)+2,3:length(y)+2) = cell2mat(hcell(1,L));
  dx = deltax(1,L);
  dy = deltay(1,L);
  U1st = zeros(length(x),length(y));
  V1st = zeros(length(x),length(y));
  ng1st = zeros(length(x),length(y));
  U2nd = zeros(length(x),length(y));
  V2nd = zeros(length(x),length(y));
  ng2nd = zeros(length(x),length(y));
  U3rd = zeros(length(x),length(y));
  V3rd = zeros(length(x),length(y));
  ng3rd = zeros(length(x),length(y));
  ngR = cell2mat(savecell(1+LL,1));
  UR = cell2mat(savecell(2+LL,1));
  VR = cell2mat(savecell(3+LL,1));
  
  savetimecount = 1;
  %Save data setting
  %ex:if want interval with 3 second to save data,
  %but delta t is 1.15 then use 3.3 second(1.15*2) to save data
  %savecell{3,:} will save the really time
  timeindex = round(savetime/(tt(1,2) - tt(1,1)));
  timeiterinti = 1;
  timeiter(1,1) = 1;
  for i = 1 :  length(tt(1,:))
      if timeiterinti < (length(tt(1,:)))
      timeiterinti = timeiterinti + timeindex;
      timeiter(1,i) = timeiterinti;
      end
  end
  for n = 2 : (length(tt(1,:))-1)%each time
    %delta t
    dt = tt(1,n+1) - tt(1,n);
     for E = 1 : 3%SSPRK 3 round
     ng(:,1) = ng (:,5);
      ng(:,2) = ng (:,4);
      U(:,1) = -U (:,5);
      U(:,2) = -U (:,4);
      V(1,:) = -V (5,:);
      V(2,:) = -V (4,:);
      ng(1,:) = ng (5,:);
      ng(2,:) = ng (4,:);
      ng(:,length(x)+2+1) = ng (:,length(x)+2-1);
      ng(:,length(x)+2+2) = ng (:,length(x)+2-2);
      ng(length(y)+2+1,:) = ng(length(y)+2-1,:);
      ng(length(y)+2+2,:) = ng(length(y)+2-2,:);
      U(:,length(x)+2+1) = -U(:,length(x)+2-1);
      U(:,length(x)+2+2) = -U(:,length(x)+2-2);
      U(length(y)+2+1,:) = U(length(y)+2-1,:);
      U(length(y)+2+2,:) = U(length(y)+2-2,:);
      V(length(y)+2+1,:) = -V(length(y)+2-1,:);
      V(length(y)+2+2,:) = -V(length(y)+2-2,:);
      V(:,length(x)+2+1) = V(:,length(x)+2-1);
      V(:,length(x)+2+2) = V(:,length(x)+2-2);
      
      h(:,1) = h(:,5);
      h(:,2) = h(:,4);
      h(1,:) = h (5,:);
      h(2,:) = h (4,:);
      h(:,length(x)+2+1) = h(:,length(x)+2-1);
      h(:,length(x)+2+2) = h(:,length(x)+2-2);
      h(length(y)+2+1,:) = h (length(y)+2-1,:);
      h(length(y)+2+2,:) = h (length(y)+2-2,:);
     for j = 3 : length(y)+2%each y points in space
      for i = 3 : length(x)+2%each x points in space
      
      UM = U(j,i-1);
      UMM = U(j,i-2);
      VM = V(j-1,i);
      VMM = V(j-2,i);
      ngxM = ng(j,i-1);
      ngxMM = ng(j,i-2);
      ngyM = ng(j-1,i);
      ngyMM = ng(j-2,i);
      
      UP =  U(j,i+1);
      UPP = U(j,i+2);
      VP =  V(j+1,i);
      VPP = V(j+2,i);
      ngxP = ng(j,i+1);
      ngxPP = ng(j,i+2);
      ngyP = ng(j+1,i);
      ngyPP = ng(j+2,i);
      
      hxM = h(j,i-1);
      hxMM = h(j,i-2);
      hxP = h(j,i+1);
      hxPP = h(j,i+2);
      hyM = h(j-1,i);
      hyMM = h(j-2,i);
      hyP = h(j+1,i);
      hyPP = h(j+2,i);
      %end
        if E == 1
        %first round
        ng1st(j-2,i-2) = ngR(j-2,i-2) - dt/12/dx*(-UPP*hxPP+8*UP*hxP-8*UM*hxM+UMM*hxMM)-dt/12/dy*(-VPP*hyPP+8*VP*hyP-8*VM*hyM+VMM*hyMM);
        U1st(j-2,i-2) = UR(j-2,i-2) - dt/12/dx*g*(-ngxPP+8*ngxP-8*ngxM+ngxMM);
        V1st(j-2,i-2) = VR(j-2,i-2) - dt/12/dy*g*(-ngyPP+8*ngyP-8*ngyM+ngyMM);
        if j == length(y)+2 && i == length(x)+2
        U(3:(length(x)+2),3:(length(y)+2)) = U1st;
        V(3:length(x)+2,3:length(y)+2) = V1st;
        ng(3:length(x)+2,3:length(y)+2) = ng1st;
        end
        elseif E == 2
        %Second round
        ng2nd(j-2,i-2) = 3/4*ngR(j-2,i-2)+1/4*ng1st(j-2,i-2) - 1/4*dt/12/dx*(-UPP*hxPP+8*UP*hxP-8*UM*hxM+UMM*hxMM)- 1/4*dt/12/dy*(-VPP*hyPP+8*VP*hyP-8*VM*hyM+VMM*hyMM);
        U2nd(j-2,i-2) = 3/4*UR(j-2,i-2)+1/4*U1st(j-2,i-2) - 1/4*dt/12/dx*g*(-ngxPP+8*ngxP-8*ngxM+ngxMM);
        V2nd(j-2,i-2) = 3/4*VR(j-2,i-2)+1/4*V1st(j-2,i-2) - 1/4*dt/12/dy*g*(-ngyPP+8*ngyP-8*ngyM+ngyMM);
        if j == length(y)+2 && i == length(x)+2
        U(3:length(x)+2,3:length(y)+2) = U2nd;
        V(3:length(x)+2,3:length(y)+2) = V2nd;
        ng(3:length(x)+2,3:length(y)+2) = ng2nd;
        end
        elseif E == 3 
        %Third round
        %sponge layer
          Ls = 5*deltax;
          if abs(x(i-2)-X(1200)) <= Ls 
          alphax = 0.5 - 0.5*cos(abs(x(i-2)-X(1200))*pi/Ls);
          else
          alphax = 1;
          end
          if abs(y(j-2)-Y(1200)) <= Ls
          alphay = 0.5 - 0.5*cos(abs(y(j-2)-Y(1200))*pi/Ls);
          elseif abs(y(j-2)-Y(1)) <= Ls
          alphay = 0.5 - 0.5*cos(abs(y(j-2)-Y(1200))*pi/Ls);
          else
          alphay = 1;
          end
        ng3rd(j-2,i-2) = alphax*alphay*(1/3*ngR(j-2,i-2)+2/3*ng2nd(j-2,i-2) - 2/3*dt/12/dx*(-UPP*hxPP+8*UP*hxP-8*UM*hxM+UMM*hxMM)-2/3*dt/12/dy*(-VPP*hyPP+8*VP*hyP-8*VM*hyM+VMM*hyMM));
        U3rd(j-2,i-2) = alphax*alphay*(1/3*UR(j-2,i-2)+2/3*U2nd(j-2,i-2) - 2/3*dt/12/dx*g*(-ngxPP+8*ngxP-8*ngxM+ngxMM));
        V3rd(j-2,i-2) = alphax*alphay*(1/3*VR(j-2,i-2)+2/3*V2nd(j-2,i-2) - 2/3*dt/12/dy*g*(-ngyPP+8*ngyP-8*ngyM+ngyMM));
        if j == length(y)+2 && i == length(x)+2
        ng(3:length(x)+2,3:length(y)+2) = ng3rd;
        U(3:length(x)+2,3:length(y)+2) = U3rd;
        V(3:length(x)+2,3:length(y)+2) = V3rd;
        ngR = ng3rd;
        UR = U3rd;
        VR = V3rd;
        end
        else
        end
       end
      end
    end
        %information of ng3rd, U3rd, and V3rd is next time condition.
        %if intermediate iteration infor doesn't want,then ng(:,now) and
        %will be replaced ng3rd & U3rd & V3rd
        if n ~= (length(tt(1,:))-1)
        ng(3:length(x)+2,3:length(y)+2) = ng3rd;
        U(3:length(x)+2,3:length(y)+2) = U3rd;
        V(3:length(x)+2,3:length(y)+2) = V3rd;
        ngR = ng3rd;
        UR = U3rd;
        VR = V3rd;
        end
        %{
        %if starttime isn't 0,then must not save initial condition
        if starttime ~= 0
           savecell{LL+1,1} = 0;
           savetimecount = 0;
        end
        %}
        %if n equals timeiter, then save data
        %{
        %if n equals timeiter, then save data
        %for p = 1 : length(timeiter)
            if n == timeiter(1,p)
              savetimecount = savetimecount + 1;
              savecell{LL+1,savetimecount} = ng3rd;
              savecell{LL+2,savetimecount} = U3rd;
              savecell{LL+3,savetimecount} = V3rd;
              savecell{LL+4,savetimecount} = tt(1,n);
            end
           end
        end
        %}
  end
  %save end time's data
  savetimecount = savetimecount + 1;
  savecell{LL+1,savetimecount} = ng3rd;
  savecell{LL+2,savetimecount} = U3rd;
  savecell{LL+3,savetimecount} = V3rd;
  savecell{LL+4,savetimecount} = endtime; 
end

%plot coastline
for j = 1 : length(y)
    for i = 1 : length(x)
        if hcell{1,1}(j,i) == 0
           geo(j,i) = 1;
        else 
           geo(j,i) = 0;
        end
    end
end

for K =  1 : length(savecell(1,:))
   for j = 1 : length(y)
    for i = 1 : length(x)
        if hcell{1,1}(j,i) == 0
           savecell{1,K}(j,i) = NaN;
        end
    end
   end
end

%%plot code
figure(3)
for i = 1:20:length(savecell(1,:))
contour(x,y,geo,'k');
colormap default
hold on
[xx, yy] = meshgrid(x, y);
P = surf(xx, yy, savecell{1,i});
zlim([-100 100])
ylim([Y(1) Y(1200)])
xlim([X(1) X(1200)])
xlabel('x(m)','Fontsize',14)
ylabel('y(m)','Fontsize',14)
shading interp
view(2)
c=colorbar('Location','southoutside');
c.Ticks = [-8:1:8];
caxis([-8 8])
c.Label.String = '\eta(m)';
c.Label.FontSize = 14;
colormap jet
hold on
contour(x,y,hcell{1,1},'w')
hold off
drawnow
end
