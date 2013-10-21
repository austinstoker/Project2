%Fluids Project #2 
% Austin Stoker
% Oct 21 2013
close all
clearvars
clc

mycase=4;
DoPlots=false;

a=1/2;
U=1;
Gam = 4*pi*a*U;
Gam_hat=Gam/(2*pi);

if mycase==1 
    load circle16Edges.mat
    load circle16PointCoords.mat
    pointCoords=circle16PointCoords;
    edges = circle16Edges;
elseif mycase==3
    load circle128Edges.mat
    load circle128PointCoords.mat
    pointCoords=circle128PointCoords;
    edges = circle128Edges;
elseif mycase==4
    load circle512Edges.mat
    load circle512PointCoords.mat
    pointCoords=circle512PointCoords;
    edges = circle512Edges;
end

numPoints=size(pointCoords,1);
numCells=size(edges,1);


xfun = @(x,y) x;
phiHalfBodyfun = @(x,y) U*y-U*(a^2)*y/(x^2+y^2);


phifun = @(x,y) U*(sqrt(x^2+y^2)-((a)^2)/sqrt(x^2+y^2))*y+Gam_hat*log(sqrt(x^2+y^2)/a);

V0 = get_dfdx_hat(pointCoords,edges,xfun);

dfdx = get_dfdx_hat(pointCoords,edges,phifun)./V0;
dfdy = -1*get_dfdy_hat(pointCoords,edges,phifun)./V0;

u=dfdy;
v=-dfdx;

figure(1)
%surf(pointCoords(:,1),pointCoords(:,2),u);
figure(2);
quiver(pointCoords(:,1),pointCoords(:,2),u/10,v/10,'autoscale','off');
    
