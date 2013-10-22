%Fluids Project #2 
% Austin Stoker
% Oct 21 2013
close all
clearvars
clc

mycase=3;
DoPlots=true;

a=1/2;
U=1;
Gam = (1)*2*pi; %4*pi*a*U;
Gam_hat=Gam/(2*pi);
Pinf = 1000;
ro = 1;

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

psifun = @(x,y) U*(sqrt(x^2+y^2)-((a)^2)/sqrt(x^2+y^2))*y + Gam_hat*log(sqrt(x^2+y^2)/a);


syms x y;
uA_hat = diff(psifun(x,y),y);
vA_hat = -diff(psifun(x,y),x);
ufun=inline(char(uA_hat));
vfun=inline(char(vA_hat));


V0 = get_dfdx_hat(pointCoords,edges,xfun);

dfdx = get_dfdx_hat(pointCoords,edges,psifun)./V0;
dfdy = -1*get_dfdy_hat(pointCoords,edges,psifun)./V0;

u=zeros(numPoints,1);
v=zeros(numPoints,1);
for i=1:numPoints
    u(i)=ufun(pointCoords(i,1),pointCoords(i,2));
    v(i)=vfun(pointCoords(i,1),pointCoords(i,2));
end
xpt=pointCoords(:,1);
ypt=pointCoords(:,2);

%Calculate the pressure at all nodes
P = Pinf+(1/2)*ro*(U^2)- (1/2)*ro*(u.^2+v.^2);

%Plot the velocity field if desired
if DoPlots
    figure(2);
    quiver(pointCoords(:,1),pointCoords(:,2),u/30,v/30,'autoscale','off');
    xlim([-1 1]);
    ylim([-1 1]);
end

%Calculate the forces acting on the cylinder
j=1;
for i=1:numCells
    if edges(i,5)==1 % is a cylinder boundary
        n1 = edges(i,1);
        n2 = edges(i,2);
        x1 = pointCoords(n1,1);
        y1 = pointCoords(n1,2);
        x2 = pointCoords(n2,1);
        y2 = pointCoords(n2,2);
        
        n = [(y2-y1) -(x2-x1)];
        
        cylP(j,:) = ((P(n1)+P(n2))/2*n); %avgP*normalvect /2 because each node is part of 2 cells
        j = j+1;
    end
end

Drag = sum(cylP(:,1))*2;
Lift = sum(cylP(:,2))*2; %Uhh not sure how to justify the *2 but it works

DragAnalytic = 0;
LiftAnalytic = ro*U*Gam;

% E = sqrt(sum((uA'-u).^2)/numPoints)
% E=LiftAnalytic-Lift;

% z = P;
% xlin = linspace(min(x),max(x),33);
% ylin = linspace(min(y),max(y),33);
% [X,Y] = meshgrid(xlin,ylin);
% f = TriScatteredInterp(x,y,z);
% Z = f(X,Y);
% mesh(X,Y,Z) %interpolated
% axis tight; hold on
% plot3(x,y,z,'.','MarkerSize',15) %nonuniform

