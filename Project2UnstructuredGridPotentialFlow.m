%Fluids Project #2 
% Austin Stoker
% Oct 21 2013
close all
clearvars
clc

%program options
GridStart=1;
GridEnd=1;
mycase=3;
DoPlots=true;


%geometery and fluid properties
a=1/2;
U=1;
Gam = (1)*2*pi; %4*pi*a*U;
Gam_hat=Gam/(2*pi);
Pinf = 1000;
ro = 1;


%psiHalfBodyfun = @(x,y) U*y-U*(a^2)*y/(x^2+y^2); % uncomment and replace
%as psifun to analyize a half body

xfun = @(x,y) x;
psifun = @(x,y) U*(sqrt(x^2+y^2)-((a)^2)/sqrt(x^2+y^2))*y + Gam_hat*log(sqrt(x^2+y^2)/a);
syms x y;
uA_hat = diff(psifun(x,y),y);
vA_hat = -diff(psifun(x,y),x);
ufun=inline(char(uA_hat));
vfun=inline(char(vA_hat));
%TODO try not to use inline use @(x,y)

%anaytic vorticity and divergence
dudxA_hat = diff(ufun(x,y),x);
dudxA_fun=inline(char(dudxA_hat));

dvdxA_hat = diff(vfun(x,y),x);
dvdxA_fun=inline(char(dvdxA_hat));

dudyA_hat = diff(ufun(x,y),y);
dudyA_fun=inline(char(dudyA_hat));

dvdyA_hat = diff(vfun(x,y),y);
dvdyA_fun=inline(char(dvdyA_hat));









for k=GridStart:GridEnd
    %load the grid indicated by mycase
    if k==1 
        load circle16Edges.mat
        load circle16PointCoords.mat
        pointCoords=circle16PointCoords;
        edges = circle16Edges;
    elseif k==2
        load circle128Edges.mat
        load circle128PointCoords.mat
        pointCoords=circle128PointCoords;
        edges = circle128Edges;
    elseif k==3
        load circle512Edges.mat
        load circle512PointCoords.mat
        pointCoords=circle512PointCoords;
        edges = circle512Edges;
    end

    numPoints=size(pointCoords,1);
    numCells=size(edges,1);

    %find the volume of each super cell around each node
    V0 = get_dfdx_hat(pointCoords,edges,xfun);

    %preform all the differentiations needed for vorticity and velocity
    %divergence
    dudx = get_dfdx_hat(pointCoords,edges,ufun)./V0;
    dvdx = get_dfdx_hat(pointCoords,edges,vfun)./V0;
    dudy = -1*get_dfdy_hat(pointCoords,edges,ufun)./V0;
    dvdy = -1*get_dfdy_hat(pointCoords,edges,vfun)./V0;

    %assemble vorticity and velocity divergence
    vorticity = dvdx-dudy;
    vel_divergence = dudx+dvdy;

    %populate u and v based on the analytic solution
    u=zeros(numPoints,1);
    v=zeros(numPoints,1);
    for i=1:numPoints
        u(i)=ufun(pointCoords(i,1),pointCoords(i,2));
        v(i)=vfun(pointCoords(i,1),pointCoords(i,2));
    end

    %Calculate the pressure at all nodes
    P = Pinf+(1/2)*ro*(U^2)- (1/2)*ro*(u.^2+v.^2);

    %Plot the velocity field if desired
    if DoPlots
        figure(2);
        quiver(pointCoords(:,1),pointCoords(:,2),u/30,v/30,'autoscale','off');
        xlim([-1 1]);
        ylim([-1 1]);
    end

    %Calculate the forces acting on the cylinder at each edge segment
    j=1;
    for i=1:numCells
        if edges(i,5)==1 % is a cylinder boundary
            n1 = edges(i,1);
            n2 = edges(i,2);
            x1 = pointCoords(n1,1);
            y1 = pointCoords(n1,2);
            x2 = pointCoords(n2,1);
            y2 = pointCoords(n2,2);

            %create a scaled normal vector
            ds = [(y2-y1) -(x2-x1)]; %x and y are swapped because dy creates force in the x direction

            cylP(j,:) = ((P(n1)+P(n2))/2*ds); %avgP*normalvect
            j = j+1;
        end
    end

    %sum the forces around the edge to get lift and drag
    Drag = sum(cylP(:,1))*2;
    Lift = sum(cylP(:,2))*2; %Uhh not sure how to justify the *2 but it works

    %the analytic solutions for Lift and drag on a potential flow rotating
    %cylinder
    DragAnalytic = 0;
    LiftAnalytic = ro*U*Gam;

    %evaluate analytic solutions for vorticity and vel_div at each point
    vorticityA=zeros(numPoints,1);
    vel_divergenceA=zeros(numPoints,1);
    for i=1:numPoints
        x1 = pointCoords(i,1);
        y1 = pointCoords(i,2);
        vorticityA(i) = dvdxA_fun(x1,y1)-dudyA_fun(x1,y1);
        vel_divergenceA(i) = dudxA_fun(x1,y1)+dvdyA_fun(x1,y1);
    end


    % Error analysis

    E_rms_vorticity(k) = sqrt(sum((vorticity-vorticityA).^2)/numPoints);
    E_rms_vel_divergence(k) = sqrt(sum((vel_divergence-vel_divergenceA).^2)/numPoints);

end