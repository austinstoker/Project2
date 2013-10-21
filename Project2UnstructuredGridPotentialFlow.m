%Fluids Project #2 
% Austin Stoker
% Oct 21 2013
close all
clearvars
clc

mycase=3;
DoPlots=false;

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


nodeVolumes = zeros(numCells,1);
if DoPlots
    figure(1)
    hold on;
end


for i=1:numCells
    use4=false;
    n1 = edges(i,1);
    n2 = edges(i,2);
    n3 = edges(i,3);
    n4 = edges(i,4);
    %isInterior = edges(i,5);
    
    x1 = pointCoords(n1,1);
    x2 = pointCoords(n2,1);
    x3 = pointCoords(n3,1);
    
    y1 = pointCoords(n1,2);
    y2 = pointCoords(n2,2);
    y3 = pointCoords(n3,2);
    
    p1 = pointCoords(n1,:);
    p2 = pointCoords(n2,:);
    p3 = pointCoords(n3,:);
    
    nodeVolumes(n3) = nodeVolumes(n3)+1/2*(x1+x2)*(y2-y1);
    
    
    if n4~=0
        x4 = pointCoords(n4,1);
        y4 = pointCoords(n4,2);
        p4 = pointCoords(n4,:);
        use4=true;
        
        nodeVolumes(n4) = nodeVolumes(n4)-1/2*(x1+x2)*(y2-y1);
        if DoPlots
            plot([p1(1) p4(1)],[p1(2) p4(2)])
            plot([p4(1) p2(1)],[p4(2) p2(2)])
            plot([p2(1) p3(1)],[p2(2) p3(2)])
            plot([p3(1) p1(1)],[p3(2) p1(2)])
        end
    else
        use4=false;
        
        nodeVolumes(n1) = nodeVolumes(n1)+1/2*(x1+x2)*(y2-y1);
        nodeVolumes(n2) = nodeVolumes(n2)+1/2*(x1+x2)*(y2-y1);
        if DoPlots
            plot([p1(1) p2(1)],[p1(2) p2(2)])
            plot([p2(1) p3(1)],[p2(2) p3(2)])
            plot([p3(1) p1(1)],[p3(2) p1(2)])
        end
    end
end
sum(nodeVolumes)/3
    
