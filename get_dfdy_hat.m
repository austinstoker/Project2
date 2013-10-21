function [ dfdy ] = get_dfdy_hat( pointCoords,edges,fun )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%DoPlots=false;
numPoints=size(pointCoords,1);
numCells=size(edges,1);
dfdy = zeros(numPoints,1);

for i=1:numCells
    %use4=false;
    n1 = edges(i,1);
    n2 = edges(i,2);
    n3 = edges(i,3);
    n4 = edges(i,4);
    %isInterior = edges(i,5);
    
    x1 = pointCoords(n1,1);
    x2 = pointCoords(n2,1);
    %x3 = pointCoords(n3,1);
    
    y1 = pointCoords(n1,2);
    y2 = pointCoords(n2,2);
    %y3 = pointCoords(n3,2);
    
    %p1 = pointCoords(n1,:);
    %p2 = pointCoords(n2,:);
    %p3 = pointCoords(n3,:);
    
    dfdy(n3) = dfdy(n3)+1/2*(fun(x1,y1)+fun(x2,y2))*(x2-x1);
    
    
    if n4~=0
        %x4 = pointCoords(n4,1);
        %y4 = pointCoords(n4,2);
        %p4 = pointCoords(n4,:);
        %use4=true;
        
        dfdy(n4) = dfdy(n4)-1/2*(fun(x1,y1)+fun(x2,y2))*(x2-x1);
%         if DoPlots
%             plot([p1(1) p4(1)],[p1(2) p4(2)])
%             plot([p4(1) p2(1)],[p4(2) p2(2)])
%             plot([p2(1) p3(1)],[p2(2) p3(2)])
%             plot([p3(1) p1(1)],[p3(2) p1(2)])
%         end
    else
        %use4=false;
        
        dfdy(n1) = dfdy(n1)+1/2*(fun(x1,y1)+fun(x2,y2))*(x2-x1);
        dfdy(n2) = dfdy(n2)+1/2*(fun(x1,y1)+fun(x2,y2))*(x2-x1);
%         if DoPlots
%             plot([p1(1) p2(1)],[p1(2) p2(2)])
%             plot([p2(1) p3(1)],[p2(2) p3(2)])
%             plot([p3(1) p1(1)],[p3(2) p1(2)])
%         end
    end
end


end

