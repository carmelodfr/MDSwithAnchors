function [ Coord ] = srf( Coord )
%srf is a function that shift rotate and flip coordinates Coord using the
%last three nodes in the Matrix

% Copyright (c) 2020
% Author: Carmelo Di Franco 
% Email: carmelodfr@gmail.com
% This code is licensed under MIT license (see LICENSE.txt for details)

 

%shift 
Coord(:,1) = Coord(:,1) - Coord(end,1);
Coord(:,2) = Coord(:,2) - Coord(end,2);
%rotate
angle = -atan2(Coord(end-1,2),Coord(end-1,1));
Rx = [cos(angle)  sin(angle) ;
          -sin(angle) cos(angle) ];
Coord = Coord*Rx;
%flip
if (Coord(end-2,2) < 0) 
    Coord(:,2) = -  Coord(:,2);
end
end

