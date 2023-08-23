%Analytical solution for scaled strains for a triangular inclusion
%
%Outputs:
% Gxx = scaled horizontal normal strain
% Gyy = scaled vertical normal strain
% Gxy = scaled shear strain

%Inputs
% x = x-coordinate of points where the strains are evaluated
% y = y-coordinate of points where the strains are evaluated
% o = minimum x-coordinate of the triangle
% p = maximum x-coordinate of the triangle
% r = minimum y-coordinate of the triangle
% s = maximum y-coordinate of the triangle
% theta = angle of the hypotenuse of the triangle with the x-axis (in radians)

% Pavan Cornelissen
% Department of Geoscience and Engineering
% Faculty of Civil Engineering and Geosciences
% Delft University of Technology
% The Netherlands
% E-mail: p.cornelissen@tudelft.nl
%

function [Gxx,Gyy,Gxy] = strain_triangle(x,y,o,p,r,s,theta)
Gxx = 0.5*log(((y-r).^2+(x-r*cot(theta)).^2)./((y-s).^2+(x-s*cot(theta)).^2))*sin(theta)*cos(theta) ...
    + atan2((r-s).*(x-p),(x-p).^2+(y-r).*(y-s)) ...
    - atan2((s-r).*(y*cot(theta)-x),(x.^2+y.^2+r.*s*csc(theta)^2 - (r+s).*(y+x*cot(theta))))*sin(theta)^2;

Gyy = 0.5*sin(theta)*cos(theta)*log(((x-p).^2+(y-p*tan(theta)).^2)./((x-o).^2+(y-o*tan(theta)).^2)) ...
    - atan2((o-p).*(y-r),(y-r).^2+(x-p).*(x-o)) ...
    - atan2((p-o).*(y-x*tan(theta)),(x.^2+y.^2+o.*p*sec(theta)^2-(p+o).*(x+y*tan(theta))))*cos(theta)^2;

Gxy = 0.5*log(((x-p).^2+(y-s).^2)./((x-p).^2+(y-r).^2)) + 0.5*log(((y-r).^2+(x-r*cot(theta)).^2)./((y-s).^2+(x-s*cot(theta)).^2))*sin(theta)^2 ...
    + atan2((s-r).*(y*cot(theta)-x),(x.^2+y.^2+r.*s*csc(theta)^2-(r+s).*(y+x*cot(theta))))*sin(theta)*cos(theta);
end