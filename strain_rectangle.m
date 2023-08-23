%Analytical solution for scaled strains for a rectangular inclusion
%
%Outputs:
% Gxx = scaled horizontal normal strain
% Gyy = scaled vertical normal strain
% Gxy = scaled shear strain

%Inputs
% x = x-coordinate of points where the strains are evaluated
% y = y-coordinate of points where the strains are evaluated
% p = minimum x-coordinate of the rectangle
% q = maximum x-coordinate of the rectangle
% r = minimum y-coordinate of the rectangle
% s = maximum y-coordinate of the rectangle

% Pavan Cornelissen
% Department of Geoscience and Engineering
% Faculty of Civil Engineering and Geosciences
% Delft University of Technology
% The Netherlands
% E-mail: p.cornelissen@tudelft.nl
%

function [Gxx,Gyy,Gxy] = strain_rectangle(x,y,p,q,r,s)
Gxx = atan((y-s)./(x-q)) - atan((y-s)./(x-p)) - atan((y-r)./(x-q)) + atan((y-r)./(x-p));
Gyy = atan((x-q)./(y-s)) - atan((x-p)./(y-s)) - atan((x-q)./(y-r)) + atan((x-p)./(y-r));
Gxy = 0.5*log(((x-q).^2+(y-s).^2).*((x-p).^2+(y-r).^2)./(((x-q).^2+(y-r).^2).*((x-p).^2+(y-s).^2)));
end