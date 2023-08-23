%Computation of the volumetric strain in- and around a rectangular and triangular inclusion. 
% This produces Figure 3 in Cornelissen et al. (2023)
%
%Cornelissen, P., Meulenbroek, B.J., Jansen, J.D. (2023). On the derivation
% of closed-form expressions for displacements, strains and stresses inside 
% poroelastic reservoirs. Journal of Geophysical Research - Solid Earth. 
%
% Pavan Cornelissen
% Department of Geoscience and Engineering
% Faculty of Civil Engineering and Geosciences
% Delft University of Technology
% The Netherlands
% E-mail: p.cornelissen@tudelft.nl
%
%Rectangular inclusion
%Geometry
p = -1; %minimum x-coordinate of the rectangle
q = 1; %maximum x-coordinate of the rectangle
r = -0.5; %minimum y-coordinate of the rectangle
s = 0.5; %maximum y-coordinate of the rectangle

n = 200; %number of points in x- and y-direction
[x,y] = meshgrid(linspace(p-1,q+1,n),linspace(r-1,s+1,n)); %points where the strains are evaluated
[Gxx,Gyy,~] = strain_rectangle(x,y,p,q,r,s); %compute the normal strains

%plotting
subplot(1,2,1), surf(x,y,(Gxx+Gyy)/pi), shading interp, view(0,90), cbar = colorbar;
xlim([min(x(:)),max(x(:))]);
ylim([min(y(:)),max(y(:))])
xlabel('$x$ (m)','Interpreter','latex'), ylabel('$y$ (m)','interpreter','latex')
cbar.Label.String = '$(G_{xx}+G_{yy})/\pi$';
cbar.Label.Interpreter = 'latex';
set(gca,'FontSize',12)

%Triangular inclusion
%Geometry
theta = 70/180*pi; %angle of the hypotenuse of the triangle with the x-axis (in radians)
o = 0; %minimum x-coordinate of the triangle
p = 1/tan(theta); %maximum x-coordinate of the triangle
r = 0; %minimum y-coordinate of the triangle
s = 1; %maximum y-coordinate of the triangle

n = 1000; %number of points in x- and y-direction
[x,y] = meshgrid(linspace(o-0.5,p+0.5,n),linspace(r-0.5,s+0.5,n)); %points where the strains are evaluated
[Gxx,Gyy,~] = strain_triangle(x,y,o,p,r,s,theta); %compute the normal strains

%plotting
subplot(1,2,2), surf(x,y,(Gxx+Gyy)/pi), shading interp, view(0,90), cbar = colorbar;
xlim([min(x(:)),max(x(:))]);
ylim([min(y(:)),max(y(:))])
xlabel('$x$ (m)','Interpreter','latex'), ylabel('$y$ (m)','interpreter','latex')
cbar.Label.String = '$(G_{xx}+G_{yy})/\pi$';
cbar.Label.Interpreter = 'latex';
set(gca,'FontSize',12)
% figuresize(32,10,'cm',0);
% print('Output/figure_volumetric_strain','-dpdf','-r300');
% print('Output/figure_volumetric_strain','-dmeta','-r300');