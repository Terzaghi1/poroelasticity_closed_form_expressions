% Computation of the strain discontinuity across inclusion boundaries of
% various orientations
%
% This produces Figure 4 in Cornelissen et al. (2023)
% Cornelissen, P., Meulenbroek, B.J., Jansen, J.D. (2023). On the derivation
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
%Constants
mu = 6500; %shear modulus (MPa)
nu = 0.15; %Poisson's ratio (-)
lambda = 2*mu*nu/(1-2*nu); %Lame's first parameter (MPa)
alpha = 0.9; %Biot's coefficient (-)
pf = 1; %incremental pore pressure (MPa)
D = (1-2*nu)*alpha*pf/(2*pi*(1-nu)*mu); %Strain scaling coefficient (-)

%For angles 1 to 89 degrees, we compute the jump in strains across the
%hypotenuse of a triangular inclusion
dip = (1:89)'; %angles for a triangular inclusion

r = 0; %minimum y-coordinate of the triangle
s = 1; %maximum y-coordinate of the triangle
o = 0; %minimum x-coordinate of the triangle

f = 1e-6; %scaling factor to determine how close to the boundary to evaluate the strains

%Initialize output data
Delta_Gxx_Required = zeros(length(dip),1);
Delta_Gyy_Required = zeros(length(dip),1);
Delta_Gxy_Required = zeros(length(dip),1);

Delta_Gxx_CornelissenEtAl = zeros(length(dip),1);
Delta_Gyy_CornelissenEtAl = zeros(length(dip),1);
Delta_Gxy_CornelissenEtAl = zeros(length(dip),1);

Delta_Gxx_WuEtAl = zeros(length(dip),1);
Delta_Gyy_WuEtAl = zeros(length(dip),1);
Delta_Gxy_WuEtAl = zeros(length(dip),1);

for i = 1:length(dip)
    theta = dip(i)/180*pi; %convert degrees to radians
    p = (s-r)/tan(theta); %maximum x-coordinate of the triangle
    u = 1/tan(theta); %x-component of vector parallel to hypotenuse
    v = s-r; %y-component of vector parallel to hypotenuse
    n = [-v, u]; %vector in normal directon of hypotenuse
    n = n/norm(n); %unit vector normal to the hypotenuse
    nx = n(:,1); %x-component of unit normal vector
    ny = n(:,2); %y-component of unit normal vector

    %Required jump in strains to comply with continuity of normal traction    %   
    %First, we compute beta_x and beta_y from Equations 31-32 (Cornelissen
    %et al., 2023)
    beta_x = alpha*pf*((lambda+mu)*nx*ny^2 - ((lambda+2*mu)*ny^2+mu*nx^2)*nx)/(((lambda+2*mu)*nx^2+mu*ny^2)*((lambda+2*mu)*ny^2+mu*nx^2) - ((lambda+mu)*nx*ny)^2);
    beta_y = alpha*pf*((lambda+mu)*nx^2*ny - ((lambda+2*mu)*nx^2+mu*ny^2)*ny)/(((lambda+2*mu)*nx^2+mu*ny^2)*((lambda+2*mu)*ny^2+mu*nx^2) - ((lambda+mu)*nx*ny)^2);
    %Compute the jump in scaled strains by from Equations 27-29 (Cornelissen et al., 2023) 
    Delta_Gxx_Required(i) = 2*beta_x*nx/D;
    Delta_Gyy_Required(i) = 2*beta_y*ny/D;
    Delta_Gxy_Required(i) = (beta_x*ny + beta_y*nx)/D;

    x_in = 0.5./tan(theta) - f*nx; %x-coordinate of points just inside the inclusion
    x_out = 0.5./tan(theta) + f*nx ; %x-coordinate of point just outside the inclusion
    y_in = 0.5*(s-r) - f*ny; %y-coordinate of points just inside the inclusion
    y_out = 0.5*(s-r) + f*ny; %y-coordinate of points just inside the inclusion

    [Gxx_out,Gyy_out,Gxy_out] = strain_triangle(x_out,y_out,o,p,r,s,theta);
    [Gxx_in,Gyy_in,Gxy_in] = strain_triangle(x_in,y_in,o,p,r,s,theta);

    Delta_Gxx_CornelissenEtAl(i) = Gxx_out - Gxx_in;
    Delta_Gyy_CornelissenEtAl(i) = Gyy_out - Gyy_in;
    Delta_Gxy_CornelissenEtAl(i) = Gxy_out - Gxy_in;

    Delta_Gxx_WuEtAl(i) = Gxx_out - (Gxx_in - pi);
    Delta_Gyy_WuEtAl(i) = -Gxx_out - (-Gxx_in + pi);
    Delta_Gxy_WuEtAl(i) = Gxy_out - Gxy_in;
end

%For a horizontal boundary (dip = 0 degrees) and a vertical boundary (dip =
%90 degrees), we could use the short edges of a triangular inclusion or use
%a rectangular inclusion. We opt to use the latter method.

%Square inclusion
p = 0; %minimum x-coordinate of the rectangle
q = 1; %maximum x-coordinate of the rectangle
r = 0; %minimum y-coordinate of the rectangle
s = 1; %maximum y-coordinate of the rectangle

%Horizontal boundary (dip = 0 degrees)
%Required jump in strains to comply with continuity of normal traction
nx = 0;
ny = -1;
%First, we compute beta_x and beta_y from Equations 31-32 (Cornelissen
%et al., 2023)
beta_x = alpha*pf*((lambda+mu)*nx*ny^2 - ((lambda+2*mu)*ny^2+mu*nx^2)*nx)/(((lambda+2*mu)*nx^2+mu*ny^2)*((lambda+2*mu)*ny^2+mu*nx^2) - ((lambda+mu)*nx*ny)^2);
beta_y = alpha*pf*((lambda+mu)*nx^2*ny - ((lambda+2*mu)*nx^2+mu*ny^2)*ny)/(((lambda+2*mu)*nx^2+mu*ny^2)*((lambda+2*mu)*ny^2+mu*nx^2) - ((lambda+mu)*nx*ny)^2);
%Compute the jump in scaled strains by from Equations 27-29 (Cornelissen et al., 2023)
Delta_Gxx_Required = [2*beta_x*nx/D; Delta_Gxx_Required];
Delta_Gyy_Required = [2*beta_y*ny/D; Delta_Gyy_Required];
Delta_Gxy_Required = [(beta_x*ny + beta_y*nx)/D; Delta_Gxy_Required];

x_in = 0.5*(q-p); %x-coordinate of points just inside the inclusion
x_out = x_in; %x-coordinate of point just outside the inclusion
y_in = r+f; %y-coordinate of points just inside the inclusion
y_out = r-f; %y-coordinate of points just inside the inclusion

[Gxx_out,Gyy_out,Gxy_out] = strain_rectangle(x_out,y_out,p,q,r,s);
[Gxx_in,Gyy_in,Gxy_in] = strain_rectangle(x_in,y_in,p,q,r,s);

Delta_Gxx_CornelissenEtAl = [Gxx_out - Gxx_in; Delta_Gxx_CornelissenEtAl];
Delta_Gyy_CornelissenEtAl = [Gyy_out - Gyy_in; Delta_Gyy_CornelissenEtAl];
Delta_Gxy_CornelissenEtAl = [Gxy_out - Gxy_in; Delta_Gxy_CornelissenEtAl];

Delta_Gxx_WuEtAl = [Gxx_out - (Gxx_in - pi); Delta_Gxx_WuEtAl];
Delta_Gyy_WuEtAl = [-Gxx_out - (-Gxx_in + pi); Delta_Gyy_WuEtAl];
Delta_Gxy_WuEtAl = [Gxy_out - Gxy_in; Delta_Gxy_WuEtAl];

%Vertical boundary (dip = 90 degrees)
%Required jump in strains to comply with continuity of normal traction
nx = -1;
ny = 0;
%First, we compute beta_x and beta_y from Equations 31-32 (Cornelissen
%et al., 2023)
beta_x = alpha*pf*((lambda+mu)*nx*ny^2 - ((lambda+2*mu)*ny^2+mu*nx^2)*nx)/(((lambda+2*mu)*nx^2+mu*ny^2)*((lambda+2*mu)*ny^2+mu*nx^2) - ((lambda+mu)*nx*ny)^2);
beta_y = alpha*pf*((lambda+mu)*nx^2*ny - ((lambda+2*mu)*nx^2+mu*ny^2)*ny)/(((lambda+2*mu)*nx^2+mu*ny^2)*((lambda+2*mu)*ny^2+mu*nx^2) - ((lambda+mu)*nx*ny)^2);
%Compute the jump in scaled strains by from Equations 27-29 (Cornelissen et al., 2023)
Delta_Gxx_Required = [Delta_Gxx_Required; 2*beta_x*nx/D];
Delta_Gyy_Required = [Delta_Gyy_Required; 2*beta_y*ny/D];
Delta_Gxy_Required = [Delta_Gxy_Required; (beta_x*ny + beta_y*nx)/D];

x_in = p+f; %x-coordinate of points just inside the inclusion
x_out = p-f; %x-coordinate of point just outside the inclusion
y_in = 0.5*(s-r); %y-coordinate of points just inside the inclusion
y_out = y_in; %y-coordinate of points just inside the inclusion

[Gxx_out,Gyy_out,Gxy_out] = strain_rectangle(x_out,y_out,p,q,r,s);
[Gxx_in,Gyy_in,Gxy_in] = strain_rectangle(x_in,y_in,p,q,r,s);

Delta_Gxx_CornelissenEtAl = [Delta_Gxx_CornelissenEtAl; Gxx_out - Gxx_in];
Delta_Gyy_CornelissenEtAl = [Delta_Gyy_CornelissenEtAl; Gyy_out - Gyy_in];
Delta_Gxy_CornelissenEtAl = [Delta_Gxy_CornelissenEtAl; Gxy_out - Gxy_in];

Delta_Gxx_WuEtAl = [Delta_Gxx_WuEtAl; Gxx_out - (Gxx_in - pi)];
Delta_Gyy_WuEtAl = [Delta_Gyy_WuEtAl; -Gxx_out - (-Gxx_in + pi)];
Delta_Gxy_WuEtAl = [Delta_Gxy_WuEtAl; Gxy_out - Gxy_in];

%plotting
dip = [0; dip; 90];

clf;
lw = 1.2; %linewidth for plotting
ms = 20; %marker size for plotting

plot(dip,Delta_Gxx_Required,'Color',[0 0.4470 0.7410],'LineWidth',lw), hold on;
plot(dip,Delta_Gyy_Required,'Color',[0.8500 0.3250 0.0980],'LineWidth',lw);
plot(dip,Delta_Gxy_Required,'Color',[0.9290 0.6940 0.1250],'LineWidth',lw);

scatter(dip,Delta_Gxx_CornelissenEtAl,ms,[0 0.4470 0.7410]);
scatter(dip,Delta_Gyy_CornelissenEtAl,ms,[0.8500 0.3250 0.0980]);
scatter(dip,Delta_Gxy_CornelissenEtAl,ms,[0.9290 0.6940 0.1250]);

scatter(dip,Delta_Gxx_WuEtAl,ms+4,[0 0.4470 0.7410],'*');
scatter(dip,Delta_Gyy_WuEtAl,ms+4,[0.8500 0.3250 0.0980],'*');
scatter(dip,Delta_Gxy_WuEtAl,ms+4,[0.9290 0.6940 0.1250],'*');

xlabel('$\theta$ ($^\circ$)','interpreter','latex','FontSize',11);
ylabel('Scaled jump in strains (-)','FontSize',11)
legend('$\Delta G_{xx}$ (Eqs. 27-29)','$\Delta G_{yy}$ (Eqs. 27-29)', '$\Delta G_{xy}$ (Eqs. 27-29)', ...     ...
    '$\Delta G_{xx}$ (Eqs. 38-40, 43-45)', '$\Delta G_{yy}$ (Eqs. 38-40, 43-45)', '$\Delta G_{xy}$ (Eqs. 38-40, 43-45)', ...
    '$\Delta G_{xx}$ (Wu et al., 2021)', '$\Delta G_{yy}$ (Wu et al., 2021)', '$\Delta G_{xy}$ (Wu et al., 2021)',...
    'interpreter','latex','NumColumns',3,'FontSize',11,'Location','southoutside');
legend box off;

figuresize(18,12,'cm',0);
print('Output/figure_strain_jump','-dpdf','-r300');
print('Output/figure_strain_jump','-dmeta','-r300');




