% mehhhhh
close all; clear; clc;


% SECTION 1: Configuration
% Variables can be changed by the user ++++++++++++++++++++++++++++++++++++
% In order to find the concentration, copy and paste this function in the command window: 
% C(find(y_ == meting_y), find(x_ == meting_x), find(z_ == meting_z))

x_min = -20;                                                               % minima and maxima of the domain (m)
x_max = 100;
y_min = -50;
y_max = 50;
z_min = 0;
z_max = 30;

step_y = 1;                                                                % step size of y (m)
stability = 1; % set from 1-6                                              % stability parameter (1-6)

measured_x = 63;                                                           % x position of measurement (m)
measured_y = 38;                                                           % y position of measurement (m)
measured_z = 0;                                                            % z position of measurement (m)

stack_x = 34;                                                              % x position of stack (m)
stack_y = 13;                                                              % y position of stack (m)
H = 4;                                                                     % stack height (m)

Q = 1;                                                                     % mass emitted per unit time (g/s)
windspeed = 4.7;                                                           % wind speed (m/s)
wind_direction=90;                                                         % wind_angle in degrees
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% stability of the atmosphere
stability_str={'Very unstable','Moderately unstable','Slightly unstable',...
               'Neutral','Moderately stable','Very stable'};

dy= step_y;                                                                % resolution of the model in x direction
dx=(x_max - x_min)/(y_max - y_min) * step_y;                               % resolution of the model in y direction
dz=(z_max - z_min)/(y_max - y_min) * step_y;                               % resolution of the model in z direction
x_ = (x_min - 0.5*x_max):dx:(x_max + 0.5*x_max);                           % domain, x-grid (m)
y_ = (y_min - 0.5*y_max):dy:(y_max + 0.5*y_max);                           % y-grid (m)
z_ = z_min:dz:z_max;                                                       % z-grid (m)

% Determination of the wind direction
HS_direction = atan2(measured_y - stack_y, measured_x - stack_x) * 180 / pi;

% -------------------------------------------------------------------------

% SECTION 2: Act on the configuration information
% Decide which stability profile to use
stability_str=stability_str{stability};

% decide what kind of run to do, plan view or y-z slice, or time series
[x,y,z]=meshgrid(x_,y_,z_);                                                % x, y and z defined at all positions on the grid

%--------------------------------------------------------------------------

% SECTION 3: Main loop

x1=x-stack_x;                                                                   % shift the coordinates so that stack is centre point
y1=y-stack_y; 
% components of u (wind speed) in x and y directions
wx=windspeed.*sin((wind_direction-180).*pi./180);
wy=windspeed.*cos((wind_direction-180).*pi./180);

% Need angle between point x, y and the wind direction, so use scalar product:
dot_product=wx.*x1+wy.*y1;
% product of magnitude of vectors:
magnitudes=windspeed.*sqrt(x1.^2+y1.^2); 

% angle between wind and point (x,y)
subtended=acos(dot_product./magnitudes);
% distance to point x,y from stack
hypotenuse=sqrt(x1.^2+y1.^2);

% distance along the wind direction to perpendicular line that intersects
% x,y
downwind=cos(subtended).*hypotenuse;

PC = [122.8 0.94470 24.1670 2.5334;90.673 0.93198 18.3330 1.8096; 61.141 0.91465 12.5 1.0857;34.459 0.86974 8.3330 0.72382;24.26 0.83660 6.25 0.54287; 15.209 0.81558 4.1667 0.36191];
a=PC(stability,1); b=PC(stability,2); c=PC(stability,3); d=PC(stability,4);

sig_z=a.*(downwind./1000).^b;
sig_z(find(sig_z(:)>5000))=5000;

theta=0.017453293.*(c-d.*log(downwind./1000));
sig_y=465.11628.*downwind./1000.*tan(theta);

% Now calculate distance cross wind.
crosswind=sin(subtended).*hypotenuse;

ind=find(downwind>0);
C=zeros(size(downwind));

% calculate sigmas based on stability and distance downwind
    
C(ind)=Q./(2.*pi.*windspeed.*sig_y(ind).*sig_z(ind)).*exp(-crosswind(ind).^2./(2.*sig_y(ind).^2)).* (exp(-(z(ind)-H).^2./(2.*sig_z(ind).^2))+ exp(-(z(ind)+H).^2./(2.*sig_z(ind).^2)) );

C_PV = squeeze(C(:, :, find(z_ == measured_z)));                           % concentration of plan view


B = imtranslate(C, [mean(x_)-stack_x, mean(y_)-stack_y, mean(z_)-H]);
A = imrotate3(B, HS_direction, [0 0 1], 'crop');
C_2 = imtranslate(A, [-mean(x_)+stack_x, -mean(y_)+stack_y, -mean(z_)+H]); % translation and rotation of concentration for height slice
C_HS = squeeze(C_2(find(y_ == stack_y), :, :))';                           % concentration of height slice

point = [measured_x - stack_x, measured_y - stack_y, measured_z - H];
rotm = eul2rotm([deg2rad(HS_direction) 0 0]);
measurement = point * rotm;
measurement2 = measurement + [stack_x, stack_y, H];                        % translation and rotation of measured x-value
measured_x2 = measurement2(1);

%--------------------------------------------------------------------------

% SECTION 4: Post process / output
% Plot height slice
figure;
hold on
pcolor(x_,z_,C_HS.*1e6);shading flat       
title({'Gaussian Plume Model Height Slice', stability_str});
xlabel('x (meters)');
ylabel('z (meters)');
caxis manual
caxis([0 max(C_PV,[],'all').*1e6]);
axis([x_min x_max z_min z_max])
h = colorbar;
ylabel(h,'\mu g m^{-3}');
plot(stack_x, H, 'or', 'MarkerFaceColor', 'r')
plot(measured_x2, measured_z, 'ok', 'MarkerFaceColor', 'k')
hold off

% Plot plan view
figure;
hold on
pcolor(x_,y_,C_PV.*1e6);shading flat
title({'Gaussian Plume Model Plan View', stability_str});
xlabel('x (meters)');
ylabel('y (meters)');
caxis manual
caxis([0 max(C_PV,[],'all').*1e6]);
axis([x_min x_max y_min y_max])
h = colorbar;
ylabel(h,'\mu g m^{-3}');
plot(stack_x, stack_y, 'or', 'MarkerFaceColor', 'r')
plot(measured_x, measured_y, 'ok', 'MarkerFaceColor', 'k')
hold off

% Plot 3D model
figure;
isosurface(x_,y_,z_,C.*1e6,(1/2)*sum(C,'all')/numel(C).*1e6);              % 25%
hold on
isosurface(x_,y_,z_,C.*1e6,mean(C,'all').*1e6);                            % 50%
isosurface(x_,y_,z_,C.*1e6,(3/2)*sum(C,'all')/numel(C).*1e6);              % 75%
alpha(0.5);
shading flat;
grid on;
title({'Gaussian Plume Model 3D', stability_str});
xlabel('x (meters)');
ylabel('y (meters)');
zlabel('z (meters)');
caxis manual
caxis([0 max(C_PV,[],'all').*1e6]);
axis([x_min x_max y_min y_max z_min z_max])
h = colorbar;
ylabel(h,'\mu g m^{-3}');
plot3(stack_x, stack_y, H, 'or', 'MarkerFaceColor', 'r')
plot3(measured_x, measured_y, measured_z, 'ok', 'MarkerFaceColor', 'k')
hold off
