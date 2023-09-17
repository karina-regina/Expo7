%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSIAN PLUME MODEL                                                    %
% PAUL CONNOLLY (UNIVERSITY OF MANCHESTER, 2014)                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not change these variables                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SECTION 0: Definitions (don't modify this section)
close all; clear; clc;


% stability of the atmosphere
stability_str={'Very unstable','Moderately unstable','Slightly unstable',...
               'Neutral','Moderately stable','Very stable'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

stability = 4; % set from 1-6                                              % stability parameter (1-6)

measured_x = 56;                                                           % x position of measurement (m)
measured_y = 38;                                                           % y position of measurement (m)
measured_z = 0;                                                            % z position of measurement (m)

stack_x = 34;                                                              % x position of stack (m)
stack_y = 13;                                                              % y position of stack (m)
H = 4;                                                                     % stack height (m)

Q = 1;                                                                     % mass emitted per unit time (g/s)
windspeed = 4.7;                                                           % wind speed (m/s)

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dy= step_y;                                                                % resolution of the model in x direction
dx=(x_max - x_min)/(y_max - y_min) * step_y;                               % resolution of the model in y direction
dz=(z_max - z_min)/(y_max - y_min) * step_y;                               % resolution of the model in z direction
x_ = (x_min - 0.5*x_max):dx:(x_max + 0.5*x_max);                           % domain, x-grid (m)
y_ = (y_min - 0.5*y_max):dy:(y_max + 0.5*y_max);                           % y-grid (m)
z_ = z_min:dz:z_max;                                                       % z-grid (m)

% Determination of the wind direction
a = atand(abs((measured_y - stack_y) / (measured_x - stack_x)));
if (stack_x <= measured_x & stack_y <= measured_y)
    wind_direction = -90 - a;
    HS_direction = a;
elseif (stack_x > measured_x & stack_y <= measured_y)
    wind_direction = 90 + a;
    HS_direction = -a;
elseif (stack_x > measured_x & stack_y > measured_y)
    wind_direction = 90 - a;
    HS_direction = a;
elseif (stack_x <= measured_x & stack_y > measured_y)
    wind_direction = -90 + a;
    HS_direction = -a;
end

% -------------------------------------------------------------------------

% SECTION 2: Act on the configuration information
% Decide which stability profile to use
stability_str=stability_str{stability};

% decide what kind of run to do, plan view or y-z slice, or time series
[x,y,z]=meshgrid(x_,y_,z_);                                                % x, y and z defined at all positions on the grid

%--------------------------------------------------------------------------

% SECTION 3: Main loop
C=gauss_func(Q,windspeed,wind_direction,x,y,z, ...                         % concentration of 3D model
             stack_x,stack_y,H,stability);


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

% Functions ---------------------------------------------------------------
function C=gauss_func(Q,u,dir1,x,y,z,xs,ys,H,STABILITY)
u1=u;
x1=x-xs;                                                                   % shift the coordinates so that stack is centre point
y1=y-ys; 

% components of u (wind speed) in x and y directions
wx=u1.*sin((dir1-180).*pi./180);
wy=u1.*cos((dir1-180).*pi./180);

% Need angle between point x, y and the wind direction, so use scalar product:
dot_product=wx.*x1+wy.*y1;
% product of magnitude of vectors:
magnitudes=u1.*sqrt(x1.^2+y1.^2); 

% angle between wind and point (x,y)
subtended=acos(dot_product./magnitudes);
% distance to point x,y from stack
hypotenuse=sqrt(x1.^2+y1.^2);

% distance along the wind direction to perpendicular line that intersects
% x,y
downwind=cos(subtended).*hypotenuse;

% Now calculate distance cross wind.
crosswind=sin(subtended).*hypotenuse;

ind=find(downwind>0);
C=zeros(size(downwind));

% calculate sigmas based on stability and distance downwind
[sig_y,sig_z]=calc_sigmas(STABILITY,downwind);
    
C(ind)=Q./(2.*pi.*u1.*sig_y(ind).*sig_z(ind)) ...
    .*exp(-crosswind(ind).^2./(2.*sig_y(ind).^2)).* ...
    (exp(-(z(ind)-H).^2./(2.*sig_z(ind).^2))+ ...
    exp(-(z(ind)+H).^2./(2.*sig_z(ind).^2)) );
end


function [sig_y,sig_z]=calc_sigmas(CATEGORY,x)

a=zeros(size(x));
b=zeros(size(x));
c=zeros(size(x));
d=zeros(size(x));

switch CATEGORY
    case 1 % very unstable
        % vertical
        a(:)=122.800;
        b(:)=0.94470;

        % cross wind
        c(:)=24.1670;
        d(:)=2.5334;
    case 2 % moderately unstable
        % vertical
        a(:)=90.673;
        b(:)=0.93198;
          
        % cross wind
        c(:)=18.3330;
        d(:)=1.8096;
    case 3 % slightly unstable
        % vertical
        a(:)=61.141;
        b(:)=0.91465;
        
        % cross wind
        c(:)=12.5;
        d(:)=1.0857;
    case 4 % neutral
        % vertical
        a(:)=34.459;
        b(:)=0.86974;
             
        % cross wind
        c(:)=8.3330;
        d(:)=0.72382;
    case 5 % moderately stable
        % vertical
        a(:)=24.26;
        b(:)=0.83660;
        
        % cross wind
        c(:)=6.25;
        d(:)=0.54287;
    case 6 % very stable
        % vertical
        a(:)=15.209;
        b(:)=0.81558;
        
        % cross wind
        c(:)=4.1667;
        d(:)=0.36191;
    otherwise
        return;
end

sig_z=a.*(x./1000).^b;
sig_z(find(sig_z(:)>5000))=5000;

theta=0.017453293.*(c-d.*log(x./1000));
sig_y=465.11628.*x./1000.*tan(theta);
end

