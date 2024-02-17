%% For a cylinder of radii r and R, and height L, plot the undeformed and deform configurations
% when the deformation gradient is given by a complex function as shown in
% the class.
% Changes: f(r) to be f(x_3)

clear all, close all, clc

I = eye(3); %Identity matrix

n0 = [1;1;0]; % Use this vector to compute change of length for each point in this direction.
n1 = [0;1;1]; % Use this vector to compute change of angle at each point between n0 and n1

xmin = -2;
xmax = 2;

ymin = -2;
ymax = 2;

zmin = 0;
zmax = 1.5;

%f(r)
radial_function = 'sinusoidal';%'quadratic'; %'linear'; %'sinusoidal'; %

%phi(x_3)
twist_angle = 'exponential'; %'sinusoidal'; %'constant'; %

a = 0.25; %Inner radius
b = 0.75; %Outer radius
L = 1; %Height
lambda = 1.25; %compression/stretching value

r0 = a; %minimum radius of the inner layer
R0 = b; %minimum radius of the outer layer

%Inner surface
R=a*ones(50,1);
[X,Y,Z] = cylinder(R);

%Outer surface
r=b*ones(50,1);
[x,y,z] = cylinder(r);

%Plot both surfaces
figure(1)
subplot(1,2,1)
surf(X,Y,Z)
hold on
surf(x,y,z)
axis([xmin xmax ymin ymax zmin zmax])
title('Undeformed cylinder')
set(gca, 'FontSize',18,'FontWeight','bold')

xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
disp('Press a key !') % Press a key here. You can see the message 'Paused: Press any key' in
% the lower left corner of the
% pause;

%% Compute the new points

%Create a vector that goes from 0 to 1 to make several drawings of the transformation

step = 1;

for k=1:length(step)
for i=1:length(X)
    for j=1:length(X(1,:))

        %Get the point
        xp = [X(i,j); Y(i,j); Z(i,j)];
        r = sqrt(xp(1)^2+xp(2)^2);

        %Compute the functions f(r) and phi(r)
        switch radial_function
            case 'linear'
                f = R0+2*xp(3);
            case 'quadratic'
                f = R0+xp(3)*xp(3);
            case 'sinusoidal'
                f= r0+r0*sin(pi*xp(3)/L);
        end

        switch twist_angle
            case 'constant'
                phi = 30*pi/180; %30 degrees -> rad
            case 'sinusoidal'
                phi = 30*sin(pi*xp(3)/L)*pi/180; %30 degrees -> rad
            case 'quadratic'
                phi = 30*xp(3)^2*pi/180; %
            case 'exponential'
                phi = 30*exp(2*xp(3))*pi/180;
        end
        
        F = [f*cos(phi) -f*sin(phi) 0
             f*sin(phi)  f*cos(phi) 0
             0           0           lambda];

        Fm1=inv(F);

        yp = F*xp;
        y1(i,j) = yp(1);
        y2(i,j) = yp(2);
        y3(i,j) = yp(3);
        EL(i,j,1:3,1:3) = F'*F-I; % Lagrange strain tensor
        el(i,j,1:3,1:3) = I-Fm1'*Fm1; % Eulerian strain tensor
        dl(i,j) = norm(F*n0); % Change in length in n0-direction.
        cos_theta(i,j) = (F*n0)'*(F*n1)/(norm(F*n0)*norm(F*n1));
    end
end

for i=1:length(x)
    for j=1:length(x(1,:))
        
        %Get the point
        xp = [x(i,j); y(i,j); z(i,j)];
        r = sqrt(xp(1)^2+xp(2)^2);
        
        %Compute the functions f(r) and phi(r)
        switch radial_function
            case 'linear'
                f = R0+2*xp(3);
            case 'quadratic'
                f = R0+xp(3)*xp(3);
            case 'sinusoidal'
                f = R0+R0*sin(pi*xp(3)/L);
        end

        switch twist_angle
            case 'constant'
                phi = 30*pi/180; %30 degrees -> rad
            case 'sinusoidal'
                phi = 30*sin(pi*xp(3)/L)*pi/180; %30 degrees -> rad
            case 'quadratic'
                phi = 30*xp(3)^2*pi/180; %
            case 'exponential'
                phi = 30*exp(2*xp(3))*pi/180;
        end

        F = [f*cos(phi) -f*sin(phi) 0
             f*sin(phi)  f*cos(phi) 0
             0           0           lambda];

        Fm1=inv(F);

        yp = F*xp;
        ym1(i,j) = yp(1);
        ym2(i,j) = yp(2);
        ym3(i,j) = yp(3);
        ELm(i,j,1:3,1:3) = F'*F-I; % Lagrange strain tensor
        elm(i,j,1:3,1:3) = I-Fm1'*Fm1; % Eulerian strain tensor
        DL(i,j) = norm(F*n0); % Change in length in n0-direction.
        cos_THETA(i,j) = (F*n0)'*(F*n1)/(norm(F*n0)*norm(F*n1));

    end
end

subplot(1,2,2)
hold off
%Plot strain measure in color
%surf(y1,y2,y3,eIm(:,:,1,2))

%Plot change of length along n0
%surf(y1,y2,y3,dL)

%Plot origin and after deformation
surf(y1,y2,y3,ELm(:,:,1,2))
surf(y1,y2,y3,elm(:,:,1,2))

hold on
surf(ym1,ym2,ym3,elm(:,:,1,1))



hold on
%Plot strain measure in color
%surf(ym1,ym2,ym3,eIm(:,:,1,2))

%Plot change of length along n0
%surf(ym1,ym2,ym3,DL)
%surf(y1,y2,y3,ELm(:,:,1,1))

%Plot change of angle between n0 and n1
%surf(ym1,ym2,ym3,cos_THETA)

axis([xmin xmax ymin ymax zmin zmax])
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')

set(gca, 'FontSize', 18, 'FontWeight', 'bold')
title('Deformed cylinder')
drawnow
pause(0.1)

end
