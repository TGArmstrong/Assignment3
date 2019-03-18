function [] = EFieldElectrons(Etype,PlotJ)
%   E: Type of electric field (see below)
%   PlotJ: If non-zero, generate a plot of current density

clf;
cla;
xlength = 200E-9; % Length of region in the x-direction, m
ylength = 100E-9; % Length of region in the y-direction, m
N = 10000; % Total number of electrons
px = 0; % x position of electrons, m
py = 0; % y position of electrons, m
vx = 0; % x-velocity of electrons, m/s
vy = 0; % y-velocity of electrons, m/s
dx = 1E-9;
dy = 1E-9;

m = 0.26*(9.11E-31); % Electron mass, kg
T = 300; % Temperature, K
q = 1.602E-19; % Electric field, C
dt = 1E-15; % Time step, s
reflecty = zeros(N,1); % For reflecting boundary
reflectx = zeros(N,1);
scatter = zeros(N,1); % For scattering
Pscat = 1 - exp(-dt/0.2E-12); % Scattering probability
Time = 1000; % Total number of time steps
maxTraj = 20; % Maximum number of trajectories to plot

boxes = 0; % Turn boxes on/off
Lb = 20e-9;
Wb = 10e-9;

px = -(xlength/2) + xlength*rand(N,1);
py = -(ylength/2) + ylength*rand(N,1);
vx = sqrt(1.380E-23*T/m)*randn(N,1);
vy = sqrt(1.380E-23*T/m)*randn(N,1);

% Type of electric field
Evx = zeros(round(xlength/dx),round(ylength/dy));
Evy = zeros(round(xlength/dx),round(ylength/dy));
if (Etype == 0)
    % Constant E-field
    Vx = 0.1; % Voltage driving electric field, in V
    Vy = -1;
    Evx(:,:) = Vx*q/m/xlength*dt; 
    Evy(:,:) = Vy*q/m/ylength*dt;
elseif (Etype == 1)
    % Electric field calculated with FD solver
    [Evx, Evy] = SaddleFieldFD(xlength,ylength,round(xlength/dx),round(ylength/dy),1);
    Evx = Evx*q/m*dt;
    Evy = Evy*q/m*dt;
elseif (Etype == 2)
    boxes = 1;
    [Evx, Evy] = BoxFieldFD(xlength,ylength,round(xlength/dx),round(ylength/dy),2*Lb,2*Wb,1,1e-2);
    Evx = Evx*q/m*dt;
    Evy = Evy*q/m*dt;
end

if (boxes == 1) % Ensure that no electrons begin inside the barriers
    while (sum(px(abs(px) < Lb & abs(py) > ylength/2-Wb) + py(abs(py) > ylength/2-Wb & abs(px) < Lb)) ~= 0)
        px(abs(px) < Lb & abs(py) > ylength/2-Wb) = -(xlength/2) + xlength*rand(1,1);
        py(abs(py) > ylength/2-Wb & abs(px) < Lb) = -(xlength/2) + ylength*rand(1,1);
    end
end
colours = rand(min(maxTraj, N), 3);
xold = px(1:min(maxTraj,N));
yold = py(1:min(maxTraj,N));
wrapped = zeros(min(maxTraj,N),1);
figure(1+2*Etype);
for t = 1:Time
   
   xind = 1+floor((px+xlength/2)./dx);
   yind = 1+floor((py+ylength/2)./dy);
   xold = px(1:min(maxTraj,N));
   yold = py(1:min(maxTraj,N));
   % Update velocities with E-field
   vx = vx + Evx(sub2ind(size(Evx),xind,yind));
   vy = vy + Evy(sub2ind(size(Evy),xind,yind));
   % Update positions
   px = px + vx.*dt;
   py = py + vy.*dt;
   
   % Wrapping boundary
   wrapped = (abs(px(1:min(maxTraj,N))) > xlength/2);
   px(px < -xlength/2) = px(px < -xlength/2) + xlength;
   px(px > xlength/2) = px(px > xlength/2) - xlength;
   
   % Reflecting boundary and boxes
   if (boxes == 0)
       reflecty = -1*(abs(py) > ylength/2);
       vy = vy.*(2*reflecty+1);
   end
   if (boxes == 1)
       reflecty = -1*((abs(py) > ylength/2) | (abs(py) > ylength/2-Wb & abs(px) < Lb));
       reflectx = -1*(abs(px) < Lb & abs(py) > ylength/2-Wb);
       vy = vy.*(2*reflecty+1);
       vx = vx.*(2*reflectx+1);
   end
   py = py - vy.*reflecty.*dt;
   px = px - vx.*reflectx.*dt;

   % Scattering
   scatter = rand(N,1) < Pscat;
   vx = vx + (sqrt(1.380E-23*T/m)*randn(N,1) - vx).*scatter;
   vy = vy + (sqrt(1.380E-23*T/m)*randn(N,1) - vy).*scatter;
   
   if (boxes == 1) % Draw the boxes
       plot([-Lb Lb], [ylength/2-Wb ylength/2-Wb], 'k');
       hold on;
       plot([-Lb Lb], [-ylength/2+Wb -ylength/2+Wb], 'k');
       hold on;
       plot([-Lb -Lb], [ylength/2 ylength/2-Wb], 'k');
       hold on;
       plot([-Lb -Lb], [-ylength/2 -ylength/2+Wb], 'k');
       hold on;
       plot([Lb Lb], [ylength/2 ylength/2-Wb], 'k');
       hold on;
       plot([Lb Lb], [-ylength/2 -ylength/2+Wb], 'k');
       hold on;
   end
   
   for i = 1:min(maxTraj,N)
    if (wrapped(i) == 0)
       plot([xold(i) px(i)], [yold(i) py(i)], 'color', colours(i,:));
    end
   end
   xlim([-xlength/2 xlength/2]);
   ylim([-ylength/2 ylength/2]);
   hold on;
   %pause(0.0001);
end
Tav = 0.5*m/(1.380E-23)*sum(vx.^2 + vy.^2)/N;
title('Electron trajectories');
xlabel('x (m)');
ylabel('y (m)');

figure(2+2*Etype);
xlim([-xlength/2 xlength/2]);
ylim([-ylength/2 ylength/2]);
phist = histogram2(px, py, [20 10], 'displaystyle', 'tile');
title('Electron Density Map');
xlabel('x (m)');
ylabel('y (m)');
colorbar;

figure(3+2*Etype);
Tmap = zeros(200,100); 
for i = 1:xlength*1E9
    for j = 1:ylength*1E9
        thiscell = px > (i-1)*1E-9-xlength/2 & px < i*1E-9-xlength/2 & py > (j-1)*1E-9-ylength/2 & py < j*1E-9-ylength/2;
        if (sum(thiscell) <= 0)
            Tmap(i,j) = 0;
        else
            Tmap(i,j) = 0.5*m/(1.380E-23)*sum(vx(thiscell).^2 + vy(thiscell).^2)./sum(thiscell);
        end
    end
end
X = linspace(-xlength/2, xlength/2, xlength*1E9);
Y = linspace(-ylength/2, ylength/2, ylength*1E9);
surfc(Y, X, Tmap);
title('Temperature map');
zlabel('Temperature (K)');
xlabel('x (m)');
ylabel('y (m)');

if PlotJ ~= 0
   figure(4+2*Etype);
   Jx = 1e19*q*vx;
   Jy = 1e19*q*vy;
   Jhist = histogram2(Jx, Jy, [50 50], 'displaystyle', 'tile');
   title('Plot of Current Density');
   xlabel('J_x (A/m)');
   ylabel('J_y (A/m)');
   colorbar;
end