function [Ex,Ey] = BoxFieldFDPlots(L,W,xmesh,ymesh,Lb,Wb,V0,sig0)
%BoxFieldFDPlots Calculates the potential at each point in space, with the x
%boundaries held at V0 and the y boundaries held at zero. There are also
%boxes of dimensions Lb,Wb at the centre of the top and bottom edges, with
%conductivity sig0. This function also draws plots of the voltage and
%electric field.
%   L: x length, m
%   W: y length, m
%   xmesh,ymesh: number of mesh steps
%   Lb,Wb: dimensions of boxes, m
%   V0: x-boundary potential, V
%   sig0: Conductivity of boxes, Ohm^-1

dx = L/xmesh;
dy = W/ymesh;

n = 0;
nxp = 0;
nxm = 0;
nyp = 0;
nym = 0;
G = spalloc(xmesh*ymesh,xmesh*ymesh,xmesh*ymesh*10);
B = zeros(1,xmesh*ymesh);

x = linspace(0, L, xmesh);
y = linspace(0, W, ymesh);
inBoxes = @(x,y) (x < L/2+Lb & x > L/2-Lb).*(y < Wb | y > W-Wb);
sigma = @(x,y) sig0.*inBoxes(x,y) + ~inBoxes(x,y);

for i = 1:xmesh
   for j = 1:ymesh
      n = j + (i-1)*ymesh;
      
      nxp = j + i*ymesh;
      nxm = j + (i-2)*ymesh;
      nyp = (j+1) + (i-1)*ymesh;
      nym = (j-1) + (i-1)*ymesh;
      
      if (i == 1)
         G(n,n) = (1/dx^2);
         B(n) = (V0/dx^2);
      elseif (i == xmesh)
         G(n,n) = (1/dx^2);
         B(n) = (V0/dx^2);
      elseif (j == 1)
         G(n,n) = (1/dy^2);
      elseif (j == ymesh)
         G(n,n) = (1/dy^2);
      else
         G(n,n) = (-2*(1/dx^2 + 1/dy^2))*sigma(x(i),y(j));
         G(n,nxp) = (1/dx^2);
         G(n,nxm) = (1/dx^2);
         G(n,nyp) = (1/dy^2);
         G(n,nym) = (1/dy^2);
      end
   end
end
V = G\B.';

Vcalc = zeros(xmesh,ymesh);
for i=1:xmesh
    for j = 1:ymesh
        n = j + (i-1)*ymesh;
        Vcalc(i,j) = V(n); 
    end
end
figure(9);
surf(y,x,Vcalc);
title('Voltage Plot');
xlabel('y (m)');
ylabel('x (m)');
zlabel('V (V)');
[Ex,Ey] = gradient(Vcalc);
Ex = Ex/dx;
Ey = Ey/dy;
figure(10);
quiver(y,x,Ex,Ey);
xlabel('y (m)');
ylabel('x (m)');
title('Electric Field Vector Plot');
xlim([0 W]);
ylim([0 L]);

end

