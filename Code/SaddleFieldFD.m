function [Ex,Ey] = SaddleFieldFD(L,W,xmesh,ymesh,V0)
%SaddleFieldFD Calculates the potential at each point in space, with the x
%boundaries held at V0 and the y boundaries held at zero
%   L: x length, m
%   W: y length, m
%   xmesh,ymesh: number of mesh steps
%   V0: x-boundary potential, V

dx = L/xmesh;
dy = W/ymesh;

n = 0;
nxp = 0;
nxm = 0;
nyp = 0;
nym = 0;
G = spalloc(xmesh*ymesh,xmesh*ymesh,xmesh*ymesh*10);
B = zeros(1,xmesh*ymesh);

for i = 1:xmesh
   for j = 1:ymesh
      n = j + (i-1)*ymesh;
      
      nxp = j + i*ymesh;
      nxm = j + (i-2)*ymesh;
      nyp = (j+1) + (i-1)*ymesh;
      nym = (j-1) + (i-1)*ymesh;
      
      if (i == 1)
         G(n,n) = 1/dx^2;
         B(n) = V0/dx^2;
      elseif (i == xmesh)
         G(n,n) = 1/dx^2;
         B(n) = V0/dx^2;
      elseif (j == 1)
         G(n,n) = 1/dy^2;
      elseif (j == ymesh)
         G(n,n) = 1/dy^2;
      else
         G(n,n) = -2*(1/dx^2 + 1/dy^2);
         G(n,nxp) = 1/dx^2;
         G(n,nxm) = 1/dx^2;
         G(n,nyp) = 1/dy^2;
         G(n,nym) = 1/dy^2;
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
[Ex,Ey] = gradient(Vcalc);
Ex = Ex/dx;
Ey = Ey/dy;

end

