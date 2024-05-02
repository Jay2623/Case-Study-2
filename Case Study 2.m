Lx = 1;
Ly = 1;
Lz = 1;
Nx = 50;
Ny = 50;
Nz = 50;
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
dz = Lz / (Nz - 1);
dt = 0.01;
t_final = 1;
alpha = 0.1;
T_initial = 0;
T_left = 100;
T_right = 0;
T_bottom = 0;
T_top = 0;
T_front = 0;
T_back = 0;
T = ones(Nx, Ny, Nz) * T_initial;
T(1,:,:) = T_left;
T(end,:,:) = T_right;
T(:,1,:) = T_bottom;
T(:,end,:) = T_top;
T(:,:,1) = T_front;
T(:,:,end) = T_back;
t = 0;
while t < t_final
T_new = T;
for i = 2:Nx-1
for j = 2:Ny-1
for k = 2:Nz-1
T_new(i,j,k) = T(i,j,k) + alpha * dt * ((T(i+1,j,k) - 2*T(i,j,k) +
T(i-1,j,k)) / dx^2 + ...
(T(i,j+1,k) - 2*T(i,j,k) + T(i,j-1,k)) / dy^2 + ...
(T(i,j,k+1) - 2*T(i,j,k) + T(i,j,k-1)) / dz^2);
end
end
end
T = T_new;
t = t + dt;
end
[X,Y,Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
slice(X, Y, Z, T, Nx/2, Ny/2, Nz/2);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Temperature Distribution');
colorbar;