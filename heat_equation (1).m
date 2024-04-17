clc;
close all;
clear;

%% Constants

a = 1; %thermal diffusivity

dt = 0.0001;    %step time
sim_time = 0.5;  %simulation time
total_itr = floor(sim_time/dt); %number of iterations

grid_size = 100;
actual_size = 2;
h = actual_size/grid_size; % step in space

if (a*dt/h^2 > 1/4)
    error("     a*dt/h^2 = %d\n" + ...
          "     Which is not less than 1/4\n",a*dt/h^2);
end

%% Initial Value

u = zeros(grid_size);

%% Conditions

top_condition = 10*sin(linspace(0,pi,grid_size));
left_condition = zeros(grid_size,1);
right_condition = zeros(grid_size,1);
bottom_condition = 5*sin(linspace(0,pi,grid_size));

u(:,end) = top_condition;
u(1,:) = left_condition;
u(:,1) = bottom_condition;
u(end,:) = right_condition;

%% Solve loop
[X,Y] = meshgrid(linspace(0,actual_size,grid_size),linspace(0,actual_size,grid_size));
surf(X,Y,u);

for i=1:total_itr
    delta_u = my_laplacian(u, h);
    
    % Only changing the center to keep rand conditions
    u_center = u(2:end-1,2:end-1);
    
    % Explicit Euler 
    u(2:end-1,2:end-1) = (u_center + dt * a * delta_u);
    
    % Plot simulation
    if(mod(i,100)==0)
        surf(X,Y,u);
        title(sprintf("%d / %d",i,total_itr));
        pause(0.1)
    end

end


%% Function definitions
function L = my_laplacian(u, step)
    uCenter = u(2:end-1, 2:end-1);
    uLeft   = u(1:end-2, 2:end-1);
    uRigth  = u(3:end,   2:end-1);
    uTop    = u(2:end-1, 3:end);
    uBottom = u(2:end-1, 1:end-2);
    
    L = (uLeft + uRigth + uTop + uBottom - 4*uCenter)/(step^2);
end
