clear;

% domain
xL = -2;
xR = 2;
yB = -2.5;
yT = 2.5;

% diffusion coefficient
lambda = 0.0015;

% source term
f = @(t,x,y) 0;

% initial conditions
T_start = @(x,y) 20;

% time interval
t_start = 0;
t_final = 1500;

% boundary conditions
T_bc = @(t,x,y) min(20+8*t/6, 100);

% size of matrix
Nx = 80;
Ny = 100;

% count used to allocate array of times
t_count = 1;    

% space discretization
x = linspace(xL, xR, Nx);
y = linspace(yB, yT, Ny);
dx = (xR-xL)/(Nx-1);
dy = (yT-yB)/(Ny-1);

% time-step 
dt = 5;

% holds number of plots
plot_num = 1;

% pre-allocate arrays for previous calculated time step, new time step, 
% and exact solutions 
c_old = zeros(Ny,Nx);
c_new = zeros(Ny,Nx);
c_exa = zeros(Ny,Nx);

% pre-allocate arrays for temperature at center of potato and time values
c_center = zeros(1, (t_final - t_start)/dt);
t_vals = zeros(1, (t_final - t_start)/dt);

% assign values to array that will hold solution at previous time steps
for i = 1:Nx
    for j = 1:Ny
        c_old(j,i) = T_start(x(i),y(j));
    end
end

% initialize t to start time
t = t_start;

% Create sparse matrix and allocate memory for right-hand side
A = sparse(Nx*Ny,Nx*Ny);
RHS = zeros(Nx*Ny,1);

% initializing values for solving internal points of the matrix
% C - internal points in the center of matrix
% R - internal points on the right of matrix
% L - internal points on the left of matrix
% B - internal points on the bottom of matrix
% T - internal points on the top of matrix
C = 1 + lambda*((2*dt)/(dx^2) + (2*dt)/(dy^2));
R = -(lambda*dt)/(dx^2);
L = -(lambda*dt)/(dx^2);
B = -(lambda*dt)/(dy^2);
T = -(lambda*dt)/(dy^2);

% assigning intial values to matrix
for i = 1:Nx
    for j = 1:Ny
        
        % boundary points
        if i == 1 || i == Nx || j == 1 || j == Ny
            A(((j-1)*Nx + i),j) = 0;
            A(((j-1)*Nx + i),((j-1)*Nx + i)) = 1;
            
        % internal points
        else
            A(((j-1)*Nx + i),((j-1)*Nx + i)) = C;
            A(((j-1)*Nx + i),((j-1)*Nx + i) - 1) = L;
            A(((j-1)*Nx + i),((j-1)*Nx + i) + 1) = R;
            A(((j-1)*Nx + i),((j-1)*Nx + i) - Nx) = B;
            A(((j-1)*Nx + i),((j-1)*Nx + i) + Nx) = T;
        end
    end
end

% create a figure for the subplots to go on
figure('rend','painters','pos',[400 100 600 1000])    

% entering loop to solve using impilicit method
% condition to run until pre designated final time
while t < t_final
    
    % create subplots of potato temperature at 0s 200s 400s and 600s
    if (t == 0) || (t == 200) || (t == 400) || (t == 600)
        subplot(4,2,plot_num)
        contourf(x,y,c_old, 100, 'LineColor','none');
        title(['      Temperature of potato at t = ' num2str(t)]);
        colorbar;
        caxis([20,100]);
        hold on
        plot_num = plot_num + 1;
    end
    
    % make sure t_final is never passed 
    % shorten dt if next time step were to "step" out of bounds
    if t + dt > t_final
        dt = t_final-t;
        
        % initializing values for solving internal points of the matrix
        C = 1 + lambda*((2*dt)/(dx^2) + (2*dt)/(dy^2));
        R = -(lambda*dt)/(dx^2);
        L = -(lambda*dt)/(dx^2);
        B = -(lambda*dt)/(dy^2);
        T = -(lambda*dt)/(dy^2);
        
        % need to recalculate the matrix since dt has changed
        for i = 1:Nx
            for j = 1:Ny
                
                % boundary points
                if i == 1 || i == Nx || j == 1 || j == Ny
                    A(((j-1)*Nx + i),j) = 0;
                    A(((j-1)*Nx + i),((j-1)*Nx + i)) = 1;
                    
                % internal points
                else
                    A(((j-1)*Nx + i),((j-1)*Nx + i)) = C;
                    A(((j-1)*Nx + i),((j-1)*Nx + i) - 1) = L;
                    A(((j-1)*Nx + i),((j-1)*Nx + i) + 1) = R;
                    A(((j-1)*Nx + i),((j-1)*Nx + i) - Nx) = B;
                    A(((j-1)*Nx + i),((j-1)*Nx + i) + Nx) = T;
                end
            end
        end
    end
    
    % calculating the rhs
    for i = 1:Nx
        for j = 1:Ny
            
            % boundary points
            if i == 1 || i == Nx || j == 1 || j == Ny
                RHS(((j-1)*Nx + i)) = T_bc(t+dt,x(i),y(j));
                
            % internal points
            else
                RHS(((j-1)*Nx + i)) = c_old(j,i) + dt*f(t+dt,x(i),y(j));
            end
        end
    end
    
    % solve system of equations
    c_new = reshape(A\RHS,Nx,Ny)';
    
    % place the current time in array of time values
    t_vals(t_count) = t;
    
    % place the temperature solutions at the center (middle of potato)
    % in our array of center temperatures at the current time starting
    % with the start temperature
    c_center(t_count) = c_old(Ny/2, Nx/2);
    
    % assigning c_old to c_new so that it calculates using 
    % previous solution
    c_old = c_new;
    
    % go to next time step
    t = t+dt;
    
    % increment t_count so that we get the next values of time and temp
    % for our two arrays that we fill
    t_count = t_count + 1;    
end

% plot a 5th subplot of the temperature vs time at the center of potato
subplot(4,2,[5,8]);
plot(t_vals, c_center, '-')
title('Temperature vs Time at Center of Potato');
text(890,65.03,'\leftarrow 65.03^{\circ}C at t=890', 'fontsize', 14)
xlabel('Time');
ylabel('Temperature');
 