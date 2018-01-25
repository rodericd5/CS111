clear;

% Domain
xL = 0;
xR = 12;
yB = 0;
yT = 3;

% Diffusion coefficient
D = 0.2;

% Velocity field
v_x = -0.8;
v_y = -0.4;

% Initial condition
T_start = @(x,y) 0;

% Time interval
t_start = 0;
t_final = 10;

% Exact solution
T_exact = @(t,x,y) 0;

% Robin boundary condition 
g = @(t,x,y) 0;

% Setting up number of grid
Nx = 160;
Ny = 40;

% Variables holding counts for plots and beach numbers
plot_num = 1;
t_count = 1;
t_count_minus_one = 1;

% Space discretization    
x = linspace(xL, xR, Nx);
y = linspace(yB, yT, Ny);
dx = (xR-xL)/(Nx-1);
dy = (yT-yB)/(Ny-1);

% Time-step 
dt = .1;

% pre allocating arrays for time values and the
% associated oil concentration levels at each beach 
% which will determine when to close them
t_vals = zeros(1, floor((t_final - t_start)/dt)+2);
beach_close = zeros(1, floor((t_final - t_start)/dt)+1);
beach_close_two = zeros(1, floor((t_final - t_start)/dt)+1);
beach_close_three = zeros(1, floor((t_final - t_start)/dt)+1);

% initializing beach locations
beach_plot_point = 4;
beach_plot_point_two = 6;
beach_plot_point_three = 8;

% intializing initial concentration of oil at each beach
beach_close(1) = 0;
beach_close_two(1) = 0;
beach_close_three(1) = 0;

% Allocating arrays for finding concentration using schemes below
c_old = zeros(Ny,Nx);
c_new = zeros(Ny,Nx);

% Initializing c_old to T_start values
for i = 1:Nx
    for j = 1:Ny
        c_old(j,i) = T_start(x(i),y(j));
    end
end

% Setting the time
t = t_start;

% Create sparse matrix and allocate memory for right-hand side
A = sparse(Nx*Ny,Nx*Ny);
RHS = zeros(Nx*Ny,1);

% Values used to calculate internal points in matrix A
C = 1 + D*((2*dt)/(dx^2) + (2*dt)/(dy^2));
R = -(D*dt)/(dx^2);
L = -(D*dt)/(dx^2);
B = -(D*dt)/(dy^2);
T = -(D*dt)/(dy^2);

% Values used to calculate robin BC points in matrix A
ac = C + 2*v_y*dt/dy;
ar = R;
al = L;
at = T + B;            
            
% Figure for plotting oil concentration at 3 different times
figure('rend','painters','pos',[400 100 600 1000])    
    
    % Making sure not to overstep t_final
    while t < t_final
        
        % Plotting concentration of oil at times 1,3, and 7
        if round(t,1) == 1 || round(t,3) == 3 || round(t,7) == 7
            subplot(3,1,plot_num)
            contourf(x,y,c_old, 100, 'LineColor','none');
            title(['      Oil concentration at t = ' num2str(t)]);
            colorbar;
            xlabel('x-axis');
            ylabel('y-axis');
            caxis([0,.025]);
            hold on
            plot_num = plot_num + 1;
        end
        
        % assigning intial values to matrix
        for i = 1:Nx
            for j = 1:Ny
                
                % assigning p to (j-1)*Nx + i for readability
                p = (j-1)*Nx + i;
                
                % boundary points
                if i == 1 || i == Nx || j == Ny
                    A(p,p) = 1;
                    
                % Robin boundary condition points
                elseif j == 1
                    A(p,p) = ac;
                    A(p,p - 1) = al;
                    A(p,p + 1) = ar;
                    A(p,p + Nx) = at;
                    
                % internal points
                else
                    A(p,p) = C;
                    A(p,p - 1) = L;
                    A(p,p + 1) = R;
                    A(p,p - Nx) = B;
                    A(p,p + Nx) = T;
                end
            end
        end
        
    % make sure t_final is never passed 
    % shorten dt if next time step were to "step" out of bounds
    if t + dt > t_final
        dt = t_final-t;
            
        % initializing values for solving internal points of the matrix
        C = 1 + D*((2*dt)/(dx^2) + (2*dt)/(dy^2));
        R = -(D*dt)/(dx^2);
        L = -(D*dt)/(dx^2);
        B = -(D*dt)/(dy^2);
        T = -(D*dt)/(dy^2);
            
        % initializing values for solving Robin boundary points
        ac = C + 2*v_y*dt/dy;
        ar = R;
        al = L;
        at = T + B;
            
        % need to recalculate the matrix since dt has changed
        for i = 1:Nx
            for j = 1:Ny  
                p = (j-1)*Nx + i;
                
                % boundary points
                if i == 1 || i == Nx || j == Ny
                    A(p,p) = 1;
                    
                % Robin condition boundary points
                elseif j == 1
                    A(p,p) = ac;
                    A(p,p - 1) = al;
                    A(p,p + 1) = ar;
                    A(p,p + Nx) = at;
                   
                % internal points
                else
                    A(p,p) = C;
                    A(p,p - 1) = L;
                    A(p,p + 1) = R;
                    A(p,p - Nx) = B;
                    A(p,p + Nx) = T;
                end
            end
        end
    end
   
    % calculating rhs
    for i = 1:Nx
        for j = 1:Ny
            p = (j-1)*Nx + i;
            
            % boundary points
            if i == 1 || i == Nx || j == Ny
                RHS(p) = T_exact(t+dt, x(i), y(j));
            
            % Robin condition boundary points
            elseif j == 1
                RHS(p) = c_old(j,i) + dt*f(t+dt,x(i),y(j))-v_x*dt*(c_old(j,i+1)-c_old(j,i))/dx-v_y*dt*(c_old(j+1,i)-c_old(j,i))/dy-(2*dt/dy)*g(t+dt,x(i),y(1));
            
            % internal points
            else          
                RHS(p) = c_old(j,i) + dt*f(t+dt,x(i),y(j))-v_x*dt*(c_old(j,i+1)-c_old(j,i))/dx-v_y*dt*(c_old(j+1,i)-c_old(j,i))/dy;
            end
        end
    end
    
    
    % Solve system of equations
    c_new = reshape(A\RHS,Nx,Ny)';

    % assigning c_old to c_new so that it calculates using 
    % previous solution
    c_old = c_new;
      
    % place the current time in array of time values
    t_vals(t_count) = t;
    
    % increment t_count so that we get the next values of time and temp
    % for our two arrays that we fill
    t_count = t_count + 1;
        
    % place the solutions of oil concentration in each array
    % representing its respective beach
    beach_close(t_count) = c_new(1,floor(beach_plot_point/dx+1));
    beach_close_two(t_count) = c_new(1, floor((beach_plot_point_two)/dx+1));
    beach_close_three(t_count) = c_new(1, floor((beach_plot_point_three)/dx+1));        
      
    % go to next time step
    t = t+dt;
        
end    
    
    % plot oil concentration at beach 1   
    figure;
    plot(t_vals, beach_close,'-')
    title('Temperature vs Oil Concentration at Beach 1');
    xlabel('Time');
    ylabel('Concentration');
    
    % plot oil concentration at beach 2
    figure;
    plot(t_vals, beach_close_two,'-')
    title('Temperature vs Oil Concentration at Beach 2');
    xlabel('Time');
    ylabel('Concentration');
    
    % plot oil concentration at beach 3
    figure;
    plot(t_vals, beach_close_three,'-')
    title('Temperature vs Oil Concentration at Beach 3');
    xlabel('Time');
    ylabel('Concentration');
    
% simple function to determine source term
function oilSpill = f(t,x,y)

xs = 10;
rs = 0.1;
e = 0.1;

    if t > 0.5
        oilSpill = 0;
    else
        oilSpill = 1/2*(1-tanh((sqrt((x-xs)^2+y^2)-rs)/e));
    end
end