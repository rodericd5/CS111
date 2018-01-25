clear;

% domain
xL = -1;
xR = 3;
yB = -1.5;
yT = 1.5;

% diffusion coefficient
D = 0.7;

% velocity field
v_x = -0.8;
v_y = -0.4;

% source term
f = @(t,x,y) exp(-t)*((-sin(x)*cos(y))+(v_x*cos(x)*cos(y))-(v_y*sin(x)*sin(y))+(2*D*sin(x)*cos(y)));

% initial conditions
T_start = @(x,y) sin(x)*cos(y)*exp(0);

% time interval
t_start = 0;
t_final = 1;

% exact solution
T_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% boundary conditions
g = @(t,x,y) D*(-sin(x)*sin(y)*exp(-t)) - v_y*sin(x)*cos(y)*exp(-t);

% starting matrix size and number of times to loop
Nx = 20;
Ny = 15;

% number of times to run through
num_splits = 4;

% will run through this 3 times and increase matrix size then calculate
% error
for k = 1:num_splits 
    
    % space discretization
    x = linspace(xL, xR, Nx);
    y = linspace(yB, yT, Ny);
    dx = (xR-xL)/(Nx-1);
    dy = (yT-yB)/(Ny-1);

    % time-step as given
    dt = .5*dx;
    
    % pre-allocate arrays for previous calculated time step, new time step, 
    % and exact solutions
    c_old = zeros(Ny,Nx);
    c_new = zeros(Ny,Nx);
    c_exa = zeros(Ny,Nx);

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
    C = 1 + D*((2*dt)/(dx^2) + (2*dt)/(dy^2));
    R = -(D*dt)/(dx^2);
    L = -(D*dt)/(dx^2);
    B = -(D*dt)/(dy^2);
    T = -(D*dt)/(dy^2);

    % initializing values for solving boundary points that abide by
    % the Robin boundary condition
    ac = C + 2*v_y*dt/dy;
    ar = R
    al = L
    at = T + B;            
            

    while t < t_final
        
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
            ar = R
            al = L
            at = T + B;
            
            % need to recalculate the matrix since dt has changed
            for i = 1:Nx
                for j = 1:Ny  
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
        end
    
        % calculating the rhs
        for i = 1:Nx
            for j = 1:Ny
                p = (j-1)*Nx + i;
                
                % boundary points
                if i == 1 || i == Nx || j == Ny
                    RHS(p) = T_exact(t+dt, x(i), y(j));
                    
                % Robin boundary condition points
                elseif j == 1
                    RHS(p) = c_old(j,i) + dt*f(t+dt,x(i),y(j))-v_x*dt*(c_old(j,i+1)-c_old(j,i))/dx-v_y*dt*(c_old(j+1,i)-c_old(j,i))/dy-(2*dt/dy)*g(t+dt,x(i),y(1));
                
                % internal points
                else
                    RHS(p) = c_old(j,i) + dt*f(t+dt,x(i),y(j))-v_x*dt*(c_old(j,i+1)-c_old(j,i))/dx-v_y*dt*(c_old(j+1,i)-c_old(j,i))/dy;
                end
            end
        end
    

    
        % solve system of equations
        c_new = reshape(A\RHS,Nx,Ny)';
    
        % assigning c_old to c_new so that it calculates using 
        % previous solution
        c_old = c_new;
        
        % go to next time step
        t = t+dt;
    
    end
    
    % intialize c_exa to the exact solution
    for i = 1:Nx
        for j = 1:Ny
            c_exa(j,i) = T_exact(t,x(i),y(j));
        end
    end
    
    % find error based on the max absolute deviation for each size matrix
    error(k) = max(max(abs(c_exa - c_new)));
    
    % create list of matrix dimensions used 
    Nx_all(k) = Nx;
    Ny_all(k) = Ny;

    % resize matrix for next run through    
    Nx=Nx*2;
    Ny=Ny*2;
    
end

    % calculating order of error
    for i = 2:num_splits
        order(i) = log(error(i-1)/error(i))/log(2);
    end
    
% print error and order of accuracy
fprintf('Resolution \tError \t\t\t Order\n');
    for i = 1:num_splits
        if Nx_all(i) >= 100 || Ny_all(i) >= 100
            fprintf('%gx%g \t %g \t\t %g\n', Nx_all(i), Ny_all(i), error(i), order(i));
        else
        fprintf('%gx%g \t\t %g \t\t %g\n', Nx_all(i), Ny_all(i), error(i), order(i));
        end
    end