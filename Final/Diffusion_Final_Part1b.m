clear;

% domain
xL = -1;
xR = 1;
yB = -0.5;
yT = 1.7;

% diffusion coefficient
lambda = 0.75;

% source term
f = @(t,x,y) (2*lambda - 1)*(sin(x)*cos(y)*exp(-t));

% initial conditions
T_start = @(x,y) sin(x)*cos(y)*exp(0);

% time interval
t_start = 0;
t_final = 1;

% exact solution
T_exact = @(t,x,y) sin(x)*cos(y)*exp(-t);

% starting matrix size and number of times to loop
Nx = 25;
Ny = 30;
num_splits = 3;
    
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

    % create sparse matrix and allocate memory for right-hand side
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
            
    % entering loop to solve using impilicit method
    % condition to run until pre designated final time
    while t < t_final
    
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
                    RHS(((j-1)*Nx + i)) = T_exact(t+dt,x(i),y(j));
                    
                % internal points
                else
                    RHS(((j-1)*Nx + i)) = c_old(j,i) + dt*f(t+dt,x(i),y(j));
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
fprintf('Resolution \t\tError \t\t Order\n');
    for i = 1:num_splits
        if Nx_all(i) >= 100 || Ny_all(i) >= 100
            fprintf('%gx%g \t %g \t\t %g\n', Nx_all(i), Ny_all(i), error(i), order(i));
        else
        fprintf('%gx%g \t\t %g \t\t %g\n', Nx_all(i), Ny_all(i), error(i), order(i));
        end
    end