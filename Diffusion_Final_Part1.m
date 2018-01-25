clear;

% domain
xL = 0;
xR = 1;

% diffusion coefficient
D = 0.6;

% boundary conditions
cL = @(t) 0;
cR = @(t) 0;

% source term
f = @(t,x) 0;

% initial conditions
x0 = 0.5;
tau = 0.001;
c_start = @(x) sqrt(4*pi*D*tau)/sqrt(4*pi*D*(tau)) * exp(-(x-x0)^2/(4*D*(tau)));

% time interval
t_start = 0;
t_final = 0.005;

% exact solution
tau = 0.001;
c_exact = @(t,x) sqrt(4*pi*D*tau)/sqrt(4*pi*D*(t+tau)) * exp(-(x-x0)^2/(4*D*(t+tau)));

% space discretization
Nx = 100;
x = linspace(xL, xR, Nx);
dx = (xR-xL)/(Nx-1);

% time-step is unconditionally stable for implicit method
% so I chose a time-step that would be unstable for 
% explicit method to demonstrate
dt = 5*dx*dx/2/D;

% pre-allocate arrays for previous calculated time step, new time step, 
% and exact solutions
c_old = zeros(1,Nx);
c_new = zeros(1,Nx);
c_exa = zeros(1,Nx);

% assign values to array that will hold solution at previous time steps
for i = 1:Nx
    c_old(i) = c_start(x(i));
end

% initialize t to start time
t = t_start;

% create sparse matrix and allocate memory for right-hand side
A = sparse(Nx,Nx);
RHS = zeros(Nx,1);

% calculate the matrix before the while-loop to save time
% internal points
for i = 2:Nx-1
    A(i,i) = 1+2*dt*D/dx/dx;
    A(i,i-1) = -dt*D/dx/dx;
    A(i,i+1) = -dt*D/dx/dx;
end

% boundary points
A(1,1) = 1;
A(Nx,Nx) = 1;
    
% condition to run until pre designated final time
while t < t_final
    
    % make sure t_final is never passed 
    % shorten dt if next time step were to "step" out of bounds
    if t + dt > t_final
        dt = t_final-t;
            
        % need to recalculate the matrix since dt has changed
        % internal points
        for i = 2:Nx-1
            A(i,i) = 1+2*dt*D/dx/dx;
            A(i,i-1) = -dt*D/dx/dx;
            A(i,i+1) = -dt*D/dx/dx;
        end
    end
    
    % internal points
    for i = 2:Nx-1
        RHS(i) = c_old(i) + dt*f(t+dt,x(i));
    end
    
    % boundary points
    RHS(1) = cL(t+dt);
    RHS(Nx) = cR(t+dt);
    
    % solve system of equations
    c_new = A\RHS;
    
    % assigning c_old to c_new so that it calculates using 
    % previous solution
    c_old = c_new;
    
    % go to next time step
    t = t+dt;
    
    % intialize c_exa to the exact solution
    for i = 1:Nx
        c_exa(i) = c_exact(t,x(i));
    end
    
    % plot comparing the exact solution to calculated solution
    % exact appears as a line and calculated as circles
    plot(x,c_exa,'LineWidth',2);
    hold on
    plot(x,c_new,'o-','LineWidth',1);
    hold off
    xlabel('x');
    ylabel('c');
    axis([xL xR 0 1]);
    pause(dt);
    
end