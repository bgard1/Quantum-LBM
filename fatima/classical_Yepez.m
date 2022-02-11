clear all; close all;


flBurgers = 0;   %if not, solve diffusion equation

%initial solutions for diffusion equation
flDelta = 0;
flGaussian = 1;
flExp = 0;

%D1Q2 lattice constants
w=[1/2 1/2]; % weight coeffecients
cx = [1 -1]; %lattice velocities
csq = 1; %square of sound speed
ux = 0.;
D = 0.5;

%domain params
L = 64;     %domain length
dx = 1;     %spacing
nx = L/dx +1; % number of nodes
x = linspace(0,L,nx);

sigma0 = (L)/10;
mu0 = ceil(L/2) ;%12;

maxT  = 40; % total number of iterations 
tPlot = 1;    % cycles for graphical output

%number operator
n1 = [0 0 0 0; ...
      0 1 0 0; ...
      0 0 0 0; ...
      0 0 0 1];
n2 = [0 0 0 0; ...
      0 0 0 0; ...
      0 0 1 0; ...
      0 0 0 1]; 
 
%collision operator  
if (flBurgers)
    C = [1 0 0 0; ...
        0 1/sqrt(2) 1/sqrt(2) 0; ...
        0 -1/sqrt(2) 1/sqrt(2) 0; ...
        0 0 0 -1];
else
    C = [1 0 0 0; ...
        0 (1+1i)/2 (1-1i)/2 0; ...
        0 (1-1i)/2 (1+1i)/2 0; ...
        0 0 0 1];
end
    
%initial density
if (flGaussian)
    rho = 1/4 * exp(-(x-mu0).^2 / sigma0^2) + 1/2;
else
    if(flExp)
        rho = 1/4 * sin(2*pi*x/L) + 1/2;
    else
        rho = zeros(1,nx);
        rho(L/2 + 1) = 1; 
    end
end
%initialize distribution functons
for k=1:2
    cu = cx(k)*ux / csq; 
    feq(k,:) = rho .* w(k) .* (1 + cu); 
end
f = feq;  % distribution function f initialized from the equilibrium function 

% MAIN LOOP
for cycle = 1:maxT
    clf; 
    
    %initialize each pair of qubits
    q0 = [sqrt(1-f(1,:)), sqrt(f(1,:))];  %right qubit
    q1 = [sqrt(1-f(2,:)), sqrt(f(2,:))];  %left qubit

    %initial combined state
    for k =1:nx
        initial_state(k,:) = [sqrt( (1-f(1,k))*(1-f(2,k)) ), ...
                         sqrt( f(1,k)*(1-f(2,k)) ), ...
                         sqrt( f(2,k)*(1-f(1,k)) ), ...
                         sqrt( f(1,k)*f(2,k) ) ];
    end
    

    %post-collision state (4xlattice_sites)
    post_collision_state = C * initial_state';

   
    %post-collision distribution
    for k=1:nx
        post_collision_distribution(1,k) = post_collision_state(:,k)' * n1 * post_collision_state(:,k);
        post_collision_distribution(2,k) = post_collision_state(:,k)' * n2 * post_collision_state(:,k);
    end

    % STREAMING STEP 
    f(1,2:nx) = post_collision_distribution(1,1:nx-1);
    f(2,1:nx-1) = post_collision_distribution(2,2:nx);

     % periodic BC
    if (flExp)
        f(1,1) = 1/2 - f(2,1); % x = 0, rho = 1/2.
        f(2,nx) = 1/2 - f(1,nx); % x = L, rho = 1/2.
    else
        f(1,1) = f(1,nx);
        f(2,nx) = f(2,1);
    end
    % MACROSCOPIC VARIABLES (moments)
    rho = sum(f);
    
    %analytical solution
    if(flGaussian)
        rhoAnalytical = 1/4 * sigma0 / sqrt(sigma0^2 + 4*D*(cycle)) * exp(-(x-mu0).^2 / (sigma0^2 + 4*D*(cycle))) + 1/2;
    else
        if(flExp)
            kWave = 2*pi/L;
            rhoAnalytical =1/4 *exp(-D*kWave^2*cycle) * sin(2*pi*x/L)+ 1/2;
        else
            rhoAnalytical = 1/sqrt(4*pi*D*cycle) * exp(- (x-L/2).^2 / (4*D*cycle));
        end
    end

      % VISUALIZATION
    if (mod(cycle,tPlot)==0)
        plot(x,rho);
%         if(flGaussian || flExp)
            hold on;
            plot(x, rhoAnalytical, '*r');
            legend('classical Yepez', 'Analytical');
%             axis([0 L 0.5 0.75 ]);
%         end
        xlabel('x', 'FontSize', 16);
        ylabel('\rho', 'FontSize', 16); 
        if(flDelta)
            axis([0 L 0 1 ]);
        end
        drawnow
    end
                
end
set(gca,'FontSize',16);
title(['Solution  t = ', num2str(maxT)]);