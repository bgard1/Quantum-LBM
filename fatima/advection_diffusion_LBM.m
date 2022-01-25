% January 25th 2022 by Fatima Ezahra Chrit
% D1Q2 Lattice Boltzmann code for: advection-diffusion equation
% (with special case of heat diffusion in slab with const temp)
% (rho is the diffused scalar, advection velocity u is externally imposed)
% c2 - c1

clear all; close all;

flHeatDiffusion = 0;    %flag to solve heat diffusion equation (Dirichlet at T(x=0)=Tw and zero flux at x=L)
Tw = 1;  % wall temp in heat diffusion equation
mu0 = 13;  % location of point source

%domain params
L = 63;    %domain length
dx = 1;     %spacing
nx = L/dx +1; % number of nodes
x = linspace(0,L,nx);

%D1Q2 lattice constants
w=[1/2 1/2]; % weight coeffecients
cx = [1 -1]; %lattice velocities
csq = 1; %square of sound speed

%simulation parameters
if (flHeatDiffusion)
    ux = 0.; %pure diffusion
else
    ux = 0.2; %Advection x-velocity 
end
tau = 1; 
omega = 1/tau;
D = 0.5; %csq * (tau - 1/2);  % diffusion coefficient
maxT  = 100; % total number of iterations 
tPlot = 1;    % cycles for graphical output

% initialisation 
if (flHeatDiffusion)
    rho = zeros(1,nx);
else
    rho = 0.1 * ones(1,nx);
%     rho(mu0-1) = 0.15; 
    rho(mu0) = 0.2;
%     rho(mu0+1) = 0.15;
  
end

%plot initial dstribution
figure; plot(x,rho,'k'); title(['Solution  t = 0']);
xlabel('x'); ylabel('\rho');
% return;

%initialize distribution functons
for i=1:2
    cu = cx(i)*ux / csq; 
    feq(i,:) = rho .* w(i) .* (1 + cu); %+ cu.*cu/2 - ux.^2/(2*csq)); %commented is the non-linear term in feq
end
f = feq;  % distribution function f initialized from the equilibrium function 
rho = sum(f); 

% MAIN LOOP
for cycle = 1:maxT
    clf;     
    
     % COLLISION STEP 
    for i=1:2
        cu = cx(i) * ux / csq; 
        feq(i,:) = rho .* w(i) .* (1 + cu); % + cu.*cu/2 - ux.^2/(2*csq)); %commented is the non-linear term in feq
        fout(i,:) = f(i,:) - omega .* (f(i,:)-feq(i,:)); 
    end 
    
    % STREAMING STEP 
    f(1,2:nx) = fout(1,1:nx-1);
    f(2,1:nx-1) = fout(2,2:nx);
    
    if (flHeatDiffusion)
        %Heat diffusion BC
        f(1,1) = Tw - f(2,1);       %Dirichlet BC
        f(1,nx) = f(1,nx-1);        %zero flux
        f(2,nx) = f(2,nx-1);        %zero flux
    else
        % periodic BC
        f(1,1) = fout(1,nx);
        f(2,nx) = fout(2,1);
    end
    
    %     MACROSCOPIC VARIABLES (moments)
    rho = sum(f);
    
%   analytical solution
    if (flHeatDiffusion)
        rhoAnalytical = Tw *(1-erf(x/(sqrt(4*D*cycle))));
    end

%     % VISUALIZATION
    if (mod(cycle,tPlot)==0)
        plot(x,rho);
        hold on;
        if(flHeatDiffusion)
            hold on
            plot(x, rhoAnalytical, '*r');
            legend('LBM', 'Analytical');
        end
%         title('Solution at t = ' + cycle);
        xlabel('x');
        ylabel('\rho');
%         axis([0 64 0.1 0.12]);
%         axis equal off;
        drawnow
    end
end

%data from Budinski figure (keep same params)
dataBudinski_t20 = [-0.03971	0.10001
4.0081	0.10096
4.93231	0.10143
5.95011	0.10246
6.95622	0.10347
7.97402	0.1053
8.99182	0.10706
10.92214	0.11168
12.94604	0.11588
14.98165	0.11792
16.72478	0.11649
18.93587	0.11242
20.97147	0.10747
23.08897	0.10347
26.03709	0.10076
27.9674	0.10015
32.95113	0.10001
44.93078	0.10001
58.19731	0.10008];

dataBudinski_t40 = [0.05388	0.10008
4.0081	0.10062
5.95011	0.1015
8.99182	0.10347
12.03353	0.10673
14.98165	0.11019
17.09915	0.11201
19.02946	0.11269
21.8021	0.1114
24.01318	0.10937
26.9613	0.10571
30.00301	0.10279
33.04472	0.10096
36.99894	0.10008
41.98266	0.10001
48.15967	0.09994
58.10372	0.09994];

dataBudinski_t60 = [0.05388	0.10001
4.93231	0.10062
7.97402	0.1015
11.01573	0.103
13.96385	0.1051
17.00555	0.10754
19.95367	0.10951
22.99538	0.11045
26.9613	0.1091
30.00301	0.10692
33.04472	0.10449
35.99284	0.10246
40.04065	0.10076
44.10016	0.10015
50.08999	0.09994
57.09761	0.09987];

dataBudinski_t100 = [0.05388	0.10001
4.0081	0.10015
7.97402	0.10049
14.05744	0.10178
19.02946	0.10381
22.99538	0.10571
27.05489	0.10741
31.00911	0.10808
35.06862	0.10747
39.03454	0.10584
42.98877	0.10388
47.04828	0.10225
51.10779	0.10103
56.07981	0.10028
61.14543	0.10001];

if(~flHeatDiffusion)
    switch maxT
        case 20
            plot(dataBudinski_t20(:,1)+1, dataBudinski_t20(:,2),'*-r');
            legend('classical','Budinski');
        case 40
            plot(dataBudinski_t40(:,1)+1, dataBudinski_t40(:,2),'*-r');
            legend('classical','Budinski');
        case 60
            plot(dataBudinski_t60(:,1)+1, dataBudinski_t60(:,2),'*-r');
            legend('classical','Budinski');
        case 100
            plot(dataBudinski_t100(:,1)+1, dataBudinski_t100(:,2),'*-r');
            legend('classical','Budinski');
    end
end    
title(['Solution  t = ', num2str(maxT)]);
