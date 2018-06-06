% Lab 3: 2019862s.

% PDE solver.
function pdeFileNeumannBCs()
m=0;             % Symmetry parameter.
x=0:0.1:400;     % Spatial coordinates.
t=0:0.1:110;      % Time span.

% Solve the system of PDEs by invoking the PDE function,
% Initial Conditions, and Boundary Conditions.
sol=pdepe(m,@diffusionPDE,@diffusionIC,@diffusionBC,x,t); 

% Approximate the u component solution for all t and all x. 
u1=sol(:,:,1);
% Approximate the v component solution for all t and all x.
u2=sol(:,:,2);

% Calculating the wave speed from the numerical solutions.
front_locu=zeros;
for i=9:10
    % Calculate v for t=80 and t=90 for all x. Store if v>0.001.
    front_vecu=u2(i,:)>0.001;
    % Multiply the v vector by the spatial position x
    % to obtain the location of the leading edge (before v=0).
    front_locu(i-8)=max(front_vecu.*x);
end
% Calculate the speed via a finite difference approximation.
speedu=diff(front_locu);
% Find the average wave speed of v by dividing through 
% the time interval.
cv=sum(speedu)/10;
disp(cv);

% Plot of u and v as functions of x for all t.
figure
for j=1:12
    plot(x,u1(j,:),'k');
    hold on
    plot(x,u2(j,:),'k--');
end
title('Travelling waves of u and v. Neumann BCs, a=0.5');
xlabel('Distance x');
ylabel('u(x,t) and v(x,t)');
axis([0, 400, -0.1, 1.1]); legend('u(x,t)','v(x,t)');
% ------------------------------------------------------------
% PDE function - specify parameters, flux, and source terms.
function [c,flux,source] = diffusionPDE(~,~,u,DuDx)
% Parameters:
c=[1;1];   % Coefficients of partial u and v w.r.t. time.
Du=1;    % Diffusion coefficient in u equation.
Dv=2;      % Diffusion coefficient in v equation.
r=1.5;     % Ratio of growth rates

% Interspecies coefficients:
a=0.5;     % u decreases due to interaction with v.
b=0.7;     % v decreases due to interaction with u.

% Define the flux (diffusion term) of the PDEs.
flux=[Du;Dv].*DuDx;

% Define the source (reaction term) of the PDEs.
source=[u(1)*(1-u(1)-a*u(2)); r*u(2)*(1-u(2)-b*u(1))];
% ------------------------------------------------------------
% Initial conditions
function u0 = diffusionIC(x)
% Initiate a vector for the ICs.
u0=[1;1];
% If statement for the value of x.
if x >= 0 && x <= 10    % For x=[0,10],
   u0(2) = 1;           % v = 1.
else                    % For x=(10,400],
   u0(2) = 0;           % v = 0.
end
% ------------------------------------------------------------
% Boundary conditions
function [pl,ql,pr,qr] = diffusionBC(~,~,~,ur,~)
% Neumann boundary conditions. 
% No specifications about
% the values of the functions on the left, this is
% already taken care of by the ICs.
pl=[0;0];  
% Both derivatives w.r.t. x are zero on the left. Neumann.
ql=[1;1];
% The values of u and v on the right are calculated
% numerically by the pdepe sovler, no need to specify BCs.
pr=[0;0];
% Both derivatives w.r.t. x are zero on the right. Neumann.
qr=[1;1];

% % Mixed boundary conditions.
% % Both u and v are already described by ICs on the left.
% pl = [0;0];
% % Both derivatives w.r.t. x are zero on the left. Neumann.
% ql = [1;1];
% % Set u=1 and v=0 on the right boundary. Dirichlet BC.
% pr = [ur(1)-1;ur(2)];
% % Both derivatives w.r.t. x are zero on the right.
% qr = [0;0];