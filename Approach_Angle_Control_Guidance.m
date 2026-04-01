clc
clear
close all


%% Initial and Final of Dynamic Components
x_i = 0;  %km
y_i = 0;  %km
x_f = 50e3; %km
y_f = 50e3; %km
gam_i = deg2rad(45); %rad
gam_f = deg2rad(90); %rad

V = 300;
g0 = 9.8;
%% No Fly zone 1:
xc1  = 30e3;
yc1  = 15e3;
r1   = 10e3;

%% No Fly zone 2:
xc2  = 20e3;
yc2  = 35e3;
r2   = 10e3;

%% No Fly zone 3:
xc3  = 35e3;
yc3  = 30e3;
r3   = 5e3;

%% No Fly zone 4:
xc4  = 5e3;
yc4  = 15e3;
r4   = 10e3;

%% Boundaries of dynamic components
u_max = +deg2rad(80);
u_min = -deg2rad(80);

%% Trust region defined
eps = 1*1e3;

%% Split Node
N   =  100;                  %Node
tf  =   sqrt(x_f^2+y_f^2)/V; %Seconds
t0  =   0;                   %Seconds
dt  =   tf / N;
dt2 =   dt/2;

%% number node of slack variable 
% ns1x = N;
% ns2x = N;
% ns1y = N;
% ns2y = N;
% ns1g = N;
% ns2g = N;
ns1NFZ = N+1;
ns2NFZ = N+1;
% ns     = ns1x + ns2x + ns1y + ns2y + ns1g + ns2g + ns1NFZ + ns2NFZ;
ns     = ns1NFZ + ns2NFZ;
% ns     = ns1NFZ;
% ns     = ns1x + ns2x + ns1y + ns2y + ns1g + ns2g;
xi     = 1e3;       % Penaltial weight

%% number node of dynamic variable
nx      = N+1;  
ny      = N+1;
ngam    = N+1;
nu      = N+1;
neta    = N+1;
ntf     = 1;
% nvar    = nx + ny + ngam + nu + neta + ntf + ns;    %Within eta and slack
% nvar    = nx + ny + ngam + nu + neta + ntf;       %Within eta and without slack
% nvar    = nx + ny + ngam + nu + 1;                %Without eta
nvar    = nx + ny + ngam + nu + ntf + ns;

%% Position Dynamic components 
Nx      = 0;
Ny      = Nx + nx;
Ngam    = Ny + ny;
Nu      = Ngam + ngam;
Neta    = Nu + nu;
% Ntf     = Neta + neta + ntf;                  %Within eta
Ntf     = Nu + nu + ntf;                      %Without eta

%% position Slack variables
Ns      = Ntf;
% Ns1x    = Ns;
% Ns2x    = Ns1x  +   ns1x;
% Ns1y    = Ns2x  +   ns2x;
% Ns2y    = Ns1y  +   ns1y;
% Ns1g    = Ns2y  +   ns2y;
% Ns2g    = Ns1g  +   ns1g;
% Ns1NFZ  = Ns2g  +   ns2g;
% Ns2NFZ  = Ns1NFZ  +   ns1NFZ;

Ns1NFZ    = Ns;
Ns2NFZ    = Ns1NFZ + ns1NFZ;

%% Max iteration and tolerance
max_scp_iter = 30;
tol = 1e-4;

%% Random of Dynamic Components (Linear Guess)
x       = linspace(x_i,x_f,nx);
y       = linspace(y_i,y_f,ny);
gam     = linspace(gam_i,gam_f,ngam);
u       = zeros(nu,1)';
eta     = zeros(neta,1)';

%% history
x_hist      = zeros(nx,max_scp_iter);
y_hist      = zeros(ny,max_scp_iter);
gam_hist    = zeros(ngam,max_scp_iter);
u_hist      = zeros(nu,max_scp_iter);
eta_hist    = zeros(neta,max_scp_iter);
tf_hist     = zeros(1,max_scp_iter);
Res_hist    = [];

%% Store first history of random initial dynamic components
x_hist(:,1)     = x;
y_hist(:,1)     = y;
gam_hist(:,1)   = gam;
u_hist(:,1)     = u;
eta_hist(:,1)   = eta;
tf_hist(:,1)    = tf;

%% Calculate Residual of Stationary condition for initial guess
Res_x_Vec       = abs((x(2:end) - x(1:end-1))/(dt2) - fx(V, gam(1:end-1)) - fx(V, gam(2:end)));
Res_y_Vec       = abs((y(2:end) - y(1:end-1))/(dt2) - fy(V, gam(1:end-1)) - fy(V, gam(2:end)));
Res_gam_Vec     = abs((gam(2:end) - gam(1:end-1))/(dt2) - fgam(V, u(1:end-1)) - fgam(V, u(2:end)));
Matrix_Res = [Res_x_Vec Res_y_Vec Res_gam_Vec];
Residual = max(max(Matrix_Res, [], 2));
Res_hist = [Res_hist; Residual];
fprintf('SCP iter %d residual = %g\n', 1, Residual);


%% Calculate 'Constraint inequality'
lb = -Inf.*ones(nvar, 1);
ub = +Inf.*ones(nvar, 1);

lb(Ns1NFZ+1:Ns1NFZ+ns1NFZ) = 0;
lb(Ns2NFZ+1:Ns2NFZ+ns2NFZ) = 0;

lb(Nx+1:Nx+nx)    = x_i;
ub(Nx+1:Nx+nx)    = x_f;

lb(Ny+1:Ny+ny)    = y_i;
ub(Ny+1:Ny+ny)    = y_f;
% lb(Nu+1:Nu+nu)       = u_min;       ub(Nu+1:Nu+nu)       = u_max;

%% Cost function:
% f = zeros(nvar,1); 
% f(Neta+1:Neta+neta,1) = 1*dt;
f(Ns1NFZ+1 : Ns1NFZ+ns1NFZ) = +xi;  % Within Slack variable
f(Ns2NFZ+1 : Ns2NFZ+ns2NFZ) = +xi;  % Within Slack variable
% f(Ns1x+1 : Ns1x+ns1x) = +xi;
% f(Ns2x+1 : Ns2x+ns2x) = +xi;
% f(Ns1y+1 : Ns1y+ns1y) = +xi;
% f(Ns2y+1 : Ns2y+ns2y) = +xi;
% f(Ns1g+1 : Ns1g+ns1g) = +xi;
% f(Ns2g+1 : Ns2g+ns2g) = +xi;
H = zeros(nvar,nvar);
for i = 1:nu
   H(Nu+i,Nu+i) = 2*dt;
end

%% Perform the calculation in each loop (Sequential Programming)
for iter = 2:max_scp_iter %First interation is initial guess for all dynamic components
    dt      = tf / N;
    dt2     = dt / 2;
    fprintf('START AT ITERATION %d\n', iter);

    %% Constraint equivality:
    Aeq = [];
    beq = [];
    A_zone = [];
    b_zone = [];
    for i = 1:N
        %% Dynamic component: fx(V,gam)
        row_x      = zeros(1,nvar);
        row_x(Nx + i)       = - (1./dt2);
        row_x(Nx + i + 1)   = + (1./dt2);
        row_x(Ngam + i)     = - dfx_dgam(V,gam(i));
        row_x(Ngam + i + 1) = - dfx_dgam(V,gam(i+1));
        row_x(Ntf)          = - (fx(V,gam(i)) + fx(V,gam(i+1)))./(tf-t0);
        
        %Slack variable constraint
%         row_x(Ns1x + i)     = + (1./dt2);
%         row_x(Ns2x + i)     = - (1./dt2);
        
        beq_x_0             = - (fx(V,gam(i)) + fx(V,gam(i+1))).*t0./(tf-t0);
        beq_x_1             = - dfx_dgam(V,gam(i))*gam(i) - dfx_dgam(V,gam(i+1))*gam(i+1);
        
        Aeq = [Aeq; row_x];
        beq = [beq; beq_x_0 + beq_x_1];
        
        %% Dynamic component: fy(V,gam)
        row_y      = zeros(1,nvar);
        row_y(Ny + i)       = - (1./dt2);
        row_y(Ny + i + 1)   = + (1./dt2);
        row_y(Ngam + i)     = - dfy_dgam(V,gam(i));
        row_y(Ngam + i + 1) = - dfy_dgam(V,gam(i+1));
        row_y(Ntf)          = - (fy(V,gam(i)) + fy(V,gam(i+1)))./(tf-t0);
        
        %Slack variable constraint
%         row_y(Ns1y + i)     = + (1./dt2);
%         row_y(Ns2y + i)     = - (1./dt2);
        
        beq_y_0             = - (fy(V,gam(i)) + fy(V,gam(i+1))).*t0./(tf-t0);
        beq_y_1             = - dfy_dgam(V,gam(i))*gam(i) - dfy_dgam(V,gam(i+1))*gam(i+1);
        
        Aeq = [Aeq; row_y];
        beq = [beq; beq_y_0 + beq_y_1];
        
        %% Dynamic component: fgam(V,u)
        row_gam      = zeros(1,nvar);
        row_gam(Ngam + i)       = - (1./dt2);
        row_gam(Ngam + i + 1)   = + (1./dt2);
        row_gam(Nu   + i)       = -1/V;
        row_gam(Nu   + i + 1)   = -1/V;
        
        %Slack variable constraint
%         row_gam(Ns1g + i)     = + (1./dt2);
%         row_gam(Ns2g + i)     = - (1./dt2);
        
        Aeq = [Aeq; row_gam];
        beq = [beq; 0];
        if i>1
            %% No fly zone 1 constraint:
            A_zone_1 = zeros(1,nvar);
            b_zone_1 = zeros(nvar,1);
            A_zone_1(Nx + i)  = 2.*(x(i)-xc1);
            A_zone_1(Ny + i)  = 2.*(y(i)-yc1);
            A_zone_1(Ns1NFZ + i) = +1;
            A_zone_1(Ns2NFZ + i) = -1;

            b_zone_1 = r1.^2 - (xc1.^2 - x(i).^2 + yc1.^2 - y(i).^2);

            A_zone = [A_zone; -A_zone_1];
            b_zone = [b_zone; -b_zone_1];
            
            %% No fly zone 2 constraint:
            A_zone_2 = zeros(1,nvar);
            b_zone_2 = zeros(nvar,1);
            A_zone_2(Nx + i)  = 2.*(x(i)-xc2);
            A_zone_2(Ny + i)  = 2.*(y(i)-yc2);
            A_zone_2(Ns1NFZ + i) = +1;
            A_zone_2(Ns2NFZ + i) = -1;

            b_zone_2 = r2.^2 - (xc2.^2 - x(i).^2 + yc2.^2 - y(i).^2);

            A_zone = [A_zone; -A_zone_2];
            b_zone = [b_zone; -b_zone_2];
            
            %% No fly zone 3 constraint:
            A_zone_3 = zeros(1,nvar);
            b_zone_3 = zeros(nvar,1);
            A_zone_3(Nx + i)  = 2.*(x(i)-xc3);
            A_zone_3(Ny + i)  = 2.*(y(i)-yc3);
            A_zone_3(Ns1NFZ + i) = +1;
            A_zone_3(Ns2NFZ + i) = -1;

            b_zone_3 = r3.^2 - (xc3.^2 - x(i).^2 + yc3.^2 - y(i).^2);

            A_zone = [A_zone; -A_zone_3];
            b_zone = [b_zone; -b_zone_3];
            
            %% No fly zone 4 constraint:
            A_zone_4 = zeros(1,nvar);
            b_zone_4 = zeros(nvar,1);
            A_zone_4(Nx + i)  = 2.*(x(i)-xc4);
            A_zone_4(Ny + i)  = 2.*(y(i)-yc4);
            A_zone_4(Ns1NFZ + i) = +1;
            A_zone_4(Ns2NFZ + i) = -1;

            b_zone_4 = r4.^2 - (xc4.^2 - x(i).^2 + yc4.^2 - y(i).^2);

            A_zone = [A_zone; -A_zone_4];
            b_zone = [b_zone; -b_zone_4];
        end
    end
    
%     %% No fly zone 1 constraint at final point
%     A_zone_f_1 = zeros(1,nvar);
%     A_zone_f_1(Nx + nx)  = 2.*(x(nx)-xc1);
%     A_zone_f_1(Ny + ny)  = 2.*(y(ny)-yc1);
%     A_zone_f_1(Ns1NFZ + nx) = +1;
%     A_zone_f_1(Ns2NFZ + ny) = -1;
%     b_zone_f_1           = r1.^2 - (xc1.^2 - x(nx).^2 + yc1.^2 - y(ny).^2);
%     A_zone = [A_zone; A_zone_f_1];
%     b_zone = [b_zone; b_zone_f_1];
    
    
    
    %% Boundary conditions: x_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Nx+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; x_i];
        
    %% Boundary conditions: y_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Ny+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; y_i];
    
    %% Boundary conditions: x_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Nx+nx)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; x_f];
    
    %% Boundary conditions: y_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Ny+ny)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; y_f];
    
    %% Boundary conditions: gam_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Ngam+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; gam_i];
    
    %% Boundary conditions: gam_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Ngam+ngam)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; gam_f];
    
    %% Using solver to solve Optimal problem (Without eta)
     [var, fval, exitflag] = quadprog(H, f, A_zone, b_zone, Aeq, beq, lb, ub);
     if exitflag <= 0
        warning('Failed Solver');
        break;
     end

    %% Using solver to solve Optimal problem (Within eta)
%     Var = [x y gam u tf];
%     cvx_begin
%         cvx_solver sedumi
% 
%         variable Solution(nvar)
% %         minimize(0.5*quad_form(Solution, H) )
% 
%         minimize(f'*Solution)
%     
%         subject to
%             Aeq*Solution == beq
%             Solution(Neta+1:Neta+neta)       >= Solution(Nu+1:Nu+nu).^2 
%             A_zone*Solution                  >= b_zone             %NFZ 1
%             Solution(Ns1NFZ+1:Ns1NFZ+ns1NFZ) >= 0
%             Solution(Ns2NFZ+1:Ns2NFZ+ns2NFZ) >= 0
% %             Solution(Ns1x+1:Ns1x+ns1x)       >=0
% %             Solution(Ns2x+1:Ns2x+ns2x)       >=0
% %             Solution(Ns1y+1:Ns1y+ns1y)       >=0
% %             Solution(Ns2y+1:Ns2y+ns2y)       >=0
% %             Solution(Ns1g+1:Ns1g+ns1g)       >=0
% %             Solution(Ns2g+1:Ns2g+ns2g)       >=0
%             Solution(Ntf)                    >= 100
% 
% %             norm(Solution(Nx+1:Nx+nx)-x(1:nx)',2) <= eps            % TRUST REGION
% %             norm(Solution(Ny+1:Ny+ny)-y(1:ny)',2) <= eps             % TRUST REGION
% %             norm(Solution(Ngam+1:Ngam+ngam)-gam(1:ngam)',2) <= 2*pi   % TRUST REGION
%             
%     cvx_end
%     fval = cvx_optval;
%     var = Solution;

     %% Extract solution 
    x   = var(Nx   +1 : Nx   +nx)';
    y   = var(Ny   +1 : Ny   +ny)';
    gam = var(Ngam +1 : Ngam +ngam)';
    u   = var(Nu   +1 : Nu   +nu)';
    eta = var(Neta +1 : Neta +neta)';
    NFZ1 = var(Ns1NFZ +1 : Ns1NFZ  +ns1NFZ)';
    NFZ2 = var(Ns2NFZ +1 : Ns2NFZ  +ns2NFZ)';
    tf  = var(Ntf);
    
    %% store history
    x_hist(:,iter)    = x;
    y_hist(:,iter)    = y;
    gam_hist(:,iter)  = gam;
    u_hist(:,iter)    = u;
    eta_hist(:,iter)  = eta;
    tf_hist(:,iter)   = tf;
    J_hist(iter)      = fval;
    
    %% Calculate Residual 
    Res_x_Vec       = abs((x(2:end) - x(1:end-1))/(dt2) - fx(V, gam(1:end-1)) - fx(V, gam(2:end)));
    Res_y_Vec       = abs((y(2:end) - y(1:end-1))/(dt2) - fy(V, gam(1:end-1)) - fy(V, gam(2:end)));
    Res_gam_Vec     = abs((gam(2:end) - gam(1:end-1))/(dt2) - fgam(V, u(1:end-1)) - fgam(V, u(2:end)));
    Matrix_Res = [Res_x_Vec Res_y_Vec Res_gam_Vec];
    Residual = max(max(Matrix_Res, [], 2));
    fprintf('SCP iter %d residual = %g\n', iter, Residual);
    if Residual < tol
        fprintf('Dynamics satisfied at iteration %d\n',iter);
%         break;
    end
    Res_hist = [Res_hist; Residual];
    fprintf('FINIST AT ITERATION %d\n', iter);
end


t = linspace(t0, tf, N+1);
figure
subplot(3,2,1)
hold on
plot(t,x_hist(:,1)/1e3,'r--','LineWidth',1)
for k=2:iter-1
plot(t,x_hist(:,k)/1e3,'k','LineWidth',0.1)
end
plot(t,x_hist(:,iter-1)/1e3,'g','LineWidth',1)
xlabel('t')
ylabel('km')
title('Horizontal')

subplot(3,2,2)
hold on
plot(t,y_hist(:,1)/1e3,'r--','LineWidth',1)
for k=2:iter-1
plot(t,y_hist(:,k)/1e3,'k','LineWidth',0.1)
end
plot(t,y_hist(:,iter-1)/1e3,'g','LineWidth',1)
xlabel('t')
ylabel('km')
title('Veritical')

subplot(3,2,3)
hold on
plot(t,rad2deg(gam_hist(:,1)),'r--','LineWidth',1)
for k=2:iter-1
plot(t,rad2deg(gam_hist(:,k)),'k','LineWidth',0.1)
end
plot(t,rad2deg(gam_hist(:,iter-1)),'g','LineWidth',1)
xlabel('t')
ylabel('(deg)')
title('Fligth path angle')

subplot(3,2,4)
hold on
plot(t,rad2deg(u_hist(:,1)./g0),'r--','LineWidth',1)
for k=2:iter-1
plot(t,rad2deg(u_hist(:,k)./g0),'k','LineWidth',0.1)
end
plot(t,rad2deg(u_hist(:,iter-1)./g0),'g','LineWidth',1)
xlabel('t')
ylabel('(deg)')
title('Bank angle')

subplot(3,2,5)
plot(1:iter,J_hist,'o-','LineWidth',2)
xlabel('SCP iteration')
ylabel('Objective value')
title('Objective function convergence')

%% -- Trajectory figure -- %%
theta = linspace(0, 2*pi, 200);

x_circle_1 = xc1 + r1*cos(theta);
y_circle_1 = yc1 + r1*sin(theta);

x_circle_2 = xc2 + r2*cos(theta);
y_circle_2 = yc2 + r2*sin(theta);

x_circle_3 = xc3 + r3*cos(theta);
y_circle_3 = yc3 + r3*sin(theta);

x_circle_4 = xc4 + r4*cos(theta);
y_circle_4 = yc4 + r4*sin(theta);

figure
hold on
plot(x_hist(:,1),y_hist(:,1), 'g--', 'LineWidth', 1.5)
for k=2:iter-2
plot(x_hist(:,k),y_hist(:,k),'y')
end
plot(x_hist(:,iter),y_hist(:,iter), 'r', 'LineWidth', 1.5)
hold on
plot(x_circle_1, y_circle_1, 'b--', 'LineWidth', 0.5)
hold on
plot(x_circle_2, y_circle_2, 'b--', 'LineWidth', 0.5)
hold on
plot(x_circle_3, y_circle_3, 'b--', 'LineWidth', 0.5)
hold on
plot(x_circle_4, y_circle_4, 'b--', 'LineWidth', 0.5)
xlabel('km')
ylabel('km')
title('UAV trajectory')
grid on