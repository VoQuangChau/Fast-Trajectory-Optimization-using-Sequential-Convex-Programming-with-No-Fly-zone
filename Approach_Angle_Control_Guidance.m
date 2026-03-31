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

%% Boundaries of dynamic components
u_max = +deg2rad(80);
u_min = -deg2rad(80);

%% Trust region defined
eps = 1;

%% Split Node
N   =   100; %Node
tf  =   sqrt(x_f^2+y_f^2)/V; %Seconds
t0  =   0; %Seconds
dt  =   tf / N;
dt2 =   dt/2;

nx      = N+1;  
ny      = N+1;
ngam    = N+1;
nu      = N+1;
neta    = N+1;
% nvar    = nx + ny + ngam + nu + neta + 1;
nvar    = nx + ny + ngam + nu + 1;

Nx      = 0;
Ny      = Nx + nx;
Ngam    = Ny + ny;
Nu      = Ngam + ngam;
% Neta    = Nu + nu;
% Ntf     = Neta + neta + 1;
Ntf     = Nu + nu + 1;

max_scp_iter = 20;
tol = 1e-4;
%% Random of Dynamic Components (Linear Guess)
x       = linspace(x_i,x_f,nx)';
y       = linspace(y_i,y_f,ny)';
gam     = linspace(gam_i,gam_f,ngam)';
u       = zeros(nu,1);
eta     = zeros(neta,1);

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

% lb(Nu+1:Nu+nu)       = u_min;       ub(Nu+1:Nu+nu)       = u_max;

%% Cost function: J = f^T * var
% f = zeros(nvar,1); 
% f(Neta+1:Neta+neta,1) = 1;
H = zeros(nvar,nvar);
for i = 1:nu
   H(Nu+i,Nu+i) = 2*dt;
end

%% Perform the calculation in each loop (Sequential Programming)
for iter = 2:max_scp_iter %First interation is initial guess for all dynamic components
%     Var =   [r the phi V gam psi sig u eps_r eps_the eps_phi eps_V eps_gam eps_psi eps_sig];
    dt      = tf / N;
    dt2     = dt / 2;
    
    %% Constraint equivality:
    Aeq = [];
    beq = [];
    
    for i = 1:N
        %% Dynamic component: fx(V,gam)
        row_x      = zeros(1,nvar);
        row_x(Nx + i)       = - (1./dt2);
        row_x(Nx + i + 1)   = + (1./dt2);
        row_x(Ngam + i)     = - dfx_dgam(V,gam(i));
        row_x(Ngam + i + 1) = - dfx_dgam(V,gam(i+1));
        row_x(Ntf)          = - (fx(V,gam(i)) + fx(V,gam(i+1)))./(tf-t0);
        
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
        
        Aeq = [Aeq; row_gam];
        beq = [beq; 0];
    end
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
    
%     % Objective function
%     obj = @(x) 0.5 * x' * H * x + f' * x ;
% %     obj = @(x) f(1:Neps_r)'*x(1:Neps_r) + f(Neps_r+1:nvar)'*abs(x(Neps_r+1:nvar));
% 
%     % Nonlinear constraint (SOC)
%     nonlcon = @(z) soc_constraint(z, Var', r_trust, the_trust, phi_trust, V_trust, gam_trust, psi_trust, sig_trust, Nr, Nthe, Nphi, NV, Ngam, Npsi, Nsig, N+1);
% 
%     % Options
%     options = optimoptions('fmincon', ...
%         'Display','iter', ...
%         'Algorithm','interior-point', ...
%         'MaxIterations',100,...
%         'MaxFunctionEvaluations', 1e5);
% 
%     % Solve
%     [Solution, fval, exitflag, output] = fmincon( ...
%         obj, x0, ...
%         [], [], ...           % no A, b (inequality linear)
%         Aeq, beq, ...         % equality
%         lb, ub, ...           % bounds
%         nonlcon, ...          % nonlinear (SOC)
%         options);
    %% Using solver to solve Optimal problem
     [var, fval, exitflag] = quadprog(H, [], [], [], Aeq, beq, lb, ub);
     if exitflag <= 0
        warning('Failed Solver');
        break;
     end
     %% Extract solution 
%     var = Solution;
    x   = var(Nx   +1 : Nx   +nx)';
    y   = var(Ny   +1 : Ny   +ny)';
    gam = var(Ngam +1 : Ngam +ngam)';
    u   = var(Nu   +1 : Nu   +nu)';
%     eta = var(Neta +1 : Nta  +neta)';
    tf  = var(end);

    
    %% store history
    x_hist(:,iter)    = x;
    y_hist(:,iter)    = y;
    gam_hist(:,iter)  = gam;
    u_hist(:,iter)    = u;
%     eta_hist(:,iter)  = eta;
    tf_hist(:,iter)   = tf;
    J_hist(iter) = fval;
    
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
plot(t,rad2deg(gam_hist(:,1)),'r--','LineWidth',1)
for k=2:iter-1
plot(t,rad2deg(gam_hist(:,k)),'k','LineWidth',0.1)
end
plot(t,rad2deg(gam_hist(:,iter-1)),'g','LineWidth',1)
xlabel('t')
ylabel('(deg)')
title('Fligth path angle')

subplot(3,2,5)
plot(1:max_scp_iter,J_hist,'o-','LineWidth',2)
xlabel('SCP iteration')
ylabel('Objective value')
title('Objective function convergence')

figure
plot(x,y,'r','LineWidth',2)
xlabel('km')
ylabel('km')
title('UAV trajectory')
grid on