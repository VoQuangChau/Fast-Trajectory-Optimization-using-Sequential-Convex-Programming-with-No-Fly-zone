clc;
clear;
close all
addpath(genpath(pwd));  

%% Problem "Planetary Entry via sequential convex programming": 
% bank angle control, hypersonic entry starting at the edge of the atmosphere with the optimal control problem of minimizing the impact velocity of a given vehicle

%% Parameter Definitions
global mass Aref R0 g0 V_nor r_nor t_nor rho0 Hs

mass    = 104305;       % Mass of Object                            [kg]
Aref    = 391.22;       % Reference Area of Object                  [m^2]
rho0    = 1.225;        % the sea-level surface Atmospheric density [kg*m^-3]
Hs      = 7200;         % Scale height                              [m]
R0      = 6378e3;       % Earth Radius                              [m]
g0      = 9.8;          % Acceleration at Earth surface             [m*s^-2]
V_nor   = sqrt(R0*g0);  % Normalized velocity                       [m*s^-1]
r_nor   = R0;            % Normalized radius                        [m]
t_nor   = sqrt(R0/g0);  % Normalized time                           [s]
u_nor   = sqrt(g0/R0);  % Normalized bank angle rate                [s^-1]
L_nor   = g0;           % Normalized Lift Acceleration              
D_nor   = g0;           % Normalized Drag Acceleration
status  = 1;
%% Initial and Final of Dynamic Components
r_i   = (R0 + 100e3)./r_nor;    %Dimensionless
r_f   = (R0 + 25e3)./r_nor;     %Dimensionless
%-----------%
V_i   = 7.45e3/V_nor;           %Dimensionless
%-----------%
the_i = deg2rad(0);             %rad
the_f = deg2rad(12);            %rad
%-----------%
phi_i = deg2rad(0);             %rad
phi_f = deg2rad(70);            %rad
%-----------%
gam_i = deg2rad(-0.5);          %rad
gam_f = deg2rad(-10);           %rad
%-----------%
psi_i = deg2rad(0);
psi_f = deg2rad(90);            %rad
%-----------%
sig_i = deg2rad(0);             %rad
%-----------%
u_i   = deg2rad(0);             %rad

%% Boundaries of dynamic components
r_Min   = (R0 + 25e3)./r_nor;    %Dimensionless
r_Max   = (R0 + 100e3)./r_nor;   %Dimensionless
%-----------%
V_Min   = 10.0/V_nor;            %Dimensionless
V_Max   = 7.45e3/V_nor;          %Dimensionless
%-----------%
the_Min = deg2rad(0);            %rad
the_Max = deg2rad(180);          %rad
%-----------%
phi_Min = deg2rad(0);            %rad
phi_Max = deg2rad(90);           %rad
%-----------%
gam_Min = deg2rad(-89);          %rad
gam_Max = deg2rad(89);           %rad
%-----------%
psi_Min = deg2rad(0);
psi_Max = deg2rad(180);          %rad
%-----------%
sig_Min = deg2rad(-80);          %rad
sig_Max = deg2rad(80);           %rad
%-----------%
u_Min   = deg2rad(-10);          %rad
u_Max   = deg2rad(10);           %rad

%% Trust Region
r_trust         = 10000./R0                 ;
the_trust       = 20.*pi./180               ;
phi_trust       = 40.*pi./180               ;
V_trust         = 1000./(sqrt(R0./g0))      ;
gam_trust       = 40.*pi./180               ;
psi_trust       = 60.*pi/180                ;
sig_trust       = 20.*pi/180                ;

%% Discrete trajectory
tf = 1600;                      %seconds
t0 = 0;                         %seconds
N = 100;                         %Nodes

dt      = (tf/t_nor)/(N);
dt2     = dt/2;

xi      = 1e2;

nr      = N+1;  
nthe    = N+1;
nphi    = N+1;
nV      = N+1;
ngam    = N+1;
npsi    = N+1;
nsig    = N+1;
nu      = N+1;
ntf     = 1;
% ntf     = 0;
%% Add compensation term

nslr1      = N;
nslthe1    = N;
nslphi1    = N;
nslV1      = N;
nslgam1    = N;
nslpsi1    = N;
nslsig1    = N;

nslr2      = N;
nslthe2    = N;
nslphi2    = N;
nslV2      = N;
nslgam2    = N;
nslpsi2    = N;
nslsig2    = N;

ns         = nslr1 + nslthe1 + nslphi1 + nslV1 + nslgam1 + nslpsi1 +...
             nslr2 + nslthe2 + nslphi2 + nslV2 + nslgam2 + nslpsi2 ;
nvar       = nr + nthe + nphi + nV + ngam + npsi + nsig + nu + ns + ntf;
          
Nr      = 0;
Nthe    = Nr + nr;
Nphi    = Nthe + nthe;
NV      = Nphi + nphi;
Ngam    = NV + nV;
Npsi    = Ngam + ngam;
Nsig    = Npsi + npsi;
Nu      = Nsig + nsig;

Ns       = Nu + nu;
Nslr1    = Ns;
Nslthe1  = Nslr1   + nslr1;
Nslphi1  = Nslthe1 + nslthe1;
NslV1    = Nslphi1 + nslphi1;
Nslgam1  = NslV1   + nslV1;
Nslpsi1  = Nslgam1 + nslgam1;

Nslr2    = Nslpsi1 + nslpsi1;
Nslthe2  = Nslr2   + nslr2;
Nslphi2  = Nslthe2 + nslthe2;
NslV2    = Nslphi2 + nslphi2;
Nslgam2  = NslV2   + nslV2;
Nslpsi2  = Nslgam2 + nslgam2;

Ntf      = nvar;

max_scp_iter = 40;
tol = 1e-3;

%% Random of Dynamic Components
Var = Initial_Guess_Creation (N+1, dt, r_i, the_i, phi_i, V_i, gam_i, psi_i, sig_i,  r_f, the_f, phi_f, gam_f, psi_f);
% Var = Initial_Guess_Creation_2 (N+1,dt);
r(1:nr)         = Var(Nr+1:Nr+nr);
the(1:nthe)     = Var(Nthe+1:Nthe+nthe);
phi(1:nphi)     = Var(Nphi+1:Nphi+nphi);
V(1:nV)         = Var(NV+1:NV+nV);
gam(1:ngam)     = Var(Ngam+1:Ngam+ngam);
psi(1:npsi)     = Var(Npsi+1:Npsi+npsi);
sig(1:nsig)     = Var(Nsig+1:Nsig+nsig);
u(1:nu)         = Var(Nu+1:Nu+nu);

%% history
r_hist      = zeros(nr,max_scp_iter);
V_hist      = zeros(nV,max_scp_iter);
the_hist    = zeros(nthe,max_scp_iter);
phi_hist    = zeros(nphi,max_scp_iter);
gam_hist    = zeros(ngam, max_scp_iter);
psi_hist    = zeros(npsi, max_scp_iter);
sig_hist    = zeros(nsig, max_scp_iter);
u_hist      = zeros(nu, max_scp_iter);
tf_hist     = zeros(ntf, max_scp_iter);

Res_hist        = [];
J_hist      = zeros(max_scp_iter,1);

%% Store first history of random initial dynamic components
r_hist(:,1)      = r;
V_hist(:,1)      = V;
the_hist(:,1)    = the;
phi_hist(:,1)    = phi;
gam_hist(:,1)    = gam;
psi_hist(:,1)    = psi;
sig_hist(:,1)    = sig;
u_hist(:,1)      = u;
tf_hist(:,1)     = tf;

%% Calculate Residual of Stationary condition for initial guess
Res_r_Vec   = abs((r(2:end) - r(1:end-1))/(dt2)      - gr(V(1:end-1), gam(1:end-1)) - gr(V(2:end), gam(2:end)));
Res_V_Vec   = abs((V(2:end) - V(1:end-1))/(dt2)      - gV(r(1:end-1), V(1:end-1), gam(1:end-1)) - gV(r(2:end), V(2:end), gam(2:end)));
Res_the_Vec = abs((the(2:end) - the(1:end-1))/(dt2)  - gthe(r(1:end-1), phi(1:end-1), V(1:end-1), gam(1:end-1), psi(1:end-1)) - gthe(r(2:end), phi(2:end), V(2:end), gam(2:end), psi(2:end)));
Res_phi_Vec = abs((phi(2:end) - phi(1:end-1))/(dt2)  - gphi(r(1:end-1), V(1:end-1), gam(1:end-1), psi(1:end-1)) - gphi(r(2:end), V(2:end), gam(2:end), psi(2:end)));
Res_gam_Vec = abs((gam(2:end) - gam(1:end-1))/(dt2)  - ggam(r(1:end-1), V(1:end-1), gam(1:end-1), sig(1:end-1)) - ggam(r(2:end), V(2:end), gam(2:end), sig(2:end)));
Res_psi_Vec = abs((psi(2:end) - psi(1:end-1))/(dt2)  - gpsi(r(1:end-1), phi(1:end-1), V(1:end-1), gam(1:end-1), psi(1:end-1), sig(1:end-1)) - gpsi(r(2:end), phi(2:end), V(2:end), gam(2:end), psi(2:end), sig(2:end)));
Res_sig_Vec = abs((sig(2:end) - sig(1:end-1))/(dt2)  - gsig(u(1:end-1)) - gsig(u(2:end)));
Matrix_Res = [Res_r_Vec Res_V_Vec Res_the_Vec Res_phi_Vec Res_gam_Vec Res_psi_Vec Res_sig_Vec];
Residual = max(max(Matrix_Res, [], 2));
Res_hist = [Res_hist; Residual];
fprintf('SCP iter %d residual = %g\n', 1, Residual);

%% Calculate 'Constraint inequality'
lb = -Inf.*ones(nvar, 1);
ub = +Inf.*ones(nvar, 1);

lb(Nr+1:Nr+nr)       = r_Min;       ub(Nr+1:Nr+nr)       = r_Max;
lb(Nthe+1:Nthe+nthe) = the_Min;     ub(Nthe+1:Nthe+nthe) = the_Max;
lb(Nphi+1:Nphi+nphi) = phi_Min;     ub(Nphi+1:Nphi+nphi) = phi_Max;
lb(NV+1:NV+nV)       = V_Min;       ub(NV+1:NV+nV)       = V_Max;
lb(Ngam+1:Ngam+ngam) = gam_Min;     ub(Ngam+1:Ngam+ngam) = gam_Max;
lb(Npsi+1:Npsi+npsi) = psi_Min;     ub(Npsi+1:Npsi+npsi) = psi_Max;
lb(Nsig+1:Nsig+nsig) = sig_Min;     ub(Nsig+1:Nsig+nsig) = sig_Max;
lb(Nu+1:Nu+nu)       = u_Min;       ub(Nu+1:Nu+nu)       = u_Max;

lb(Nslr1+1:Nslr1+nslr1)       = 0;
lb(Nslr2+1:Nslr2+nslr2)       = 0;
lb(Nslphi1+1:Nslphi1+nslphi1) = 0;
lb(Nslphi2+1:Nslphi2+nslphi2) = 0;
lb(Nslthe1+1:Nslthe1+nslthe1) = 0;
lb(Nslthe2+1:Nslthe2+nslthe2) = 0;
lb(NslV1+1:NslV1+nslV1)       = 0;
lb(NslV2+1:NslV2+nslV2)       = 0;
lb(Nslgam1+1:Nslgam1+nslgam1) = 0;
lb(Nslgam2+1:Nslgam2+nslgam2) = 0;
lb(Nslpsi1+1:Nslpsi1+nslpsi1) = 0;
lb(Nslpsi2+1:Nslpsi2+nslpsi2) = 0;


%% Cost function: J = f^T * var
f = zeros(nvar,1); 
% H = zeros(nvar,nvar);
f(nr + nthe + nphi + nV,1) = 1;
% xi = 100.*1.4.^(1-20) + 0.5.*1 + 100;
xi = 1e3;
f(Ns + 1: Ns + ns) = xi;

%% Perform the calculation in each loop (Sequential Programming)
for iter = 2:max_scp_iter %First interation is initial guess for all dynamic components
    
    %% Update penalty weighting
    if iter < 5
        xi = 1e2;
    elseif iter < 10
        xi = 1e4;
    else
        xi = 1e6;
    end
    
    %% Constraint equivality:
    Aeq = [];
    beq = [];
    dt  = tf./t_nor./N;
    dt2 = dt./2;
    for i = 1:N
        %% Dynamic component: gr(V,gam)
        row_r           = zeros(1,nvar);
        
        row_r(Nr + i)        = -(2/dt);
        row_r(Nr + i + 1)    = (2/dt);
        
        row_r(NV + i)       = - dgr_dV(V(i),gam(i));
        row_r(NV + i + 1)   = - dgr_dV(V(i+1),gam(i+1));
        
        row_r(Ngam + i)     = - dgr_dgam(V(i),gam(i));
        row_r(Ngam + i + 1) = - dgr_dgam(V(i+1),gam(i+1)); 
        
        if i>1
            row_r(Nslr1 + i)    = (2/dt);
            row_r(Nslr2 + i)    = -(2/dt);
        end
        
        row_r(Ntf) = - (gr(V(i), gam(i)) + gr(V(i+1), gam(i+1)))./(tf-t0);
        beq_r_0    = 0;
        
%         beq_r_0 = gr(V(i), gam(i)) + gr(V(i+1), gam(i+1));
        beq_r_1 = dgr_dV(V(i), gam(i))*V(i) + dgr_dV(V(i+1), gam(i+1))*V(i+1);
        beq_r_2 = dgr_dgam(V(i), gam(i))*gam(i) + dgr_dgam(V(i+1), gam(i+1))*gam(i+1);
        
        Aeq = [Aeq; row_r];
        beq = [beq; beq_r_0 - beq_r_1 - beq_r_2];
        
        %% Dynamic component: gthe(r,phi,V,gam,psi)
        row_the           = zeros(1,nvar);
        
        row_the(Nthe + i)       = -(2/dt);
        row_the(Nthe + i + 1)   = (2/dt);
        
        row_the(Nr + i)       = - dgthe_dr(r(i),phi(i),V(i),gam(i),psi(i));
        row_the(Nr + i + 1)   = - dgthe_dr(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_the(Nphi + i)     = - dgthe_dphi(r(i),phi(i),V(i),gam(i),psi(i));
        row_the(Nphi + i + 1) = - dgthe_dphi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_the(NV + i)       = - dgthe_dV(r(i),phi(i),V(i),gam(i),psi(i));
        row_the(NV + i + 1)   = - dgthe_dV(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_the(Npsi + i)     = - dgthe_dpsi(r(i),phi(i),V(i),gam(i),psi(i));
        row_the(Npsi + i + 1) = - dgthe_dpsi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_the(Ngam + i)       = - dgthe_dgam(r(i),phi(i),V(i),gam(i),psi(i));
        row_the(Ngam + i + 1)   = - dgthe_dgam(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1));

        if i>1
            row_the(Nslthe1 + i)    = (2/dt);
            row_the(Nslthe2 + i)    = -(2/dt);
        end
        
        row_the(Ntf) = - (gthe(r(i),phi(i),V(i),gam(i),psi(i)) + gthe(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1)))./(tf-t0);
        beq_the_0    = 0;
        
%         beq_the_0 = gthe(r(i),phi(i),V(i),gam(i),psi(i)) + gthe(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1));
        beq_the_1 = dgthe_dr(r(i),phi(i),V(i),gam(i),psi(i))*r(i) + dgthe_dr(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1))*r(i+1);
        beq_the_2 = dgthe_dphi(r(i),phi(i),V(i),gam(i),psi(i))*phi(i) + dgthe_dphi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1))*phi(i+1);
        beq_the_3 = dgthe_dV(r(i),phi(i),V(i),gam(i),psi(i))*V(i) + dgthe_dV(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1))*V(i+1);
        beq_the_4 = dgthe_dgam(r(i),phi(i),V(i),gam(i),psi(i))*gam(i) + dgthe_dgam(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1))*gam(i+1);
        beq_the_5 = dgthe_dpsi(r(i),phi(i),V(i),gam(i),psi(i))*psi(i) + dgthe_dpsi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1))*psi(i+1);
        
        Aeq = [Aeq; row_the];
        beq = [beq; beq_the_0 - beq_the_1 - beq_the_2 - beq_the_3 - beq_the_4 - beq_the_5];
        
        %% Dynamic component: gphi(r,V,gam,psi)
        row_phi           = zeros(1,nvar);
        
        row_phi(Nphi + i)       = -(2/dt);
        row_phi(Nphi + i + 1)   = (2/dt);
        
        row_phi(Nr + i)         = - dgphi_dr(r(i),V(i),gam(i),psi(i));
        row_phi(Nr + i + 1)     = - dgphi_dr(r(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_phi(NV + i)         = - dgphi_dV(r(i),V(i),gam(i),psi(i));
        row_phi(NV + i + 1)     = - dgphi_dV(r(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_phi(Ngam + i)       = - dgphi_dgam(r(i),V(i),gam(i),psi(i));
        row_phi(Ngam + i + 1)   = - dgphi_dgam(r(i+1),V(i+1),gam(i+1),psi(i+1));
        
        row_phi(Npsi + i)       = - dgphi_dpsi(r(i),V(i),gam(i),psi(i));
        row_phi(Npsi + i + 1)   = - dgphi_dpsi(r(i+1),V(i+1),gam(i+1),psi(i+1));

        if i>1
            row_phi(Nslphi1 + i)    = (2/dt);
            row_phi(Nslphi2 + i)    = -(2/dt);
        end
%         
        row_phi(Ntf) = - (gphi(r(i),V(i),gam(i),psi(i)) + gphi(r(i+1),V(i+1),gam(i+1),psi(i+1)))./(tf-t0);
        beq_phi_0    = 0;
        
%         beq_phi_0 = gphi(r(i),V(i),gam(i),psi(i)) + gphi(r(i+1),V(i+1),gam(i+1),psi(i+1));
        beq_phi_1 = dgphi_dr(r(i),V(i),gam(i),psi(i))*r(i) + dgphi_dr(r(i+1),V(i+1),gam(i+1),psi(i+1))*r(i+1);
        beq_phi_2 = dgphi_dV(r(i),V(i),gam(i),psi(i))*V(i) + dgphi_dV(r(i+1),V(i+1),gam(i+1),psi(i+1))*V(i+1);
        beq_phi_3 = dgphi_dgam(r(i),V(i),gam(i),psi(i))*gam(i) + dgphi_dgam(r(i+1),V(i+1),gam(i+1),psi(i+1))*gam(i+1);
        beq_phi_4 = dgphi_dpsi(r(i),V(i),gam(i),psi(i))*psi(i) + dgphi_dpsi(r(i+1),V(i+1),gam(i+1),psi(i+1))*psi(i+1);
        
        Aeq = [Aeq; row_phi];
        beq = [beq; beq_phi_0 - beq_phi_1 - beq_phi_2 - beq_phi_3 - beq_phi_4];
        
        %% Dynamic component: gV(r,V,gam)
        row_V       = zeros(1,nvar);
        
        row_V(NV + i)       = -(2/dt) - dgV_dV(r(i),V(i),gam(i));
        row_V(NV + i + 1)   = (2/dt) - dgV_dV(r(i+1),V(i+1),gam(i+1));
        
        row_V(Nr + i)       = - dgV_dr(r(i),V(i),gam(i));
        row_V(Nr + i + 1)   = - dgV_dr(r(i+1),V(i+1),gam(i+1));
        
        row_V(Ngam + i)     = - dgV_dgam(r(i),V(i),gam(i));
        row_V(Ngam + i + 1) = - dgV_dgam(r(i+1),V(i+1),gam(i+1));

        if i>1
            row_V(NslV1 + i)    = (2/dt);
            row_V(NslV2 + i)    = -(2/dt);
        end
        
        row_V(Ntf) = - (gV(r(i),V(i),gam(i)) + gV(r(i+1),V(i+1),gam(i+1)))./(tf-t0);
        beq_V_0    = 0;
        
%         beq_V_0 = gV(r(i),V(i),gam(i)) + gV(r(i+1),V(i+1),gam(i+1));
        beq_V_1 = dgV_dr(r(i),V(i),gam(i))*r(i) + dgV_dr(r(i+1),V(i+1),gam(i+1))*r(i+1);
        beq_V_2 = dgV_dV(r(i),V(i),gam(i))*V(i) + dgV_dV(r(i+1),V(i+1),gam(i+1))*V(i+1);
        beq_V_3 = dgV_dgam(r(i),V(i),gam(i))*gam(i) + dgV_dgam(r(i+1),V(i+1),gam(i+1))*gam(i+1);
        
        Aeq = [Aeq; row_V];
        beq = [beq; beq_V_0 - beq_V_1 - beq_V_2 - beq_V_3];
        
        %% Dynamic component: ggam(r, V, gam, sig)
        row_gam       = zeros(1,nvar);
        
        row_gam(Ngam + i)       = -(2/dt) - dggam_dgam(r(i),V(i),gam(i),sig(i));
        row_gam(Ngam + i + 1)   = (2/dt) - dggam_dgam(r(i+1),V(i+1),gam(i+1),sig(i+1));
        
        row_gam(Nr + i)       = - dggam_dr(r(i),V(i),gam(i),sig(i));
        row_gam(Nr + i + 1)   = - dggam_dr(r(i+1),V(i+1),gam(i+1),sig(i+1));
        
        row_gam(NV + i)       = - dggam_dV(r(i),V(i),gam(i),sig(i));
        row_gam(NV + i + 1)   = - dggam_dV(r(i+1),V(i+1),gam(i+1),sig(i+1));
        
        row_gam(Nsig + i)     = - dggam_dsig(r(i),V(i),gam(i),sig(i));
        row_gam(Nsig + i + 1) = - dggam_dsig(r(i+1),V(i+1),gam(i+1),sig(i+1));
        
        
        
        if i>1
            row_gam(Nslgam1 + i)    = (2/dt);
            row_gam(Nslgam2 + i)    = -(2/dt);
        end
        
        row_gam(Ntf) = - (ggam(r(i),V(i),gam(i),sig(i)) + ggam(r(i+1),V(i+1),gam(i+1),sig(i+1)))./(tf-t0);
        beq_gam_0    = 0;
        
%         beq_gam_0 = ggam(r(i),V(i),gam(i),sig(i)) + ggam(r(i+1),V(i+1),gam(i+1),sig(i+1));
        beq_gam_1 = dggam_dr(r(i),V(i),gam(i),sig(i))*r(i) + dggam_dr(r(i+1),V(i+1),gam(i+1),sig(i+1))*r(i+1);
        beq_gam_2 = dggam_dV(r(i),V(i),gam(i),sig(i))*V(i) + dggam_dV(r(i+1),V(i+1),gam(i+1),sig(i+1))*V(i+1);
        beq_gam_3 = dggam_dgam(r(i),V(i),gam(i),sig(i))*gam(i) + dggam_dgam(r(i+1),V(i+1),gam(i+1),sig(i+1))*gam(i+1);
        beq_gam_4 = dggam_dsig(r(i),V(i),gam(i),sig(i))*sig(i) + dggam_dsig(r(i+1),V(i+1),gam(i+1),sig(i+1))*sig(i+1);
        
        Aeq = [Aeq; row_gam];
        beq = [beq; beq_gam_0 - beq_gam_1 - beq_gam_2 - beq_gam_3 - beq_gam_4];
        %% Dynamic component: gpsi(r, phi, V, gam, psi, sig)
        row_psi       = zeros(1,nvar);
        
        row_psi(Npsi + i)       = -(2/dt) - dgpsi_dpsi(r(i), phi(i), V(i),gam(i), psi(i), sig(i));
        row_psi(Npsi + i + 1)   = (2/dt) - dgpsi_dpsi(r(i+1), phi(i+1), V(i+1),gam(i+1), psi(i+1), sig(i+1));
        
        row_psi(Nr + i)         = - dgpsi_dr(r(i), phi(i), V(i),gam(i), psi(i), sig(i));
        row_psi(Nr + i + 1)     = - dgpsi_dr(r(i+1), phi(i+1), V(i+1),gam(i+1), psi(i+1), sig(i+1));
        
        row_psi(Nphi + i)       = - dgpsi_dphi(r(i), phi(i), V(i),gam(i), psi(i), sig(i));
        row_psi(Nphi + i + 1)   = - dgpsi_dphi(r(i+1), phi(i+1), V(i+1),gam(i+1), psi(i+1), sig(i+1));
        
        row_psi(NV + i)         = - dgpsi_dV(r(i), phi(i), V(i),gam(i), psi(i), sig(i));
        row_psi(NV + i + 1)     = - dgpsi_dV(r(i+1), phi(i+1), V(i+1),gam(i+1), psi(i+1), sig(i+1));
        
        row_psi(Ngam + i)       = - dgpsi_dgam(r(i), phi(i), V(i),gam(i), psi(i), sig(i));
        row_psi(Ngam + i + 1)   = - dgpsi_dgam(r(i+1), phi(i+1), V(i+1),gam(i+1), psi(i+1), sig(i+1));
        
        row_psi(Nsig + i)       = - dgpsi_dsig(r(i), phi(i), V(i),gam(i), psi(i), sig(i));
        row_psi(Nsig + i + 1)   = - dgpsi_dsig(r(i+1), phi(i+1), V(i+1),gam(i+1), psi(i+1), sig(i+1));

        if i>1
            row_psi(Nslpsi1 + i)    = (2/dt);
            row_psi(Nslpsi2 + i)    = -(2/dt);
        end
        
        row_psi(Ntf) = - (gpsi(r(i),phi(i),V(i),gam(i),psi(i),sig(i)) + gpsi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1)))./(tf-t0);
        beq_psi_0    = 0;
        
%         beq_psi_0 = gpsi(r(i),phi(i),V(i),gam(i),psi(i),sig(i)) + gpsi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1));
        beq_psi_1 = dgpsi_dr(r(i),phi(i),V(i),gam(i),psi(i),sig(i))*r(i) + dgpsi_dr(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1))*r(i+1);
        beq_psi_2 = dgpsi_dphi(r(i),phi(i),V(i),gam(i),psi(i),sig(i))*phi(i) + dgpsi_dphi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1))*phi(i+1);
        beq_psi_3 = dgpsi_dV(r(i),phi(i),V(i),gam(i),psi(i),sig(i))*V(i) + dgpsi_dV(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1))*V(i+1);
        beq_psi_4 = dgpsi_dgam(r(i),phi(i),V(i),gam(i),psi(i),sig(i))*gam(i) + dgpsi_dgam(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1))*gam(i+1);
        beq_psi_5 = dgpsi_dpsi(r(i),phi(i),V(i),gam(i),psi(i),sig(i))*psi(i) + dgpsi_dpsi(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1))*psi(i+1);
        beq_psi_6 = dgpsi_dsig(r(i),phi(i),V(i),gam(i),psi(i),sig(i))*sig(i) + dgpsi_dsig(r(i+1),phi(i+1),V(i+1),gam(i+1),psi(i+1),sig(i+1))*sig(i+1);
        
        Aeq = [Aeq; row_psi];
        beq = [beq; beq_psi_0 - beq_psi_1 - beq_psi_2 - beq_psi_3 - beq_psi_4 - beq_psi_5 - beq_psi_6];
        
        %% Dynamic component: gsig(u)
        row_sig       = zeros(1,nvar);
        
        row_sig(Nsig + i)       = -(2/dt);
        row_sig(Nsig + i + 1)   = (2/dt);
        
        row_sig(Nu + i)         = -dgsig_du(u(i));
        row_sig(Nu + i + 1)     = -dgsig_du(u(i+1));
        
%         row_sig(Ntf) = - (gsig(u(i)) + gsig(u(i+1)))./(tf-t0);
%         beq_sig_0    = 0;

        beq_sig_0   = gsig(u(i)) + gsig(u(i+1));
        beq_sig_1   = dgsig_du(u(i))*u(i) + dgsig_du(u(i+1))*u(i+1);
        
        Aeq = [Aeq; row_sig];
        beq = [beq; beq_sig_0 - beq_sig_1];
    end
    
    %% Boundary conditions: r_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Nr+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; r_i];
        
    %% Boundary conditions: the_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Nthe+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; the_i];
    
    %% Boundary conditions: phi_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Nphi+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; phi_i];
    
    %% Boundary conditions: V_i
    row_bc1 = zeros(1,nvar);
    row_bc1(NV+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; V_i];
    
    %% Boundary conditions: gam_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Ngam+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; gam_i];
    
    %% Boundary conditions: psi_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Npsi+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; psi_i];
    
    %% Boundary conditions: sig_i 
    row_bc1 = zeros(1,nvar);
    row_bc1(Nsig+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; sig_i];
    
    %% Boundary conditions: u_i
    row_bc1 = zeros(1,nvar);
    row_bc1(Nu+1)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; u_i];   

    %% Boundary conditions: r_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Nr + nr)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; r_f];
        
    %% Boundary conditions: the_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Nthe + nthe)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; the_f];
    
    %% Boundary conditions: phi_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Nphi + nphi)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; phi_f];
    
    %% Boundary conditions: gam_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Ngam + ngam) = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; gam_f];
    
    %% Boundary conditions: psi_f
    row_bc1 = zeros(1,nvar);
    row_bc1(Npsi + npsi)  = 1;
    Aeq = [Aeq; row_bc1];
    beq = [beq; psi_f];

    %% Trust region update (Update boundary)
%     if iter > 20
%         r_trust         = r_trust     .*2;
%         the_trust       = the_trust   .*2;
%         phi_trust       = phi_trust   .*2;
%         V_trust         = V_trust     .*2;
%         gam_trust       = gam_trust   .*2;
%         psi_trust       = psi_trust   .*2;
%         sig_trust       = sig_trust   .*2;
%     end
    if iter > 2
        lb(Nr+1:Nr+nr)       = max(r_Min,r(1:nr)-r_trust);       ub(Nr+1:Nr+nr)       = min(r_Max,r(1:nr)+r_trust);
        lb(Nthe+1:Nthe+nthe) = max(the_Min,the(1:nr)-the_trust); ub(Nthe+1:Nthe+nthe) = min(the_Max,the(1:nr)+the_trust);
        lb(Nphi+1:Nphi+nphi) = max(phi_Min,phi(1:nr)-phi_trust); ub(Nphi+1:Nphi+nphi) = min(phi_Max,phi(1:nr)+phi_trust);
        lb(NV+1:NV+nV)       = max(V_Min,V(1:nr)-V_trust)      ; ub(NV+1:NV+nV)       = min(V_Max,V(1:nr)+V_trust);
        lb(Ngam+1:Ngam+ngam) = max(gam_Min,gam(1:nr)-gam_trust); ub(Ngam+1:Ngam+ngam) = min(gam_Max,gam(1:nr)+gam_trust);
        lb(Npsi+1:Npsi+npsi) = max(psi_Min,psi(1:nr)-psi_trust); ub(Npsi+1:Npsi+npsi) = min(psi_Max,psi(1:nr)+psi_trust);
        lb(Nsig+1:Nsig+nsig) = max(sig_Min,sig(1:nr)-sig_trust); ub(Nsig+1:Nsig+nsig) = min(sig_Max,sig(1:nr)+sig_trust);
        lb(Nu+1:Nu+nu)       = u_Min;       ub(Nu+1:Nu+nu)       = u_Max;   
    end
    
    %% Using solver to solve Optimal problem
     [var, fval, exitflag] = quadprog([], f, [], [], Aeq, beq, lb, ub);
     if exitflag <= 0
        warning('Failed Solver');
        break;
     end

    %% Extract solution 
    
    r   = var(Nr   +1 : Nr   +nr)';
    the = var(Nthe +1 : Nthe +nthe)';
    phi = var(Nphi +1 : Nphi +nphi  )';
    V   = var(NV   +1 : NV   +nV)';
    gam = var(Ngam +1 : Ngam +ngam)';
    psi = var(Npsi +1 : Npsi +npsi)';
    sig = var(Nsig +1 : Nsig +nsig)';
    u   = var(Nu   +1 : Nu   +nu)';
    slr1   = var(Nslr1   +1:Nslr1   + nslr1)';
    slthe1 = var(Nslthe1 +1:Nslthe1 + nslthe1)';
    slphi1 = var(Nslphi1 +1:Nslphi1 + nslphi1)';
    slV1   = var(NslV1   +1:NslV1   + nslV1)';
    slgam1 = var(Nslgam1 +1:Nslgam1 + nslgam1)';
    slpsi1 = var(Nslpsi1 +1:Nslpsi1 + nslpsi1)';
    
	slr2   = var(Nslr2   +1:Nslr2   + nslr2)';
    slthe2 = var(Nslthe2 +1:Nslthe2 + nslthe2)';
    slphi2 = var(Nslphi2 +1:Nslphi2 + nslphi2)';
    slV2   = var(NslV2   +1:NslV2   + nslV2)';
    slgam2 = var(Nslgam2 +1:Nslgam2 + nslgam2)';
    slpsi2 = var(Nslpsi2 +1:Nslpsi2 + nslpsi2)';
    
    tf     = var(Ntf);
    %% store history
    r_hist(:,iter) = r;
    V_hist(:,iter) = V;
    the_hist(:,iter) = the;
    phi_hist(:,iter) = phi;
    psi_hist(:,iter) = psi;
    gam_hist(:,iter) = gam;
    sig_hist(:,iter) = sig;
    u_hist(:,iter)   = u;
    tf_hist(:,iter)  = tf;
    J_hist(iter) = fval;
    
    %% Calculate Residual 
    Res_r_Vec   = abs((r(2:end) - r(1:end-1))/(dt/2)      - gr(V(1:end-1), gam(1:end-1)) - gr(V(2:end), gam(2:end)));
    Res_V_Vec   = abs((V(2:end) - V(1:end-1))/(dt/2)      - gV(r(1:end-1), V(1:end-1), gam(1:end-1)) - gV(r(2:end), V(2:end), gam(2:end)));
    Res_the_Vec = abs((the(2:end) - the(1:end-1))/(dt/2)  - gthe(r(1:end-1), phi(1:end-1), V(1:end-1), gam(1:end-1), psi(1:end-1)) - gthe(r(2:end), phi(2:end), V(2:end), gam(2:end), psi(2:end)));
    Res_phi_Vec = abs((phi(2:end) - phi(1:end-1))/(dt/2)  - gphi(r(1:end-1), V(1:end-1), gam(1:end-1), psi(1:end-1)) - gphi(r(2:end), V(2:end), gam(2:end), psi(2:end)));
    Res_gam_Vec = abs((gam(2:end) - gam(1:end-1))/(dt/2)  - ggam(r(1:end-1), V(1:end-1), gam(1:end-1), sig(1:end-1)) - ggam(r(2:end), V(2:end), gam(2:end), sig(2:end)));
    Res_psi_Vec = abs((psi(2:end) - psi(1:end-1))/(dt/2)  - gpsi(r(1:end-1), phi(1:end-1), V(1:end-1), gam(1:end-1), psi(1:end-1), sig(1:end-1)) - gpsi(r(2:end), phi(2:end), V(2:end), gam(2:end), psi(2:end), sig(2:end)));
    Res_sig_Vec = abs((sig(2:end) - sig(1:end-1))/(dt/2)  - gsig(u(1:end-1)) - gsig(u(2:end)));
    Matrix_Res = [Res_r_Vec Res_V_Vec Res_the_Vec Res_phi_Vec Res_gam_Vec Res_psi_Vec Res_sig_Vec];
    Residual = max(max(Matrix_Res, [], 2));
    fprintf('SCP iter %d residual = %g\n', iter, Residual);
    Res_hist = [Res_hist; Residual];
    if Residual < tol
        fprintf('Dynamics satisfied at iteration %d\n',iter);
        break;
    end   
end

%% plot results

t = linspace(0,tf,N+1);

figure

subplot(5,2,1)
hold on
plot(t,r_hist(:,1)*r_nor,'r--','LineWidth',1)
for k=2:iter-1
plot(t,r_hist(:,k)*r_nor,'k','LineWidth',0.1)
end
plot(t,r_hist(:,iter-1)*r_nor,'g','LineWidth',1)
xlabel('t')
ylabel('m')
title('Altitude')

subplot(5,2,2)
hold on
plot(t,phi_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,phi_hist(:,k),'k','LineWidth',0.1)
end
plot(t,phi_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('rad')
title('Longitude')

subplot(5,2,3)
hold on
plot(t,the_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,the_hist(:,k),'k','LineWidth',0.1)
end
plot(t,the_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('m/s')
title('Latitude')

subplot(5,2,3)
hold on
plot(t,the_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,the_hist(:,k),'k','LineWidth',0.1)
end
plot(t,the_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('m/s')
title('Latitude')

subplot(5,2,4)
hold on
plot(t,V_hist(:,1)*V_nor,'r--','LineWidth',1)
for k=2:iter-1
plot(t,V_hist(:,k)*V_nor,'k','LineWidth',0.1)
end
plot(t,V_hist(:,iter-1)*V_nor,'g','LineWidth',1)
xlabel('t')
ylabel('m/s')
title('Velocity')

subplot(5,2,5)
hold on
plot(t,gam_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,gam_hist(:,k),'k','LineWidth',0.5)
end
plot(t,gam_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('rad')
title('Flight Path Angle')

subplot(5,2,6)
hold on
plot(t,psi_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,psi_hist(:,k),'k','LineWidth',0.1)
end
plot(t,psi_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('rad')
title('Heading Angle')

subplot(5,2,7)
hold on
plot(t,sig_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,sig_hist(:,k),'k','LineWidth',0.1)
end
plot(t,sig_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('rad')
title('Bank angle')

subplot(5,2,8)
hold on
plot(t,u_hist(:,1),'r--','LineWidth',1)
for k=2:iter-1
plot(t,u_hist(:,k),'k','LineWidth',0.1)
end
plot(t,u_hist(:,iter-1),'g','LineWidth',1)
xlabel('t')
ylabel('rad')
title('Control')

figure
subplot(1,2,1)
plot(1:iter,J_hist(1:iter),'o-','LineWidth',2)
xlabel('SCP iteration')
ylabel('Objective value')
title('Objective function convergence')
grid on

subplot(1,2,2)
plot(1:iter,Res_hist(1:iter),'o-','LineWidth',2)
xlabel('SCP iteration')
ylabel('Residual value')
title('Residual convergence')
grid on

figure
subplot(1,2,1)
plot(rad2deg(the(:)),rad2deg(phi(:)),'r','LineWidth',2)
xlabel('[deg]')
ylabel('[deg]')
title('Trajectory')
grid on

subplot(1,2,2)
plot(V*V_nor/1000,(r*r_nor-R0)/1000,'LineWidth',2)
xlabel('[km/s]')
ylabel('[km]')
title('Trajectory')
grid on