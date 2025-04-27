clear; clc; close all
addpath WFSim/libraries/sparse_null
addpath WFSim/bin/core                      % WFSim model directory
addpath bin/core                            % WFAMPC directory
%addpath bin\archive                        % WFAMPC directory (old needs to be replaced)

Wp.name             = 'ThreeTurbine_Ampc';
Wp.Turbulencemodel  = 'WFSim3';

%% Init
WFAMPC_initialize
%{
% WFSim settings in WFAMPC
Animate       = 0;                   % Show 2D flow fields
conv_eps      = 1e-6;                % Convergence threshold
max_it_dyn    = 1;                   % Maximum number of iterations for k > 1

options.Projection     = 0;                      % Use projection (true/false)
options.Linearversion  = 1;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 1;                      % Calculate linear solution of WFSim
options.Derivatives    = 1;                      % Compute derivatives
options.exportPressures= ~options.Projection;    % Calculate pressure fields

[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,0);

%% function InitWFSim
function [Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,plotMesh)

sys    = struct;
sol    = struct;

Projection    = options.Projection;
Linearversion = options.Linearversion;

% Create meshing and import control settings
[Wp,input]   = meshing(Wp,plotMesh,1); 

% Initial flow fields
[sol.u,sol.uu] = deal(Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));  
[sol.v,sol.vv] = deal(Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));  
[sol.p,sol.pp] = deal(Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny)); 


if Linearversion
    sol.ul = sol.u;
    sol.vl = sol.v;
    sol.pl = sol.p;
    [sol.du,sol.dv,sol.dp]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
end

% Initialize parameters are empty matrices
[Power,CT,Ueffect,a] = deal(zeros(Wp.turbine.N,Wp.sim.NN)); 

% Compute boundary conditions and matrices B1, B2
[B1,B2,bc]           = Compute_B1_B2_bc(Wp);
B2                   = 2*B2;

end


% MPC settings in WFAMPC
options.AMPC.ShowGrad     = 0;                   % Show plot of gradient for each simulation
options.AMPC.ShowBeta     = 0;                   % Show plot of beta for each simulation
options.AMPC.method       = 'grad_ratio';        % Control method used

options.AMPC.beta_lim        = [0.1;0.9];        % Limits of input beta
options.AMPC.dbeta_max       = 0.1;              % Limit on change of beta per time step
options.AMPC.Np              = 400;              % Prediction Horizon (in time steps)
options.AMPC.Nc              = 400;              % Control horizon (in time steps)
options.AMPC.Nr              = 20;               % Receding horizon (in time steps)
options.AMPC.gamma           = 5e-7;             % For Steepest-descent method
options.AMPC.imax            = 10;               % Max number of iterations within timestep
options.AMPC.iter_ls         = 1;
options.AMPC.pen             = 1e4*ones(1,options.AMPC.Nc);   % penalty on dbeta
options.AMPC.labda_end       = 100;              % end-point penalty
options.AMPC.eps_grad        = 0.01;
options.AMPC.dP_norm         = 0.002;
options.AMPC.Nrmax           = 30;               % Number of receding horizons
options.AMPC.counter         = 0;
%%
%}

%% Simulate Wind Farm towards Steady State
index                   = 0;
%options.startUniform    = 1;    % Start from a uniform flowfield (true) or a steady-state solution (false)
%max_it                  = 50;  
%%
%RunWFSim;                          % index 1 (go to steady state)
%{
index                   = index + 1;

for k=1:options.AMPC.Np
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    
    while ( eps>conv_eps && it<max_it && eps<epss )
        it   = it+1;
        epss = eps;
        
        if k>1
            max_it = max_it_dyn;
        end
        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);                      % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                      % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);            % Map solution to field
        % Normv{k} = norm(vec(sol.v(2:end-1,3:end-1)-sol.vv(2:end-1,3:end-1)));
        % Normu{k} = norm(vec(sol.u(3:end-1,2:end-1)-sol.uu(3:end-1,2:end-1)));
        % eps          = sqrt((Normv{k}+Normu{k}))/((Ny-2)*(Nx-2))/2;
        % where sol.v represents the previous velocity and sol.vv
        % represents the current velocity.
        
    end
    
    derivatives = sys.derivatives;
    x           = sol.x;
    save((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(k))),'sys','x');    
end
%}

%% First forward simulation + Backwards adjoint
%load((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(options.AMPC.Nr))))

options.startUniform    = 1;    % Start from a uniform flowfield (true) or a steady-state solution (false)
max_it                  = 1;  

RunWFSim                        % index 2 (start from steady state)
Power_greedy = Power;
%load((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(options.AMPC.Nr))))
%sol.u

%%% until now, the input is still greedy control policy, but the flow field
%%% becomes to the steady state



[grad,grad_Phi,labda,J_partial]   = get_adjoint_Nc(options.AMPC,x,Wp,index); %%%% In this function, only the length of x is used 

% This is temporarily, time and inputs need to be redefined in meshing
for kk=1:options.AMPC.Nc  
    beta(:,kk) = input{kk}.beta;
    Phi(:,kk) = input{kk}.phi;
end

grad 	                 = gradproj(grad,beta,options.AMPC.beta_lim,eps);
grad_Phi = gradproj(grad_Phi,Phi,options.AMPC.Phi_lim,eps);
%betaup      = (beta' > beta_lim(2) - eps) & (grad < 0);
%betadown    = (beta' < beta_lim(1) + eps) & (grad > 0);
%grad        = grad.*(1-betaup).*(1-betadown);

%dJmax                    = max(max(abs(grad)))';

%% Gradient calculation (Can this init not be done above, at least some can)
stop_ls                     = 0;
%grad                        = get_adjoint_Nc(options.AMPC,x,Wp,index); %%%%%%% repeat with above??????
dJmaxinit                   = max(max(abs(grad)))'; %%%%%%%% repeat with above????????
alphai                      = 1/(dJmaxinit);
alphai0                     = alphai;
BETA(:,1:options.AMPC.Nr)   = beta(:,1:options.AMPC.Nr);
GRAD(:,1:options.AMPC.Nr)   = grad(1:options.AMPC.Nr,:)';

dJmaxinit_phi                   = max(max(abs(grad_Phi)))'; %%%%%%%% repeat with above????????
alphai_phi                      = 1/(dJmaxinit_phi);
alphai0_phi                     = alphai_phi;
PHI(:,1:options.AMPC.Nr)   = Phi(:,1:options.AMPC.Nr);
GRAD_PHI(:,1:options.AMPC.Nr)   = grad_Phi(1:options.AMPC.Nr,:)';

POWER(:,1:options.AMPC.Nr)  = Power(:,1:options.AMPC.Nr); %%%%% this power will be constant, because it is start from steady state and the control policy is greedy control;
POWERTOT(:,1)               = sum(Power,2);
J(1)                        = sum(POWERTOT(:,1));
P                           = zeros(options.AMPC.Nrmax,options.AMPC.imax+1);
P2                          = zeros(options.AMPC.Nrmax,options.AMPC.imax+2);
ALPHA                       = [];
dP_normi                    = options.AMPC.dP_norm*ones(1,options.AMPC.Nrmax);
constant                    = 0;

index = 1;
    
for i = 1:options.AMPC.Nrmax
    
    % Forward simulation 1 -> Np
    load((strcat('data_WFAMpc/states/state',num2str(index),'_',num2str(options.AMPC.Nr))))
    if constant == 0
        beta        = [beta(:,options.AMPC.Nr+1:end) beta(:,end)*ones(1,options.AMPC.Nr)]; %%%%%% Because after one linear search, the first Nr columns have been implemented; in the next line search process, the initial control variable shouldn't include the first Nr time steps. 
        beta0       = beta(:,options.AMPC.Nr); %%%%%% only one column, not 1:Nr columns           when we solve the finite-horizon optimization problem at each iteration, we need a new initial control variable starting from current time step. 
    end
    
    % Write here an update of input.beta
    for kk=1:size(beta,2)
        input{kk}.beta = beta(:,kk);
    end
    %%% at this time, beta is still the greedy control policy
    constant    = 0;
    indexNr     = index;
    RunWFSim
    
    % Backward adjoint Np -> 1
    grad        = get_adjoint_Nc(options.AMPC,x,Wp,index);
%     grad        = gradproj(grad,beta,beta_lim,eps);
    beta_prev   = beta;
    Power_prev  = Power;
    ls          = 1;
    check1      = 0;
    check2      = 0;
    alphai      = alphai*8; % Een keer *8 proberen?
    P(i,ls)     = sum(sum(Power));
    P2(i,ls)    = P(i,ls);
    display(['i =',num2str(i),', Initial cost J =',num2str(P(i,ls),'%10.4e')])
    
    P2(i,ls+1)  = P(i,ls);
    
    if i > 1
        dP(i)       = abs((P(i,1)/P(i-1,1))-1);
    else
        dP(1)       = dP_normi(1) + 1;
    end
    
    while ls <= options.AMPC.imax && ~(check1 && check2) && dP(i) > dP_normi(i)
        
        beta_prev2  = beta;
        Power_prev2 = Power;
        check1      = 0;
        check2      = 0;
             
        % Determine new beta
%         if P(i,ls) <= P(i,1)
%             alphai  = alphai/2;
%             ALPHA(length(ALPHA)+1) = alphai;
%         elseif P2(i,ls+1) > P2(i,ls)
%             alphai  = alphai*2;
%             ALPHA(length(ALPHA)+1) = alphai;
%         end
        alphai  = alphai/2;
        %ALPHA(length(ALPHA)+1) = alphai;
        
%         beta    = update_beta_abs(beta_prev,grad,beta_lim,Nc,dbeta_max,alphai);
        beta    = update_beta(beta_prev,grad,alphai,beta0,options);
        %beta(:,1)
        
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        load((strcat('data_WFAMpc/states/state',num2str(indexNr),'_',num2str(options.AMPC.Nr))))
        RunWFSim
        ls          = ls + 1;
        P(i,ls)     = sum(sum(Power));
        P2(i,ls+1)  = P(i,ls);
        display(['i =',num2str(i),', line search ',num2str(ls-1),', Cost J =',num2str(P(i,ls),'%10.4e')])
        
        if P(i,ls) > P(i,1) || P(i,ls-1) > P(i,1)
            check1  = 1;
        end
        if P(i,ls) < P(i,ls-1)
            check2  = 1;
        end
        
    end
    
    if dP(i) < dP_normi(i)
%         beta        = beta(:,1)*ones(1,option.AMPC.Np);
        beta(:,1:options.AMPC.Nr)    = beta(:,1)*ones(1,options.AMPC.Nr);
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        %         dP_normi(i+1)   = dP_norm*2;
        load((strcat('data_WFAMpc/states/state',num2str(indexNr),'_',num2str(options.AMPC.Nr))))
        RunWFSim
        P(i,2)      = sum(sum(Power));
        P2(i,3)     = P(i,2);
        constant    = 1;
        alphai      = alphai0;
        display(['i =',num2str(i),', constant Power, Cost J =',num2str(P(i,1),'%10.4e')]);
    elseif P(i,ls) < P(i,ls-1) && P(i,ls-1) > P(i,1)
        disp('after #ls linear search, the power start decreasing');
        beta    = beta_prev2;
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        Power   = Power_prev2;
        index   = index-1;
        stopls  = 1;
    elseif P(i,ls) < P(i,1)
        disp('after 10 linear search, it cannot find a larger power');
        beta    = beta_prev;
        % Write here an update of input.beta
        for kk=1:size(beta,2)
            input{kk}.beta = beta(:,kk);
        end
        
        Power   = Power_prev;
        index   = indexNr;
    end
    
    % Save intermediate results
    BETA(:,1+options.AMPC.Nr*(i-1):options.AMPC.Nr*i)   = beta(:,1:options.AMPC.Nr);
    GRAD(:,1+options.AMPC.Nr*(i-1):options.AMPC.Nr*i)   = grad(1:options.AMPC.Nr,:)';
    POWER(:,1+options.AMPC.Nr*(i-1):options.AMPC.Nr*i)  = Power(:,1:options.AMPC.Nr);
    POWERTOT(:,i)                                       = sum(Power,2);
    J(i)                                                = sum(POWERTOT(:,i));
    
    if options.AMPC.ShowGrad == 1
        figure;
        plot(grad);
        figure;
        plot(beta');
    end
    
end
    
%% Figures: opmaak verbeteren!
h           = Wp.sim.h;
Nrmax       = options.AMPC.Nrmax;
Nr          = options.AMPC.Nr;
Np = options.AMPC.Np;
Power_greedy_tot = sum(Power_greedy);
POWER_tot = sum(POWER);
(POWER_tot(end) - Power_greedy_tot(end))/Power_greedy_tot(end)
figure;hold on;plot(1:Np,sum(Power_greedy))
plot(h:h:h*Nrmax*Nr,sum(POWER))
grid on
xlabel('Time [s]')
ylabel('Power [W]')
legend('Greedy Control','AMPC')

figure;hold on;plot(h:h:h*Nrmax*Nr,POWER')
plot(h:h:h*Nrmax*Nr,sum(POWER))


grid on
xlabel('Time [s]')
ylabel('Power [W]')
legend('Turbine 1','Turbine 2','Turbine 3','Total power')

figure;hold on;plot(h:h:h*Nrmax*Nr,BETA')
plot(h:h:h*Nrmax*Nr,[0.26;0.1;0.54]*ones(1,Nrmax*Nr),'--')
xlabel('Time [s]')
ylabel('\beta [-]')
grid on
legend('Turbine 1','Turbine 2','Turbine 3','Optimum T_1','Optimum T_2','Optimum T_3')

figure; plot(h:h:h*Nrmax*Nr,GRAD); hold on
plot(h:h:h*Nrmax*Nr,zeros(1,Nrmax*Nr),'--')
xlabel('Time [s]')
grid on
ylabel('\nabla_{\beta} J')

figure; plot(h:h:h*Nrmax,POWERTOT');hold on
plot(h:h:h*Nrmax,J)
xlabel('Time [s]')
ylabel('Power [W]')
legend('Turbine 1','Turbine 2','Turbine 3','Total power')
grid on

delete('data_WFAMpc/states/*.mat')
%}
