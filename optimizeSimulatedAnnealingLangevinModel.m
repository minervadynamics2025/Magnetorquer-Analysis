function [ p_best ] = optimizeSimulatedAnnealingLangevinModel( N_d, B, H_0, H_init )
% this function finds optimal M_s and a values for Langeving Model
% B    : Measured Magnetic Field
% H_init : Initial guess for Magnetizing Field
% N_d : Demagnetization factor

mu0 = 1.257*1e-6;         % Air permeability [H/m]

% Here we combine H_0 and H_init in a single vector
H = zeros(length(H_0),2);
H(:,1) = H_0;
H(:,2) = H_init;

%% Langevin Model function for Magnetic Field in mT
model = @(p, H)  mu0 * (H(:,1) + p(1) * (coth((H(:,2))/p(2)) - p(2)./(H(:,2))));

%% Cost function (least squares)
cost = @(p) sum(((B/1000 - model(p,H))).^2);

%% SA Parameters
maxIter = 5000;
T = 1.0;           % Initial temperature
T_min = 1e-6;
alpha = 0.95;      % Cooling rate

% Parameter bounds
lb = [4500, 500];
ub = [20000, 5000];

% Initial solution
p_current = lb + rand(1,2).*(ub - lb);
cost_current = cost(p_current);

p_best = p_current;
cost_best = cost_current;

%% SA Loop
for k = 1:maxIter
    
    % Generate new candidate (random perturbation)
    step_size = 0.1*(ub - lb);
    p_new = p_current + step_size.*randn(size(p_current));
    
    % Enforce bounds
    p_new = max(p_new, lb);
    p_new = min(p_new, ub);
    
    % Here we update Effective Magnetizing Field H_f for new parameters
    for i=1:length(H_0)
        H(i,2) = getHfieldLangevinModel( N_d, H_0(i), p_new(1), p_new(2) );
    end 
    %% Cost function (least squares)
    cost = @(p) sum(((B/1000. - model(p,H))).^2);
    
    cost_new = cost(p_new);
    
    % Acceptance criterion
    if cost_new < cost_current
        accept = true;
    else
        prob = exp(-(cost_new - cost_current)/T);
        accept = rand < prob;
    end
    
    if accept
        p_current = p_new;
        cost_current = cost_new;
    end
    
    % Update best solution
    if cost_current < cost_best
        p_best = p_current;
        cost_best = cost_current;
    end
    
    % Cooling
    T = T * alpha;
    
    % Stopping condition
    if T < T_min
        break;
    end
end

end