function ts = run_model(timeparam,localparam,coupparam,connparam,modparam,stim,N_pop,N_nodes)
% RUN with ZERO TRIVIAL CONDITIONS: Q(0) = Q_t(0) = 0;
% Periodic boundary conditions

% Takes in input of model parameters and timesteps.
% Returns timeseries at each point.
% Inputs
%   List of step-size inputs (timeparam)
    %   T: time domain [0, T]
    %   dt: timestep
    %   check if dt divides T
%   List of local model parameters inputs (localparam)
    %   name: name/label of populations a
    %   Q0_a: initial values (not used as of now)
    %   qmax_a: max firing rate (not used as of now)
    %   sigma_a: spread of thresholds
    %   rho_a: slope of activation at steady state
%   List of interpopulation coupling parameters (coupparam)
    %   nu_ab: synaptic coupling parameter
        % first index = postsynaptic population; second index = pre
    %   gamma_e: damping rate for e->e and e->i coupling
    %   tau_ab = transmission delay
    %   alpha, beta = postsynaptic response parameters (double exponential)
%   List of modulation parameters (modparam)
    %   nu_k: synaptic coupling parameters 
        %   first element: from input (phasic), second: from cortex (tonic)
    %   alpha: postsynaptic response parameter for phasic and tonic
    %   rho: slope of activation at steady state
    %   g_a: rate at which population a's theta changes with deltaV
    %   tau_a: transmission delay from LC to population_a
%   List of stimulation inputs (stim)
    %   stimnum: number of times of stimulation
    %   stimtk = times of stimulation - check less than Nt
    %   stimIk = stimulation intensity
    %   sigma = [sigma_x,sigma_t] = Gaussian approximation of delta impulse
% Outputs
%   ts: (N_pop+1) x Nt timeseries array of driver and LC population

twoonpi = 2/pi; piontwo = pi/2;

dt = timeparam.dt;
Nt = timeparam.T / dt;
if Nt ~= round(Nt) 
    error('Nt is not an integer');
end

% % Cannot be both phasic and tonic
% if any(modparam.g > 0) && modparam.Qtonic > 0
%     error('LC cannot be both tonic and phasic');
% end

% Check dimensions of connparam
if any(size(connparam.C) < N_nodes) || any(size(connparam.tau) < N_nodes)
    error("connparam.C and tau must be size %d by %d", N_nodes, N_nodes)
end
conn_C = connparam.C(1:N_nodes, 1:N_nodes);
conn_tau = connparam.tau(1:N_nodes, 1:N_nodes);
% Check entries of connparam.C have unit sum cols and non-neg entries
if ~(all(conn_C(:) >= 0) && all(sum(conn_C, 1) - 1 < 1e-6))
    error("connparam.C must have unit sum cols and non-neg entries");
end
% Check entries of connparam.tau are non-neg entries
if ~all(conn_tau(:) >= 0)
    error("connparam.tau must have non-neg entries");
end

conn_ab = connparam.ab(1:connparam.M, :); 
conn_Chet = connparam.Chet(1:connparam.M);
conn_tauhet = connparam.tauhet(1:connparam.M);
% Check entries of connparam.ab are integers between 1 and N_nodes
if (~isinteger(conn_ab(:)) || ~all(conn_ab(:) <= N_nodes) || ~all(conn_ab(:) > 0))
    error("connparam.ab must be ints between 1 and N_nodes")
end
% Check entires of connparam.Chet are weights between 0 and 1
if any(conn_Chet <= 0) || any(conn_Chet >= 1)
    error("connparam.Chet must be between 0 and 1 exclusive")
end
% Check entries of connparam.tauhet are non-neg entries
if ~all(conn_tauhet(:) >= 0)
    error("connparam.tauhet must have non-neg entries");
end

mod_tau = [modparam.tau_c modparam.tau_c modparam.tau_t modparam.tau_t 0];

% Check modparam.tau is non-neg
if ~all(mod_tau >= 0)
    error("modparam.tau must have non-neg entries");
end

% Check modparam.kernel entries is noneg

mod_k = modparam.kernel(1:N_nodes);
if ~all(mod_k >= 0)
    error("modparam.kernel must have non-neg entries");
end


% Initiate output timeseries array  
ts = zeros(N_nodes + N_nodes + N_pop, Nt);
% Initiate temporary recursive timeseries array
% Need at least tauj previous timesteps for convolution 
% and at least 1 previous timestep for wave equation
% Therefore we track 1 + max(tau, 1) timesteps for each iteration

num_cols_Q = 1 + max([round(conn_tau(:)/dt); 1]);
write_pos_Q = num_cols_Q;
Q = zeros(N_pop, N_nodes, num_cols_Q);

num_cols_prop_e = 1 + max([round(coupparam.tau(:)/dt); 1]);
write_pos_prop_ee = num_cols_prop_e;
prop_e = zeros(N_nodes, num_cols_prop_e);

phi = zeros(N_pop, N_pop, N_nodes);

V = zeros(N_pop, N_nodes, 2);
modfactor = zeros(N_pop, 1);

% LC modulation state
Ql = zeros(1 + max([round(modparam.tau_c/dt); round(modparam.tau_t/dt); 1]), 1);

% Create array of stimulation inputs by time
stimulationinputtime = zeros(Nt, stim.stimnum);
for k = 1:stim.stimnum
    if stim.stimt(k) > timeparam.T
        error("stimulus %d occurs after simulation time", k);
    end
    t0 = stim.stimt(k) / dt;
    for n = 1:Nt
        distt = abs(n - t0);
        stimulationinputtime(n, k) = exp(-0.5 *(distt^2) * (dt^2) / (stim.sigma^2));
    end
end
stimulationinputtime = stimulationinputtime * 1/(sqrt(2*pi) * stim.sigma);
% Normalise stimulationinputspace and stimulationinputtime so that weights
% add to 1/dx^2 and 1/dt respectively
for k = 1:stim.stimnum
    stimulationinputtime(:, k) = (1/dt) * ...
        stimulationinputtime(:, k) / sum(stimulationinputtime(:, k), 'all');
end

% Synaptic sum and product - speed up calculation of synaptic dynamics
S = (1/coupparam.alpha + 1/coupparam.beta)/(2*dt); 
P = (1/coupparam.alpha * 1/coupparam.beta)/(dt^2);

M_S = (1./modparam.alpha + 1./modparam.beta)/(2*dt);
M_P = (1./modparam.alpha * 1./modparam.beta)/(dt^2);

tol = (sigmoid(localparam.V0, localparam.theta, localparam.sigma, localparam.qmax) - localparam.Q0) ./ ...
    localparam.rho;

for count = 1:Nt

    % Compute V from phi
    Vnew = zeros(N_pop, N_nodes);
    for i = 1:N_pop
        input = coupparam.nu(i, :) * squeeze(phi(i, :, :));  % [N_nodes × 1]
        Vnew(i, :) = synaptic(squeeze(V(i, :, :)), input', S, P)';
    end
    % Update
    V(:, :, 1) = V(:, :, 2); V(:, :, 2) = Vnew;

    % Compute Ql from phi
    Ql_new = synaptic(Ql(end-1:end)', modparam.rho*modparam.nu*sum(phi(4, 5, :)), M_S, M_P);
    % Update
    Ql(1:end-1) = Ql(2:end); Ql(end) = Ql_new;

    % Compute modulation factor from Ql
    modfactor_new = propagator(modfactor, ...
        modparam.g(:).*(modparam.Qtonic+Ql(end - round(mod_tau / dt))), modparam.eta*dt);
    % Update
    modfactor = modfactor_new;

    % Compute sigma and theta
    if all(mod_k == 1)  
        sigma = (localparam.sigma).*(1./(1+modfactor));
    else
        sigma = (localparam.sigma).*(1./(1+modfactor.*mod_k));
    end
    theta = localparam.V0 - sigma ./ localparam.sigma .* (localparam.V0 - localparam.theta);
    
    % Compute Q from V, sigma and theta
    Qnew = sigmoid(localparam.V0 + Vnew, theta, sigma, localparam.qmax) - localparam.Q0;
    Qnew(abs(Vnew) <= tol) = 0;
    % Qnew(Vnew == 0) = 0;
    % Update Q and key
    write_pos_Q = write_pos_Q + 1;
    if write_pos_Q > num_cols_Q
        write_pos_Q = 1;
    end
    Q(:,:,write_pos_Q) = Qnew;

    % Compute P_e from Q
    prop_e(:, write_pos_prop_ee) = 0;
    for k = 1:N_nodes
        for l = 1:N_nodes
            pos = write_pos_Q - round(conn_tau(k, l) / dt);
            if pos < 1
                pos = pos + num_cols_Q;
            end
            prop_e(k, write_pos_prop_ee) = prop_e(k, write_pos_prop_ee) + conn_C(k, l)*...
                Q(1, l, pos);
        end
    end

    % Compute phi from Q and P_e
    phinew = zeros(N_pop, N_pop, N_nodes);
    for i = 1:N_pop
        for j = 1:N_pop
            if j == 1 % if excitatory population
                pos = write_pos_prop_ee - round(coupparam.tau(i,j)/dt);
                if pos < 1
                    pos = pos + num_cols_prop_e;
                end
                % I CHANGED THIS!!!!
                phinew(i, j, :) = propagator(...
                    squeeze(phi(i, j, :)), squeeze(prop_e(:, pos)), ...
                    coupparam.gamma * dt);
            else
                pos = write_pos_Q - round(coupparam.tau(i, j) / dt);
                if pos < 1
                    pos = pos + num_cols_Q;
                end
                phinew(i, j, :) = Q(j, :, pos);
            end
        end
    end
    % Update phi and key
    phi = phinew;
    write_pos_prop_ee = write_pos_prop_ee + 1;
    if write_pos_prop_ee > num_cols_prop_e
        write_pos_prop_ee = 1;
    end

    % Update phi_sx with stimulus
    stimulationinput = zeros(N_nodes, 1);
    for k = 1:stim.stimnum
        stimulationinput(stim.pos(k)) = stimulationinput(stim.pos(k)) + ...
            stim.stimI(k) * stimulationinputtime(count, k);
    end
    % Update phi with stimulus
    phi(4, 5, :) = phi(4, 5, :) + reshape(stimulationinput, 1, 1, []);



    ts(:, count) = [reshape(phinew(1, 1, :), [], 1); reshape(Qnew(4, :), [], 1); modfactor_new];

end


end