%% Set time parameters

timeparam = struct;
timeparam.T = 2;
timeparam.dt = 2^(-14);

%% Set number of populations

N_pop = int32(5);

%% Set inter-population couplings and parameters

% array of nus;
nu = zeros(N_pop, N_pop);

% Eyes open parameters from Abeysuriya, Rennie, Robinson 2015
nu(1, 1) = 0.00785; %ee
nu(1, 2) = -0.00988; %ei
nu(1, 4) = 0.00090; %es
nu(2, 1) = 0.00785; %ie
nu(2, 2) = -0.00988; %ii
nu(2, 4) = 0.00090; %is
nu(3, 1) = 0.00021; %re
nu(3, 4) = 6e-05; %rs
nu(4, 1) = 0.00268; %se
nu(4, 3) = -0.00131; %sr
nu(4, 5) = 0.00660; %sn

gamma = 116;

tau = zeros(N_pop, N_pop);
tau(1, 4) = 4.248046875e-02; %es
tau(2, 4) = 4.248046875e-02; %is
tau(3, 1) = 4.248046875e-02; %re
tau(4, 1) = 4.248046875e-02; %se

coupparam.nu = nu;
coupparam.gamma = gamma;
coupparam.tau = tau;

% Dendrite parameters
coupparam.alpha = 83.33333333;
coupparam.beta = 769.2307692;

clear nu gamma tau i j;

%% Set populations and parameters

localparam.qmax = [340; 340; 340; 340; 340];

localparam.sigma = [0.0038; 0.0038; 0.0038; 0.0038; 0.0038];

localparam.theta = [0.0129; 0.0129; 0.0129; 0.0129; 0.0129];

% Solve for steady states
localparam.Q0 = Q_root(coupparam.nu,localparam.theta(1),localparam.sigma(1),localparam.qmax(1));

localparam.V0 = logsigmoid(localparam.Q0, localparam.theta, localparam.sigma, localparam.qmax);

localparam.rho = (localparam.Q0)./(localparam.sigma) .* ...
    (1 - localparam.Q0./localparam.qmax);

clear i localparam.qmax localparam.sigma;

%% Set corticocortical connectivity for multiple populations

N_nodes = int32(1);
connparam.C = zeros(200, 200);
connparam.tau = zeros(200, 200);

connparam.M = int32(0);
connparam.ab = int32(zeros(500, 2));
connparam.Chet = zeros(500, 1);
connparam.tauhet = zeros(500, 1);

%% Set global modulation population and parameters

%phasic coupling
modparam.nu = coupparam.nu(4, 5);
modparam.Qtonic = 0;
modparam.eta = 20; % set timeconstant as 50ms from Vilardaga 2016, 2003

% dendrite parameters
modparam.alpha = coupparam.alpha;
modparam.beta = coupparam.beta;
modparam.rho = localparam.rho(4);
modparam.g = [0; 0; 0; 0; 0];

% timedelays
modparam.tau_c = 4.248046875e-02;
modparam.tau_t = 0;

% For cholinergic system
modparam.kernel = ones(1, 200);

%% Set stimulation parameters

stim = struct;

stim.stimnum = int32(1);

stim.pos = zeros(100, 1);
stim.pos(1) = 101;

stim.sigma = timeparam.dt;

stimt = zeros(100,1);
stimt(1) = 0.1;
stimI = zeros(100,1);
stimI(1) = 1;

stim.stimt = stimt;
stim.stimI = stimI;

stim.noise = 0;

stim.seed = 1;


clear stimt stimI;

%% Set number of nodes and distance-dependent connectivity matrix

N_nodes = 200;

% Set lambda as that fitted with tractography data
load('Data/fitting.mat', 'lambda'); lambda = lambda*1000;

% Set diameter as longest measured interregional euclidean distance
load('Data/schaefer400_left_euclidean_dist.mat', 'dist_p');
diameter = max(dist_p, [], 'all') / 1000; clear dist_p;

% Set axonal velocity as 2m/s by Lyle Muller's review
v_hom = 2;

for i = 1:N_nodes
    for j = 1:N_nodes
        angle = min(abs(i-j), N_nodes - abs(i-j)) / N_nodes * 2*pi;
        dist = diameter*sin(angle/2);
        connparam.C(i, j) = exp(-lambda*dist);
        connparam.tau(i, j) = dist/v_hom;
    end
end

connparam.C = connparam.C ./ sum(connparam.C, 1);

%% Create a single LRC of diameter length

lrc_a = 51; lrc_b = 151;

connparam_het = connparam;
connparam_het.C(lrc_a, lrc_b) = connparam.C(lrc_a, lrc_b)*1e5; 
connparam_het.C(lrc_b, lrc_a) = connparam.C(lrc_b, lrc_a)*1e5;
connparam_het.C = connparam_het.C ./ sum(connparam_het.C, 1);
connparam_het.tau(lrc_a, lrc_b) = 0; connparam_het.tau(lrc_b, lrc_a) = 0;


%% Convert N_nodes to int32

N_nodes = int32(N_nodes);



