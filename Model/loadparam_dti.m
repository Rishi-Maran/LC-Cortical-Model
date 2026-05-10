%% Set time parameters

timeparam = struct;
timeparam.T = 2;
timeparam.dt = 2^(-14);

%% Set number of populations

N_pop = int32(5);

%% Set inter-population couplings and parameters

% array of nus;
nu = zeros(N_pop, N_pop);
nu(1, 1) = 0.00785; %ee
nu(1, 2) = -0.00988; %ei
nu(1, 4) = 0.00090; %es
nu(2, 1) = 0.00785; %ie
nu(2, 2) = -0.00988; %ii
nu(2, 4) = 0.00090; %is
nu(3, 1) = 0.00021; %re
nu(3, 4) = 0.00006; %rs
nu(4, 1) = 0.00268; %se
nu(4, 3) = -0.00131; %sr
nu(4, 5) = 0.00660; %sn

gamma = 116;

tau = zeros(N_pop, N_pop);
tau(1, 4) = 4.248046875e-02; %es
tau(2, 4) = 4.248046875e-02; %is
tau(3, 1) = 4.248046875e-02; %re
tau(4, 1) = 4.248046875e-02; %se

% for i = 1:N_pop
    % for j = 1:N_pop
        coupparam.nu = nu;
        coupparam.gamma = gamma;
        coupparam.tau = tau;
    % end
% end

% Dendrite parameters
coupparam.alpha = 83.33333333;
coupparam.beta = 769.2307692;

clear nu gamma tau i j;

%% Set populations and parameters

localparam.qmax = [340; 340; 340; 340; 340];

localparam.sigma = [0.0038; 0.0038; 0.0038; 0.0038; 0.0038];

localparam.theta = [0.0129; 0.0129; 0.0129; 0.0129; 0.0129];

% localparam.Q0 = [5.16; 5.16; 16.35; 8.4; 1];
localparam.Q0 = Q_root(coupparam.nu,localparam.theta(1),localparam.sigma(1),localparam.qmax(1));

localparam.V0 = logsigmoid(localparam.Q0, localparam.theta, localparam.sigma, localparam.qmax);

localparam.rho = (localparam.Q0)./(localparam.sigma) .* ...
    (1 - localparam.Q0./localparam.qmax);

% localparam.qmax(1:2) = 7.10704 - 5.24456; localparam.qmax(3:end) = NaN;

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
modparam.nu = 0;
modparam.Qtonic = 0;
modparam.eta = 25;

% dendrite parameters
modparam.alpha = coupparam.alpha;
modparam.beta = coupparam.beta;
modparam.rho = localparam.rho(4);
modparam.g = 0.01*[1; 0; 0; 0; 0];

% timedelays
modparam.tau_c = 4.248046875e-02;
modparam.tau_t = 0;

% For cholinergic system
modparam.kernel = ones(1, 200);

%% Set stimulation parameters

stim = struct;

stim.stimnum = int32(1);

stim.pos = zeros(100, 1);

% Determine parcel of most posterior vertex of mesh
v_data = gifti('../Model/Data/fs_LR.32k.L.midthickness.surf.gii').vertices;
[~, v_idx] = min(v_data(:, 2)); clear v_data;
V_P_map = dlmread('../Model/Data/fsLR_32k_Schaefer400-lh.txt');
stim.pos(1) = V_P_map(v_idx); clear V_P_map;

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

load('Data/fitting.mat', 'lambda');
load('Data/schaefer400_left_euclidean_dist.mat', 'dist_p');
lambda = lambda*1000; dist_p = dist_p / 1000;

N_nodes = 200;

% Set diameter of circle as what is from pang paper
diameter = 0.067*2;

connparam.C = exp(-dist_p * lambda);
connparam.tau = dist_p / 2;
clear dist_p;

connparam.C = connparam.C ./ sum(connparam.C, 1);

%% Create a single LRC of diameter length
% Create a single LRC of length 0.09 according to sporns rich club paper

load('../Model/Data/struct_mat_p.mat');
load('../Model/Data/outliers.mat');

N_lrcs = size(outliers, 1)*2;
connparam.M = int32(N_lrcs);
connparam.ab(1:N_lrcs, :) = ...
    [outliers(:, [1 2]);
    outliers(:, [2 1])];
for i = 1:size(outliers, 1)
    connparam.Chet(i) = sc_p(outliers(i, 1), outliers(i, 2)) / sum(sc_p(outliers(i, 1), :));
    connparam.Chet(size(outliers, 1)+i) = sc_p(outliers(i, 2), outliers(i, 1)) / sum(sc_p(outliers(i, 2), :));
end

clear outliers sc_p;

%% Convert N_nodes to int32

N_nodes = int32(N_nodes);



