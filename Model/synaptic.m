function [V_new] = synaptic(V,input,S,P)
%WAVE_EQ Simulate the wave equation of model using finite difference scheme
%   Wave equation can be written as:
%       Q_n 
%           + (1/alpha + 1/beta)*(Q_(n+1)-Q_(n-1))/(2*dt) 
%           + (1/(alpha*beta))*(Q_(n+1)+Q_(n-1)-2*Q_n)/(dt^2)
%       = input = sum_b rho_a nu_ab phi_ab
%   Therefore:
%       phi_(n+1)
%       = [input - (P-1/gamma)*phi_(n-1) - (1-2/gamma^2)*phi_n]/[1/gamma + 1/gamma^2]
%   List of activity inputs (phi)
	%	phi: Nx x Ny x 2 array containing phi(r, t=(n-1)dt) and phi(r, t = ndt)
%	List of driving inputs (input)
	%	input: Nx x Ny array containing r^2*Laplacian(phi)(r, t = ndt) + C(phi)(r, t = ndt) + nu0*phi(r, t = ndt) + P(r, t = ndt)
%	List of homogeneous parameters (gamma)
	%	gamma = damping rate
%	List of outputs (phi_new)
	%	phi_new = Nx x Ny array contianin phi(r, t = (n+1)dt)

    V_new = (input - (P - S)*V(:, end-1) - (1 - 2*P)*V(:, end)) / (S + P);
    
end