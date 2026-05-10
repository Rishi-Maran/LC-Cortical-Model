function [phi_new] = propagator(phi,input,gamma)
% PROPAGATOR Backward Euler Scheme for First Order equation
%   (1 + 1/gamma d/dt)phi(t) = input(t)
%   When discretised:
%       phi_n + (phi_n - phi_(n-1))/(gamma*dt) = input_n
%   Therefore:
%       phi_(n) = (input_n + phi_(n-1)/(gamma*dt)) / (1 + 1/(gamma*dt))
%       = [input - (1/gamma^2-1/gamma)*phi_(n-1) - (1-2/gamma^2)*phi_n]/[1/gamma + 1/gamma^2]
%   List of activity inputs (phi)
	%	phi: Nx x Ny x 2 array containing phi(r, t=(n-1)dt) and phi(r, t = ndt)
%	List of driving inputs (input)
	%	input: Nx x Ny array containing r^2*Laplacian(phi)(r, t = ndt) + C(phi)(r, t = ndt) + nu0*phi(r, t = ndt) + P(r, t = ndt)
%	List of homogeneous parameters (gamma)
	%	gamma = timescale
%	List of outputs (phi_new)
	%	phi_new = Nx x Ny array contianin phi(r, t = (n+1)dt)



phi_new = (input + phi/gamma) / (1 + 1/gamma);


end