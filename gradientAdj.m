function [dJds] = gradientAdj(var)
% Compute gradient using the adjoint solution.  
global varNum
varNum = var;

% read in necessary data from disk
tChk = load('timeChk.dat');
Jbar = load('Jbar.dat');
s = load('mesh.dat');
T = s(end);
dt = s(end-1);
n = s(end-2);
s = s(1:end-3);
% compute time step size, length of overall time interval
[K tmp] = size(tChk);

dJds = 0;

for k = 1:K
    % load primal and adjoint
    u = load(['primal' num2str(k) '.dat']);
    w = load(['adjoint' num2str(k) '.dat']);
    
    % find inner product of adjoint and df/ds
    m = tChk(k,1);
    for i = 2:m
         dJds = dJds + dt * w(m-i+1,:) * dfds(u(i,:),s)';
    end

    
end

end

function du = dfds(u,s)
% find derivative df/ds of f(u) with respect to index varNum of s.  
global varNum

eps = 1E-60;
delta  = zeros(size(s));
delta(varNum)=eps;

du = imag(ddt(u,s+1i*delta))/eps;

end

function du = ddt(u,s)
% compute time derivative of u 
sigma = s(1);
beta = s(2);
r = s(3);
z0 = s(4);

du(1) = sigma*(u(2) - u(1));
du(2) = u(1)*(r - u(3) + z0) - u(2);
du(3) = u(1)*u(2) - beta*(u(3) - z0);

end