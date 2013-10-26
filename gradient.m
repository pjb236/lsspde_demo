function [dJds] = gradient(k)


% read in mesh data from mesh.dat
s = load('mesh.dat');
T = s(end);
dt = s(end-1);
n = s(end-2);
s = s(1:end-3);

% read in time chuck data
tChk = load('timeChk.dat');
[K tmp] = size(tChk);
m = tChk(k,1);

% read in primal, tangent and time averaged objecive function
u_hist = load(['primal' num2str(k) '.dat']);
v_hist = load(['tangent' num2str(k) '.dat']);
Jbar = load('Jbar.dat');


term1 = dt*sum(v_hist(1:end-1,3));


f = ddt(u_hist(end,:),s);
term2 = (dot(v_hist(end,:),f)/dot(f,f))*(obj(u_hist(end,:))-Jbar);

dJds = (1/T)*(term1- term2);

end

function J = obj(u)
    J = u(:,3);
end

function du = ddt(u,s)
% compute time derivative of u using a 2nd order finite difference
% approximation
sigma = s(1);
beta = s(2);
r = s(3);
z0 = s(4);

du(1) = sigma*(u(2) - u(1));
du(2) = u(1)*(r - u(3) + z0) - u(2);
du(3) = u(1)*u(2) - beta*(u(3) - z0);

end


