function adjoint(k,inhomo)
% solve for adjoint w and store to adjoint.dat
%
% INPUTS:
% w_term_[k].dat - file containing terminal value of w. 
% k - time segment number
% inhomo - strength of inhomogeneous (forcing) term in tangent equation
%
% OUTPUTS: 
% w_init_[k].dat - file containing initial value of w.  

global n RK eps s strength T dt


% coefficients for dual consistant order 3 Runge kutta
RK =  [1/2, 0, 0; -1/6, 1./3, 0; 0, -2/3, 1]; 

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

% read in primal
u_hist = load(['primal' num2str(k) '.dat']);
v_hist = load(['tangent' num2str(k) '.dat']);
Jbar = load('Jbar.dat');




eps = 1i*1E-60;
strength = 1;


% terminal condition
uT = u_hist(end,:);
fT = ddt(uT,s);
wT = load(['w_term_' num2str(k) '.dat']);

w = projection(wT,fT) - inhomo*(1/(T*norm(fT)^2))*(obj(uT) -Jbar)*fT;


% write adjoint terminal condition to disk
fid = fopen(['adjoint' num2str(k) '.dat'],'w');
for j = 1:n-1
    fprintf(fid,'%32.24f ',w(j));
end
fprintf(fid,'%32.24f\n',w(n));

% step in time, write solution to file
for i = 2:m
    u = u_hist(m - i + 1,:);
    v = v_hist(m - i + 1,:);
    
    w = timeStepAdj(w,v,u,s,strength,inhomo);
    for j = 1:n-1
    fprintf(fid,'%32.24f ',w(j));
    end
    fprintf(fid,'%32.24f\n',w(n));

end

fclose(fid);
% make terminal condition orthogonal to f(t_{k}), then save it to disk.
f = ddt(u,s);
Pwp = projection(w,f); 

save2disk(['w_init_' num2str(k) '.dat'],Pwp);


end

function w0 = timeStepAdj(w3,v0,u0,s,strength,inhomo)
% make one time step with a 3rd order, dual consistant RK scheme (See Shan
% Yang's thesis)
global RK T dt

% primal
dudt0 = ddt(u0,s);
u1 = u0 + dt*RK(1,1)*dudt0;
dudt1 = ddt(u1,s);
u2 = u1 + dt*RK(2,1)*dudt0;
u2 = u2 + dt*RK(2,2)*dudt1;

% adjoint 
dwdt_u2w3 = ddtAdj(w3,u2,s);
dwdt_u1w3 = ddtAdj(w3,u1,s);
w2 = w3 - dt * RK(3,3) * dwdt_u2w3;


dwdt_u1w2 = ddtAdj(w2,u1,s);
dwdt_u0w2 = ddtAdj(w2,u0,s);
w2 = w2 - dt * RK(3,2) * dwdt_u1w3;
w1 = w2 - dt * RK(2,2) * dwdt_u1w2;


dwdt_u0w1 = ddtAdj(w1,u0,s);
w1 = w1 - dt * RK(2,1) * dwdt_u0w2;
w0 = w1 - dt * RK(1,1) * dwdt_u0w1;

% adjoint source term

w0 = w0 - dt*strength*(projection(v0,dudt0) - inhomo*(1/T)*ddObj(u0));




end


function dw = ddtAdj(w,u,s)
% compute time derivative of u

dw = transpose(-jacobian(u)'*w');

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

function v1 = projection(v0,f)
% project vector v0 into the space orthogonal to f

v1 = v0 - (sum(v0.*f)/sum(f.*f))*f;
end

function du = jacobian(u)
% generate jacobian matrix
global n s 

sigma = s(1);
beta = s(2);
r = s(3);
z0 = s(4);

% du = zeros(n,n);
% for i = 1:n
%     v = zeros(1,n);
%     v(i) = eps;
%     du(:,i) = real((ddt(u+v,c,nu) - ddt(u,c,nu)) / eps)';
% end

    x = u(1);
    y = u(2);
    z = u(3);
    du = [-sigma, sigma, 0; r - z + z0, -1, -x; y, x, -beta]; 


end

function J = obj(u)
% objective function as a function of the primal solution
    J = u(:,3);
end

function dJdu = ddObj(u)
% derivative of the objective function with respect to the primal solution
    dJdu = [0, 0, 1];
end
function save2disk(filename,v)
% save vector v to disk
global n

fid = fopen(filename,'w');

for j = 1:n-1
    fprintf(fid,'%32.24f ',v(j));
end
fprintf(fid,'%32.24f\n',v(n));

fclose(fid);

end