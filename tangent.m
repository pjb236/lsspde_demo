function tangent(k,inhomo,var)
% solve for tangent (forward sensitivity) v and store to tangent.dat 
%
% INPUTS:
% v_init_[k].dat - file containing initial value of v. 
% k - time segment number
% inhomo - strength of inhomogeneous (forcing) term in tangent equation
% var - index number of design parameter in array "s" from "mesh.dat"
% (can be any index value of array "s" for adjoint lss)
%
% OUTPUTS: 
% v_term_[k].dat - file containing terminal value of v.  


global n s RK varNum dt

varNum = var;

% coefficients for dual consistant order 3 Runge kutta
RK =  [1/2, 0, 0; -1/6, 1./3, 0; 0, -2/3, 1]; 

% read in mesh data from mesh.dat
s = load('mesh.dat');
T = s(end);
dt = s(end-1);
n = s(end-2);
s = s(1:end-3);

% read in time segment data
tChk = load('timeChk.dat');
m = tChk(k,1);

% read in primal solution
u_hist = load(['primal' num2str(k) '.dat']);


% initial condition
u = u_hist(1,:);
f0 = ddt(u,s);
v0 = load(['v_init_' num2str(k) '.dat']);

v = projection(v0,f0); 

% write tangent initial condition to disk
fid = fopen(['tangent' num2str(k) '.dat'],'w');
for j = 1:n-1
    fprintf(fid,'%32.24f ',v(j));
end
fprintf(fid,'%32.24f\n',v(n));


% step in time, write solution to file
for i = 1:m-1
    u = u_hist(i,:);
      
    v = imag(timeStep(u + 1i*1E-60 * v, s)) / 1E-60;

%    v = timeStepTan(v,u,s,inhomo);

    % forcing term
    v = v + dt * inhomo * dfds(u_hist(i+1,:),s);


    for j = 1:n-1
    fprintf(fid,'%32.24f ',v(j));
    end
    fprintf(fid,'%32.24f\n',v(n));
    
end

fclose(fid);
% make terminal condition orthogonal to f(t_{k+1}), then save it to disk.
u = u_hist(end,:);
f = ddt(u,s);
Pvm = projection(v,f);

save2disk(['v_term_' num2str(k) '.dat'],Pvm);
end

function v1 = timeStepTan(v0,u0,s,inhomo)
% make one time step with a 3rd order, dual consistant RK scheme (See Shan
% Yang's thesis)
global RK dt


dudt0 = ddt(u0,s);
u1 = u0 + dt*RK(1,1)*dudt0;
dudt1 = ddt(u1,s);
u2 = u1 + dt*RK(2,1)*dudt0;
u2 = u2 + dt*RK(2,2)*dudt1;


dvdt0 = ddtTan(v0,u0,s);
v1 = v0 + dt*RK(1,1)*dvdt0;
dvdt1 = ddtTan(v1,u1,s);
v1 = v1 + dt*RK(2,1)*dvdt0;
v1 = v1 + dt*RK(2,2)*dvdt1;
dvdt2 = ddtTan(v1,u2,s);
v1 = v1 + dt*RK(3,2)*dvdt1;
v1 = v1 + dt*RK(3,3)*dvdt2;



end

function u1 = timeStep(u0,s)
% make one time step with a 3rd order, dual consistant RK scheme (See Shan
% Yang's thesis)
global RK dt

dudt0 = ddt(u0,s);
u1 = u0 + dt*RK(1,1)*dudt0;
dudt1 = ddt(u1,s);
u1 = u1 + dt*RK(2,1)*dudt0;
u1 = u1 + dt*RK(2,2)*dudt1;
dudt2 = ddt(u1,s);
u1 = u1 + dt*RK(3,2)*dudt1;
u1 = u1 + dt*RK(3,3)*dudt2;

end

function du = ddtTan(v,u,s)
% compute time derivative of u

du = transpose(jacobian(u)*v');

end

function du = ddt(u,s)
% compute time derivative of u using a 2nd order finite difference
% approximation
% approximation
sigma = s(1);
beta = s(2);
r = s(3);
z0 = s(4);

du(1) = sigma*(u(2) - u(1));
du(2) = u(1)*(r - u(3) + z0) - u(2);
du(3) = u(1)*u(2) - beta*(u(3) - z0);

end

function v1 = projection(v0,f)
% project the tangent solution v0(t) into the space orthogonal to f(t)

v1 = v0 - (sum(v0.*f)/sum(f.*f))*f;
end

function du = jacobian(u)
% generate jacobian matrix
global n s eps

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

function du = dfds(u,s)
% find derivative df/ds of f(u) with respect to index varNum of s.  
global varNum

eps = 1E-60;
delta  = zeros(size(s));
delta(varNum)=eps;

du = imag(ddt(u,s+1i*delta))/eps;

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
