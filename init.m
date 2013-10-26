function init
% initialize primal solution u

% variables
global T0 T1 dt m RK n

%%%%% mesh data
sigma = 10; 
beta = 8/3;
r = 28.0; % rayleigh number
z0 = 0.0; % z0
s = [sigma beta r z0];
n = 3; % number of nodes -> s and n -> mesh.dat

%%%%%

T0 = 40; % Run up time
T1 = 10; % run time
dt = 0.01; % time step size
m = round(T1/dt); % number of time step
% coefficients for dual consistant order 3 Runge kutta
RK =  [1/2, 0, 0; -1/6, 1./3, 0; 0, -2/3, 1]; 

fid = fopen('mesh.dat','w');
fprintf(fid,'%f %f %f %f %d %f %f\n',sigma,beta,r,z0,n,dt,T1);
fclose(fid);


% intial condition
u0 = zeros([1,n]); %1*(rand([1,n])- 0.5);
u0(1) = 0.01;
u0(2) = 0.01;

% run from initial condition to T0
u0 = spinup(u0,s);

% compute primal, save to disk
primal(u0,s);

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



function u = spinup(u0,s)
% run solution up to time T0
global T0 dt

m0 = round(T0/dt);

u = u0;
for i = 1:m0
    u = timeStep(u,s);
end

end

function primal(u0,s)
% solve for u and store to primal.dat
global m n dt


fid_u = fopen('primal1.dat','w');
fid_t = fopen('timeChk.dat','w');
% fid_v = fopen('tangent.dat','w');

u = u0;
eps = 1E-60;

% write initial conditions to file
v = ones(size(u0))/norm(ones(size(u0)));
v = imag(timeStep(u + 1i*eps * v, s)) / eps;

for j = 1:n-1
    fprintf(fid_u,'%32.24f ',u(j));
end
fprintf(fid_u,'%32.24f\n',u(n));

% initialize primal solver
i0 = 0;
k = 1;
Jbar = 0;
nrms = [];

% cutOffRatio gives the maximum ratio between the norm of the tangent and 
% primal solutions in each time chunk.    

cutOffRatio = 15; 

% step in time, write solution to file
for i = 1:m
    % determine time chunk length 
    nrms = [nrms; norm(v)];
    
    % end current time segment, start new one if cutOffRatio is exceeded.  
    if norm(v) > cutOffRatio 
        fprintf(fid_t,'%d %f\n',i-i0,(i-i0)*dt);

        i0 = i-1;
        k = k + 1;
        v = ones(size(u0))/norm(ones(size(u0)));
        fclose(fid_u);
        fid_u = fopen(['primal' num2str(k) '.dat'],'w');  
            % store primal on disk
            for j = 1:n-1
            fprintf(fid_u,'%32.24f ',u(j));
            % fprintf(fid_v,'%f ',v(j));

            end
            fprintf(fid_u,'%32.24f\n',u(n));
            % fprintf(fid_v,'%f\n',v(n));
    end    
    
    % step forward in time
    u = timeStep(u,s);
    
    v = imag(timeStep(u + 1i*eps * v, s)) / eps;
        
    % save primal to disk
    for j = 1:n-1
    fprintf(fid_u,'%32.24f ',u(j));
    end
    fprintf(fid_u,'%32.24f\n',u(n));
    
    % update objective function average
    Jbar = Jbar + (1/m)*obj(u);
end

% store info on final time chunk
fprintf(fid_t,'%d %f\n',m-i0+1,(m-i0)*dt);


fid_J = fopen('Jbar.dat','w');
fprintf(fid_J,'%32.24f\n',Jbar);
fclose(fid_J);

fclose(fid_u);
fclose(fid_t);

plot(nrms)

end

function J = obj(u)
% objective function as a function of the primal solution
    J = u(:,3);
end
