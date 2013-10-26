function lsspde_tan
% run LSS algorithm.  Only run once init has been executed!

global K n

% detemine time chunk data
tChk = dlmread('timeChk.dat', ' ');
[K tmp] = size(tChk);

% read in mesh data from mesh.dat
s = dlmread('mesh.dat', ' ');
T = s(end);
dt = s(end-1);
n = s(end-2);
s = s(1:end-3);

% input the index of the design variable s for which to compute sensitivity
var = 4;

% run lsspde algorithm, using the following programs:
% Tangent –n [time chunk number] –i v_initial.dat –o v_terminal.dat
% Adjoint –n [time chunk number] –i w_terminal.dat –o w_initial.dat

% solve tangent for each chunk
x0 = zeros(n*(2*K-1),1);
rhs = matvec(x0,-1,var);

% solve with minres
[x1, flag, relres, iter, resvec] =  minres(@(x)matvec(x,0,var),rhs,1e-4,20);

disp(flag)
disp(resvec)

figure(4);
semilogy(resvec)
xlabel('iteration','FontSize',16)
ylabel('||r||','FontSize',16)
set(gca,'FontSize',16)

% % % exact solve
% % % generate matrix
% % A = zeros((2*K-1)*n,(2*K-1)*n);
% % 
% % for i = 1:(2*K-1)*n
% %     e1 = zeros((2*K-1)*n,1);
% %     e1(i) = 1;
% %     A(:,i) = matvec(e1,0,0,var);
% %     
% % end
% % 
% % x1 = A\rhs;

% form solution with initial & terminal conditions in x1
residual = matvec(x1,1,var);

% compute gradient
dJds = 0;
for k = 1:K
    dJds = dJds + gradient(k);
end
fprintf('%f\n',dJds);


plotter;

end

function y = matvec(x,inhomoT,var)
global K n

% adjoint boundary conditions
w0 = zeros(1,n);
wTn = zeros(1,n);

% solve tangent for each chunk
v0 = reshape(x(1:n*K),n,K)';
wT = reshape(x(n*K+1:end),n,K-1)';
Pvm = zeros(K,n);
Pwp = zeros(K,n);

for k = 1:K
    %%% tangent solver
    save2disk(['v_init_' num2str(k) '.dat'],v0(k,:));
    tangent(k,inhomoT,var);
    Pvm(k,:) = load(['v_term_' num2str(k) '.dat']);
    %%%
    %%% adjoint solver
    if k == K
        save2disk(['w_term_' num2str(k) '.dat'],wTn);
    else
        save2disk(['w_term_' num2str(k) '.dat'],wT(k,:));
    end
    adjoint(k,0);
    Pwp(k,:) = load(['w_init_' num2str(k) '.dat']);
    %%%
end

% compute residuals
Rw = v0(2:end,:) - Pvm(1:end-1,:);
Rv = [w0; wT] - Pwp;

y = [reshape(Rv',n*K,1); reshape(Rw',n*(K-1),1)];

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

function plotter
% plot of solution

% read in mesh data from mesh.dat
s = dlmread('mesh.dat', ' ');
T = s(end);
dt = s(end-1);
n = s(end-2);
m = round(T/dt);
s = s(1:end-3);
t = (0:m)*dt;

% read in data
tChk = load('timeChk.dat');
[K tmp] = size(tChk);

u = [];
v = [];
w = [];

tiEnd = cumsum(tChk(:,1));


for k = 1:K
    uk = load(['primal' num2str(k) '.dat']);
    vk = load(['tangent' num2str(k) '.dat']);
    wk = load(['adjoint' num2str(k) '.dat']);
    if k > 1
        u = [u; uk(2:end,:)];
        v = [v; vk(2:end,:)];
        w = [w; wk(end:-1:2,:)];
    else
        u = [u; uk];
        v = [v; vk];
        w = [w; wk(end:-1:1,:)];
    end

end


% w = dlmread('adjoint.dat', ' ');



figure(1);
fntsze = 14;
plot(t,u);

xlabel('t','FontSize',fntsze);
ylabel('u','FontSize',fntsze);
% set(gca,'YTickLabel','')
set(gca,'FontSize',fntsze)

figure(2);
fntsze = 14;
plot(t,v);

xlabel('t','FontSize',fntsze);
ylabel('v','FontSize',fntsze);
% set(gca,'YTickLabel','')
set(gca,'FontSize',fntsze)



% Nv = sqrt(sum(v.^2,2))./sqrt(sum(u.^2,2));
% figure(3);
% plot(1:(m+1),Nv,cumsum(tChk(:,1)),Nv(cumsum(tChk(:,1))),'ro')
% 
% sum(tChk,1)
% t(cumsum(tChk(:,1)))

figure(3);
fntsze = 14;
plot(t,w);

xlabel('t','FontSize',fntsze);
ylabel('w','FontSize',fntsze);
% set(gca,'YTickLabel','')
set(gca,'FontSize',fntsze)

% 
% figure(4);
% plot3(u(:,1),u(:,2),u(:,3))

end