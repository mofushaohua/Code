%% For convinience, we assume the order of the tensor is always 3;
clear;
addpath('tSVD','proxFunctions','solvers','twist','data-ORL');
addpath('ClusteringMeasure', 'LRR', 'Nuclear_norm_l21_Algorithm', 'unlocbox','gspbox-0.7.0');
load('yaleB.mat'); 
cls_num = length(unique(gt));
%% Note: each column is an sample (same as in LRR)
%data preparation
 X{1} = X3; X{2} = X2; X{3} = X1;
for v=1:3
    [X{v}]=NormalizeData(X{v});
end
% Initialize...
K = length(X); N = size(X{1},2); %sample number
for k=1:K
    Z{k} = zeros(N,N); 
    W{k} = zeros(N,N);
    G{k} = zeros(N,N);
    B{k} = zeros(N,N);
    L{k} = zeros(N,N);
    E{k} = zeros(size(X{k},1),N);
    Y{k} = zeros(size(X{k},1),N); 
end

dim1 = N;dim2 = N;dim3 = K;
myNorm = 'tSVD_1';
sX = [N, N, K];
%set Default
parOP         =    false;
ABSTOL        =    1e-6;
RELTOL        =    1e-4;


Isconverg = 0;epson = 1e-5;
lambda1 = 0.0016; 
iter = 0;
mu1 = 10e-5; max_mu1 = 10e10; pho_mu1 = 1.6;
mu2 = 10e-5; max_mu2 = 10e10; pho_mu2 = 1.6;
rho = 10e-5; max_rho = 10e10; pho_rho = 1.6;
tic;
start = 1;

while(Isconverg == 0)
    start = 0;
    %-------------------update Z^k-------------------------------
    for k=1:K
        Z{k} = inv(rho*eye(N)+mu1*X{k}'*X{k}) * (X{k}'*(Y{k}+mu1*X{k}-mu1*E{k})-W{k}+rho*G{k});
    end
    %-------------------update E^k-------------------------------
    F = [];
    for k=1:K    
        tmp = X{k}-X{k}*Z{k}+Y{k}/mu1;
        F = [F;tmp];
    end  
    [Econcat] = solve_l1l2(F,lambda1/mu1);
    start = 1;
    for k=1:K
        E{k} = Econcat(start:start + size(X{k},1) - 1,:);
        start = start + size(X{k},1);
    end
    %-------------------update G---------------------------------
    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:); 
    [g, objV] = wshrinkObj(z + 1/rho*w,rho,sX,0,3);
    G_tensor = reshape(g, sX);
    %------------------update auxiliary variable---------------
    w = w + rho*(z - g);
    W_tensor = reshape(w, sX);
    for k=1:K
        G{k} = G_tensor(:,:,k);
        Y{k} = Y{k} + mu1*(X{k}-X{k}*Z{k}-E{k});
        W{k} = W_tensor(:,:,k);
    end
    %% coverge condition
    Isconverg = 1;
    for k=1:K
        if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
            history.norm_Z = norm(X{k}-X{k}*Z{k}-E{k},inf);
            fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end
    end
    if (iter>200)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu1 = min(mu1*pho_mu1, max_mu1);
    mu2 = min(mu2*pho_mu2, max_mu2);
    rho = min(rho*pho_rho, max_rho);
    fprintf("\n");
end

%-----------------------------------------------------------
for k=1:K
    W{k} = zeros(N,N);
    G{k} = zeros(N,N);
    E{k} = zeros(size(X{k},1),N); 
    Y{k} = zeros(size(X{k},1),N); 
end

lambda1 = 0.0014;
lambda2 = 0.2;
iter = 0;
mu1 = 10e-5; max_mu1 = 10e10; pho_mu1 = 1.8;
mu2 = 10e-5; max_mu2 = 10e10; pho_mu2 = 1.8;
rho = 10e-5; max_rho = 10e10; pho_rho = 1.8;
start = 1;
%construct K
for k=1:K
    Weight{k} = constructW_PKN(abs(Z{k})+abs(Z{k}')./2, 5);
    A{k}=Weight_vector(Weight{k},N);
end

for k=1:K
    B{k} = zeros(N,size(A{k},1));
    J{k} = zeros(N,size(A{k},1));
end
Isconverg=0;
while(Isconverg==0)
    fprintf('----processing iter %d--------\n', iter+1); 
    % update Z^k
    for k=1:K
        A1=mu1*X{k}'*X{k}+rho*eye(N);
        B1=mu2*A{k}'*A{k};
        C1=(B{k}+mu2*J{k})*A{k}+X{k}'*Y{k}+mu1*X{k}'*(X{k}-E{k})-W{k}+rho*G{k};
        Z{k}=sylvester(A1,full(B1),full(C1));
    end
    % update E^k
    F = [];
    for k=1:K    
        tmp = X{k}-X{k}*Z{k}+Y{k}/mu1;
        F = [F;tmp];
    end
    [Econcat] = solve_l1l2(F,lambda1/mu1);
    start = 1;
    for k=1:K 
        E{k} = Econcat(start:start + size(X{k},1) - 1,:);
        start = start + size(X{k},1);
    end
    % update J^k
    for k=1:K
        J{k}=errormin(B{k},A{k},Z{k},lambda2,mu2,1);
    end
    % update G
    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:);
    [g, objV] = wshrinkObj(z + 1/rho*w,rho,sX,0,3);
    G_tensor = reshape(g, sX);
    w = w + rho*(z - g);
    W_tensor = reshape(w, sX);
    for k=1:K
        Y{k} = Y{k} + mu1*(X{k}-X{k}*Z{k}-E{k});
        B{k} = B{k} + mu2*(J{k}-Z{k}*A{k}');
        G{k} = G_tensor(:,:,k);
        W{k} = W_tensor(:,:,k);
    end
    history.objval(iter+1)   =  objV;
        %% coverge condition
    Isconverg = 1;
    for k=1:K
        if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
            history.norm_Z = norm(X{k}-X{k}*Z{k}-E{k},inf);
            fprintf('    norm_Z %7.10f  \n  ', history.norm_Z);
            Isconverg = 0;
        end
    end
    if (iter>200)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu1 = min(mu1*pho_mu1, max_mu1);
    mu2 = min(mu2*pho_mu2, max_mu2);
    rho = min(rho*pho_rho, max_rho);
end
S1 = 0;
for k=1:K
    S1 = S1 + abs(Z{k})+abs(Z{k}');
end
C1 = SpectralClustering(S1,cls_num);
%------------nmi-------------------------
[A1 nmi1 avgent1] = compute_nmi(gt,C1);
%------------acc-------------------------
ACC1 = Accuracy(C1,double(gt));
[f1,p1,r1] = compute_f(gt,C1);
%------------ARI-------------------------
[AR1,RI1,MI1,HI1]=RandIndex(gt,C1);
fprintf('----using Z clustering: ACC %f, NMI %f, AR %f---------\n', ACC1, nmi1, AR1);

