bin_path = './bin';
addpath(genpath(bin_path));

% -- loading input network --
GeneNodeFileDir = './network/index_genes.txt';
NetworkFileDir = './network/edge_list.txt';

fprintf('Loading network files ...\n');
[net_map,Adj_raw] = PreprocessNetwork(GeneNodeFileDir,NetworkFileDir);
clear GeneNodeFileDir NetworkFileDir

fprintf('Network loaded ...\n');

% -- loading input mutation data & patient' cancer types --
load('./MutationData/mutation_matrices.mat');
[X_Net,X_NoNet,Symbol_Net,Symbol_NoNet] = ...
    A00_00_InputToNetMat(net_map,m_cdata.X,m_cdata.geneID,...
    [bin_path '/GeneIDPreprocess'],1);
ind_nan = strcmp(Symbol_NoNet,'Unknown');
X_NoNet(:,ind_nan) = [];
Symbol_NoNet(ind_nan) = [];
Symbol_input = [Symbol_Net Symbol_NoNet];
X = [X_Net X_NoNet];

n_gene = size(X,2);
GeneNumNet = size(X_Net,2);
n_disease = m_cdata.num_class;

CancerNameList = m_cdata.className;
U0 = m_cdata.F;
clear m_cdata X_Net X_NoNet Symbol_Net Symbol_NoNet net_map

fprintf('Mutation data loaded ...\n');

% -- gene aligning --
Adj_mat = speye(n_gene);
Adj_mat(1:GeneNumNet,1:GeneNumNet) = Adj_raw;
degree_vec = sum(Adj_mat, 2);
D_mat_half_inv = spdiags(degree_vec.^(-0.5),0,n_gene,n_gene);
A_mat_norm = D_mat_half_inv*Adj_mat*D_mat_half_inv;
D_mat_norm = speye(n_gene);
L_mat_norm = D_mat_norm - A_mat_norm;
clear Adj_raw Adj_mat degree_vec D_mat_half_inv

fprintf('Alignment done ...\n');

rmpath(genpath(bin_path));
clear bin_path

% -- input Disease Ontology based similarites --
load('./DOSim/CancerSimilarity.mat');
S0 = CancerSimilarity;
clear CancerSimilarity

fprintf('Disease Ontology loaded ...\n');
                        
% -- initilaize parameters --
eps_t = 1e-10;
lambda_S = 1; 
lambda_L = 1;
lambda_V = 0.01;

% -- initial values of two matrices --
S_prev = S0;
V_prev = max(X'*(pinv((U0*S_prev))'), eps_t);

Obj = zeros(4,1);
str_display = '[%d iter] cost_func: %2.3e = (%2.3e + %2.3e + %2.3e + %2.3e)\n';

merr_S = 2; merr_V = 2; merr_th = 1e-20;
iter = 0; Maxiters = 500;
while ((merr_S > merr_th) || (merr_V > merr_th)) && iter < Maxiters
    iter = iter+1;
    
    %--- update S ---
    Numer_S = (U0'*X*V_prev);
    Denum_S = (U0'*U0)*S_prev*(V_prev'*V_prev);
    Numer_S_bar = 0.5*Numer_S + 0.5*Numer_S';
    Denum_S_bar = 0.5*Denum_S+0.5*Denum_S';
    clear Numer_S Denum_S
    
    % multiplicative factor
    S_hat = S_prev.*(Numer_S_bar + S_prev.*Denum_S_bar + S_prev.*(lambda_S*S_prev) + lambda_S*S0) ...
        ./(S_prev.*Numer_S_bar + Denum_S_bar + (lambda_S*S_prev) + S_prev.*lambda_S*S0 + eps_t);
    S = max(S_hat,eps_t);
    clear Numer_S_bar Denum_S_bar S_hat
                           
    %--- update V ---               
    Numer_V = (X'*U0*S);
    Denum_V = (V_prev*S'*(U0'*U0)*S);
              
    % multiplicative factor                            
    V_hat = V_prev.*(Numer_V + lambda_L*(A_mat_norm*V_prev))./(Denum_V + lambda_L*D_mat_norm*V_prev + lambda_V + eps_t);
    V = max(V_hat,eps_t);
    clear Numer_V Denum_V V_hat
    
    merr_S = sum(sum((S - S_prev).^2))/sum(sum(S_prev));
    merr_V = sum(sum((V - V_prev).^2))/sum(sum(V_prev));
    %- diplay the objective function value
    X_hat = (U0*S)*V';
    
    Obj(1) = sum(sum((X - X_hat).^2)); % ||X-(U*S)*V'||_W^2
    clear X_hat
    Obj(2) = lambda_S*(sum(sum( (S-S0).^2 )));
    Obj(3) = lambda_V*(sum(sum(V)));
    Obj(4) = lambda_L*(trace(V'*L_mat_norm*V));
    if ~mod(iter,10)
        fprintf(str_display, iter, sum(Obj), Obj(:));
    end
    
    S_prev = S;
    V_prev = V;
end % end of while

[~, ind_gene] = sort(max(V,[],2),'descend');
Candidates_list = Symbol_input(ind_gene(1:200));

BaseResDir = './Output';
mkdir(BaseResDir);
save([BaseResDir '/Results.mat'],'S','V','Candidates_list');