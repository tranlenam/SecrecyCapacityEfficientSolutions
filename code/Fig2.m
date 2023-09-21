clear
clc
%% channel generation
nA  = 2;
nB = 4;
nE = 3;



Hb=-[0.397388134804813 + 0.564069068925790i	-0.0939337320157974 + 0.253217235750467i
-0.0216281544185800 + 0.805107628931355i	-0.673432238998531 + 0.260482934998363i
-1.19032291480353 - 0.393916696585856i	-0.972781640103893 - 0.446827692307842i
0.201667162510618 - 0.689726745440122i	-0.945047869851278 - 0.730641339097362i];

He = [-0.201455080897275 + 0.312725917599147i	-0.617804134199973 - 1.04797479362569i
-0.0558694759682173 - 0.300019080432181i	-0.385799249500664 - 0.281697157494074i
0.693547640268284 + 0.0558760342293484i	-0.506416527519376 - 0.144333099022521i];
[nB,nA] = size(Hb);
nE = size(He,1);

%% transmit power parameters
SNRdB = [5 15];
figure
for iSNRbdB=1:length(SNRdB)
    P0 = 10^(SNRdB(iSNRbdB)/10);
    PAPC = ones(nA,1)*P0/nA*1.2;

    %% run iterative algorithms of comparison
    S = eye(nA)*P0/(nA*2);


    mylambda0 = ones(nA+1,1);
    maxIter = 2000; % maximum number of iterations

    Phi = He'*((eye(size(He,1))+He*S*He')\He);
    c = real(-log(det(eye(size(He,1))+He*S*He'))+trace(Phi*S));

    %% Find optimal solution using convex solver
    solveroptions = sdpsettings('verbose',0,'solver','SDPT3');
    X = sdpvar(nA,nA,'hermitian','complex');
    F = [X>=0,(trace(X)<=P0),diag(X)<=PAPC];
    obj = logdet(eye(nB)+Hb*X*Hb')-real(trace(Phi*X));
    solution=optimize(F,-obj,solveroptions);
    if solution.problem==0
        S = value(X);
    end
    SecrecyCapacityLB = real(double(obj)) + c;

    %% Algorithm 2, CoMirror Method
    X0 = zeros(nA);
    [solCoMirror,objseqCoMirror]=Algorithm2_CoMirror_SPC_PAPC(Hb,Phi,P0,PAPC,X0,maxIter);

    %% Subgradient method
    [Sbest,objseqprimal,lambdabest,objseqdual] = SubgradMaxdetSPCPAPC(Hb,Phi,P0,PAPC,mylambda0,maxIter);


    %% plot results
    SecrecyCapacityLB_CoMirror = max(objseqCoMirror+c,0);
    semilogx(SecrecyCapacityLB*ones(length(SecrecyCapacityLB_CoMirror),1),'k')
    hold on
    semilogx(SecrecyCapacityLB_CoMirror,'b');
    % plot(objseqdual,'k')
    semilogx(max(objseqprimal+c,0),'r');
end
xlabel('Interation index')
legend('Optimal Objective (Convex Solver-SDPT3)','Algorithm 2','Subgradient method')
legend('Location','southeast')
saveas(gcf,'../results/Fig2.png')
