function [X,SecrecyCapacity,SecrecyRateSeq] = Algorithm3PBRA_SPC_PAPC(Hb,He,P0,PAPC,maxIter)
[nB,nA] = size(Hb);
nE = size(He,1);
K = eye(nB+nE);
H  = [Hb;He];
SecrecyRateSeq = zeros(maxIter,1);
solveroptions = sdpsettings('verbose',0,'solver','SDPT3');
S = sdpvar(nA,nA,'hermitian','complex');
Y = sdpvar(nA,nA,'hermitian','complex');
F = [S>=0,Y>=0,real(trace(S))<=P0,real(diag(S))<=PAPC];
for iIter = 1:maxIter
    % Optimize wrt to X, (K fixed)
    H1 = K^(-0.5)*H;
    Delta = H1'*H1-He'*He;
    Delta_bar = sqrtm(Delta);
    
    obj = logdet(Y);
    
    diagnostic = optimize([F,[eye(nA)+Delta_bar*S*Delta_bar'-Y  Delta_bar*S*He';
        He*S'*Delta_bar'   eye(nE)+He*S*He']>=0],-obj,solveroptions);
    if(diagnostic.problem==0)
        X = value(S);
    else
      X = NaN;
      SecrecyRateSeq  = NaN;
      break

    end
    SecrecyRateSeq(iIter) = ComputeSecrecyRateMinMaxObj(Hb,He,X,K);
    % Optimize wrt K (X fixed)
   
    mypsi = inv(K+H*X*H');
    myphi12 = mypsi(1:nB,nB+1:nB+nE);
    [U,D] = eig(myphi12*myphi12');
    d = real(diag(D));
    K_bar1 = -2*U*diag(1./(1+sqrt(1+4*d)))*U'*myphi12;
    K = [eye(nB) K_bar1;K_bar1' eye(nE)];
    % check convergence
    if (iIter>10)
        if(abs(SecrecyRateSeq(iIter)-SecrecyRateSeq(iIter-10))<=1e-5)
            break
        end
    end
    
end
SecrecyRateSeq(iIter+1:end)=[];
SecrecyCapacity = SecrecyRateSeq(end);
end



