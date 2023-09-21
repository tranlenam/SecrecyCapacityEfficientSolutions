function [Sbest,bestobj_seq_primal,lambdabest,obj_seq_dual] = SubgradMaxdetSPCPAPC(Hb,myPhi,P0,PAPC,mylambda0,maxIter)

% this function is to compute the optimal solution to the following problem
% max log|I+Hb*S*Hb'|-trace(myPhi*S)
% subject to: S>=0; trace(S) <= P0 ; [S]_{n,n} <= PAPC(n)
% using subgradient method

[nB,nA] = size(Hb);
mylambda = mylambda0;
stepsize = 0.001;
bestobj_primal = -Inf;
bestobj_dual = Inf;
Sbest = zeros(nA);
obj_seq_primal = zeros(maxIter,1);
bestobj_seq_primal= zeros(maxIter,1);
obj_seq_dual=zeros(maxIter,1);
lambdabest = zeros(nA+1,1);

for iIter = 1:maxIter
   
    
    R = sqrtm(myPhi+diag(mylambda(1)+mylambda(2:end)));
    
    [U, Sig, ~]=svd(R\(Hb'),'econ');
    sigma_bar = diag(Sig).^2; 
    
    A = (R')\U; %inv(R')*U 
    
    S = A*diag(max(1-1./sigma_bar,0))*A';
    
    
    subgradLambda = real([P0-trace(S);PAPC - diag(S)]);
    
    if (min(subgradLambda)>=0)
        obj_primal=real(log(det(eye(nB)+Hb*S*Hb'))-real(trace(myPhi*S)));
        obj_seq_primal(iIter) = obj_primal;
        if(obj_primal>bestobj_primal)
            Sbest = S;
            bestobj_primal = obj_primal;
        end
    end
    bestobj_seq_primal(iIter) = max(obj_seq_primal(1:iIter));
    % store the true objective
    
    obj_dual = obj_primal + mylambda(2:end)'*PAPC;
    
    obj_seq_dual(iIter) = min(obj_dual,bestobj_dual);
    if(obj_dual < bestobj_dual)
        bestobj_dual = min(obj_dual,bestobj_dual);
        lambdabest = mylambda;
    end
    mylambda = max(mylambda - stepsize*subgradLambda,0);
    
   
end






