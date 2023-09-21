function [Xprev,objseq] = Algorithm1ADCA_SPC_PAPC_IPC(H,G,P0,PAPC,IPC, S0,maxIter,q)
t = (1+sqrt(5))/2;
Xprev = S0;
V = Xprev;% extrapolated point

[nB,nA] = size(H);
nE = size(G,1);

objseq=zeros(maxIter+1,1);
objseq(1)= ComputeSecrecyRate(H,G,Xprev);
solveroptions = sdpsettings('verbose',0,'solver','SDPT3');
X = sdpvar(nA,nA,'hermitian','complex');
F = [X>=0,real(trace(X))<=P0,real(diag(X))<=PAPC,...
    real(trace(IPC{1,1}*X))<=IPC{1,2}, real(trace(IPC{2,1}*X))<=IPC{2,2}];
cbj = objseq(1); % current best objective
for iIter=1:maxIter
    Phi = G'*((eye(nE)+G*V*G')\G);
    obj = logdet(eye(nB)+H*X*H')-real(trace(Phi*X));
    diagnotics = optimize(F,-obj,solveroptions);
    if (diagnotics.problem==0)
        Xnew = value(X);
        objseq(iIter+1)=ComputeSecrecyRate(H,G,Xnew);
        if(cbj<objseq(iIter+1))
            cbj = objseq(iIter+1);
            cbsignal = Xnew; % current best input signalling
        end
        tnew = (1+sqrt(1+4*t^2))/2;

        Z = Xnew+(t-1)/tnew*(Xnew-Xprev);

        FZ = ComputeSecrecyRate(H,G,Z);
        if (FZ>=min(objseq(max(0,iIter-q)+1:iIter+1)))
            V=Z;
        else
            V=Xnew;
        end
        Xprev = Xnew;
        t = tnew;
    else
        break
    end

end

end

