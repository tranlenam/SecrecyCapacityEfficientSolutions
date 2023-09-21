function [Xbest,bestobjseqMD,soltime] = Algorithm2_CoMirror_SPC_PAPC(Hb,myPhi,P0,PAPC,X0,maxIter)
X = X0;
Xbest=X;
nA = size(Hb,2);
Rx = 1/2*(P0)^2;
bestobj = -Inf;
objseqMD = zeros(maxIter,1);
bestobjseqMD = zeros(maxIter,1);
%tic
for iMD = 1:maxIter
    if (max(real((diag(X)-PAPC)))<0)
        ek = computegradObj(Hb,myPhi,X);
        objseqMD(iMD) = computeObj(Hb,myPhi,X);
        if(objseqMD(iMD)>bestobj)
            Xbest = X;
            bestobj = objseqMD(iMD);
        end
    else
        [~,m] = max(real(diag(X)-PAPC));
        ek = -diag([zeros(1,max(m-1,0)), 1, zeros(1,max(nA-m,0))]);
    end
    bestobjseqMD(iMD)=max(objseqMD(1:iMD));
    
    tk = (Rx^0.5)/norm(ek,'fro')/sqrt(iMD);
    
    X = (projectSPC(X+tk*ek,P0,'inequality'));
    %{
    if( (iMD>200) && (abs(bestobjseqMD(iMD)-bestobjseqMD(iMD-100))<1e-5))
        bestobjseqMD(iMD:end)=[];
        break
    end
    %}
end
    function y = computeObj(Hb,myPhi,X)
        y = real(log(det(eye(size(Hb,1))+Hb*X*Hb'))-(trace(myPhi*X)));
    end
    function y = computegradObj(Hb,myPhi,X)
        y = (Hb'*((eye(size(Hb,1))+Hb*X*Hb')\Hb)-myPhi);
    end

end

