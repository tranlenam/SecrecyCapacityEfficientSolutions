function [Y] = projectSPC(X,P0,flag)
%proj_SPC this function project X onto the set trace(S)<= P0 
[U,D] = eig(X);
%[U,D] = eig(X,X,'chol');

d = real(diag(D));
if strcmp(flag,'equality') % trace(S)=P0
    d1 = projsplx(d/P0);
    Y = U*diag(d1*P0)*U';
elseif strcmp(flag,'inequality') % trace(S) <= P0
    if(sum(max(d,0))<P0)
        d1 = max(d,0);
    else
        d1 = projsplx(d/P0)*P0;
    end
    Y = U*diag(d1)*U';
else
    error('Invalid flag')
    
end

end

