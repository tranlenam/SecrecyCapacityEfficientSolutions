function C = ComputeSecrecyRate(Hb,He,X)
C = max(real(log(det(eye(size(Hb,1))+Hb*X*Hb'))...
    -log(det(eye(size(He,1))+He*X*He'))),0);
end