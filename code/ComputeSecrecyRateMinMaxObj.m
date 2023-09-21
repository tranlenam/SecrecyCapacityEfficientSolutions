function C = ComputeSecrecyRateMinMaxObj(Hr,He,R,K)
H = [Hr;He];
C = real(log(det(K+H*R*H'))-log(det(eye(size(He,1))+He*R*He'))-log(det(K)));
end