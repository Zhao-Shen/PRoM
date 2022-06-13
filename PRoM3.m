function [K,v_k,L] = PRoM3( v123,venc,cov_inv)
% v123 are measured noisy (possibly wrapped) velocity
% venc is venc
% cov_inv is the inverse of covariance matrix
% K is the wrapping integers candidates
% v_k is the velocity candidates
% L is the cost

persistent Omega h 
if isempty(Omega) || isempty(h) 
    Omega = 2*double(lcm(lcm(sym(venc(1:2))),sym(venc(3))));
    h = round(Omega./venc/2);    
end

% Building K
v = [v123(1)+venc(1):2*venc(1):Omega,v123(2)+venc(2):2*venc(2):Omega,v123(3)+venc(3):2*venc(3):Omega];
K = ceil((v-v123-1e-4)/2./venc-0.5);                                        % -1e-4 is to counter rounding error in ceil function for jumping points

%% Search
w = cov_inv*ones(3,1)/(ones(1,3)*cov_inv*ones(3,1));
v_k = mod(w'*(2*K.*venc+v123),Omega);

x = v_k - v123 - round((v_k - v123)./2./venc)*2.*venc;
L = sum(x.*(cov_inv*x));                                                    % vectorize to avoid for loop 

[L,indx] = sort(L);
v_k = (v_k < Omega/2).*v_k+(v_k > Omega/2).*(v_k-Omega);
K = K(:,indx);
v_k = v_k(indx);


end

