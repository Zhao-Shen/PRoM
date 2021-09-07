function [output] = PRoM3( r1,r2,r3,m1,m2,m3,cov_inv)
% m1 m2 m3 are integers
%% Initialize several vairables
persistent  h Omega k1List k2List k3List

R = m2/m1;
if isempty(Omega)||isempty(h)|| isempty(k1List) || isempty(k2List) || isempty(k3List)
    
    Omega = lcm(lcm(m1,m2),m3);
    h = Omega./[m1, m2, m3];
    
    % number of gcd intervals searching range
    i=1;
    for k1 = -1:h(1)
        for k2 = -1:h(2)
            for k3 = -1:h(3)
                if max([(k1-1)*m1,(k2-1)*m2,(k3-1)*m3]) < min([(k1+1)*m1,(k2+1)*m2,(k3+1)*m3])
                    k1List(i)=k1;
                    k2List(i)=k2;
                    k3List(i)=k3;
                    i=i+1;
                end
            end
        end
    end
    
    clear k1 k2 k3;
    indx1 = [];
    for i = 1: numel(k1List)
        indx1 = [indx1, find((k1List == k1List(i)+h(1)).*(k2List == k2List(i)+h(2)).*(k3List == k3List(i)+h(3)))];        
    end
    k1List(indx1) = [];
    k2List(indx1) = [];
    k3List(indx1) = [];
end
%% Search
v_Hat_Candidate = zeros(1,numel(k1List));
cost = v_Hat_Candidate;

weights = 1/(2*(R^2-R+1))*[R^2 1 (-1+R)^2]';


v_Hat_Candidate = mod(weights'*([k1List;k2List;k3List].*[m1;m2;m3]+[r1;r2;r3]),Omega);

for i = 1: numel(v_Hat_Candidate)
    %     x = [Circular_Distance(r1,v_Hat_Candidate(i),m1);...
    %         Circular_Distance(r2,v_Hat_Candidate(i),m2);
    %         Circular_Distance(r3,v_Hat_Candidate(i),m3)];
    
    x = [v_Hat_Candidate(i)-r1 - round((v_Hat_Candidate(i)-r1)/m1)*m1;...
        v_Hat_Candidate(i)-r2 - round((v_Hat_Candidate(i)-r2)/m2)*m2;...
        v_Hat_Candidate(i)-r3 - round((v_Hat_Candidate(i)-r3)/m3)*m3];
    
    cost(i) = round(x'*cov_inv*x,8);
    %     To avoid the misleading due to Matlab precision onto the cost
end
[~,indx2] = min(cost(:));
k1Hat = k1List(indx2);
k2Hat = k2List(indx2);
k3Hat = k3List(indx2);
nHat = v_Hat_Candidate(indx2);
nHat = (nHat < Omega/2).*nHat+(nHat > Omega/2).*(nHat-Omega);
output = struct('k1Hat',k1Hat,'k2Hat',k2Hat,'k3Hat',k3Hat,'nHat',nHat,'mincost',cost(indx2),'cost',cost);


end

