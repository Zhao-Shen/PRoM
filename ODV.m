function v = ODV(r1,r2,m,q1,q2)
m1 = q1*m;
m2 = q2*m;

Grid = m1/2/1000;
Omega = lcm(q1,q2)*m;

v_candidate = -Omega/2:Grid:Omega/2-Grid;

cost = 1- cos(mod(v_candidate, m1)./m1*2*pi-mod(r1,m1)./m1*2*pi) + ...
       1- cos(mod(v_candidate, m2)./m2*2*pi-mod(r2,m2)./m2*2*pi);

[~, indx]= min(cost);
v = v_candidate(indx);

end

