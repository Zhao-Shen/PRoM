close all;
clear all;
clc;

snr = 5;
p = 7;
q = 5;
R = p/q;                                % 1<= R <= 2

if gcd(p,q) ~= 1
    error('p, q are not coprime !');
end

venc = [q*(p-q);    p*(p-q);     p*q]; % to make it integer


COV = venc(1)^2/(2*pi^2*snr^4)*...
    [   1+2*snr^2           snr^2*R             snr^2*R/(R-1);
    snr^2*R             (1+2*snr^2)*R^2     -snr^2*R^2/(R-1);
    snr^2*R/(R-1)     	-snr^2*R^2/(R-1)   	(1+2*snr^2)*R^2/(R-1)^2];

disp(COV)


h = lcm(lcm(venc(1),venc(2)),venc(3))./venc;

Grid = 1/5;
x = 0:Grid:1*2*venc(1)-Grid;
y = 0:Grid:1*2*venc(2)-Grid;
z = 0:Grid:1*2*venc(3)-Grid;
[r1,r2,r3] = meshgrid(x,y,z);

k1 = zeros(size(r1));
k2 = zeros(size(r1));
k3 = zeros(size(r1));

indx = (round(r1.*r2.*r3,6)==0);

r1 = r1(indx);
r2 = r2(indx);
r3 = r3(indx);

clear PRoM3

COV_Inv = pinv(COV);

tic
parfor i =1:numel(r1)
    output1 = PRoM3( r1(i),r2(i),r3(i),2*venc(1),2*venc(2),2*venc(3),COV_Inv);
    k11(i) = output1.k1Hat;
    k12(i) = output1.k2Hat;
    k13(i) = output1.k3Hat;
    nHat1(i) = output1.nHat;
    
end
toc



%% k list
m1 = 2*venc(1);
m2 = 2*venc(2);
m3 = 2*venc(3);
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


figure(1);
colors = linspecer(21);
k = 1;
for i = 1: numel(k1List)
    indx1 = find((k11==k1List(i)).*(k12==k2List(i)).*(k13==k3List(i)));
    if numel(indx1) > 10
        disp([k1List(i),k2List(i),k3List(i),numel(indx1)])
        scatter3(r1(indx1),r2(indx1),r3(indx1),'o','MarkerFaceColor',colors(k,:), 'MarkerEdgeColor','none'); hold on
        k = k+1;
    end
    
end


title('$k^\star (\tilde{r})$ ','Interpreter','Latex','FontSize',36)
xlabel('$\tilde{r}_1$','Interpreter','latex','FontSize',36)
ylabel('$\tilde{r}_2$','Interpreter','latex','FontSize',36)
zlabel('$\tilde{r}_3$','Interpreter','latex','FontSize',36)
axis image

xlim([0,2*venc(1)]);
ylim([0,2*venc(2)]);
zlim([0,2*venc(3)]);


legend('$k^\star = [-1,-1,-1]^T$','$k^\star = [-1,-1,0]^T$','$k^\star = [-1,0,-1]^T$', '$k^\star = [-1,0,0]^T$','$k^\star = [0,-1,0]^T$',...
    '$k^\star = [0,0,-1]^T$','$k^\star = [0,0,0]^T$','$k^\star = [0,1,0]^T$','$k^\star = [0,1,1]^T$','$k^\star = [1,1,0]^T$',...
    '$k^\star = [1,2,1]^T$','$k^\star = [2,2,0]^T$','$k^\star = [2,2,1]^T$','$k^\star = [2,3,1]^T$','$k^\star = [3,2,1]^T$',...
    '$k^\star = [3,3,1]^T$','$k^\star = [3,3,2]^T$','$k^\star = [4,3,1]^T$','$k^\star = [4,4,2]^T$','$k^\star = [5,4,1]^T$',...
    '$k^\star = [5,4,2]^T$','Interpreter','Latex','FontSize',36,'Location','northeastoutside')
legend('boxoff')

view([1 1 1])
set(gca,'FontSize',36);



