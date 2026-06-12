%第三章PAA与PWA的对比
clc;
clear;
%gamma = 0.2; 
%alpha = 0.2;
k = 0.5; 
t=0.1;
gamma = 0.1; 
alpha = 0.7;
A = 2 - gamma^2;
B = 4 - gamma^2 - alpha .* gamma^2;
q1 = 1; 
mu_low=gamma*(alpha*A+1-alpha)./A;
mumin=gamma./A;
qc=0;



mus = (mumin:0.01:1)';
mincb1=[];
mincb2=[];
mincb3=[];
i=1;
while i<=length(mus)
    mu= mus(i);
    q2=mu.*q1;
    qn=(1-k-t).*qc+(k+t).*q1;
    qs=(1-k-t).*qc+(k+t).*q2;

    %c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));
    %c_bmax=((2-gamma.^2).*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) ))./((1-gamma).*(2-  gamma.^2));

cbmax1=(2.*mu+gamma.*(-1+alpha.*(-1+gamma.^2)-gamma.*mu)).*(k + t)./((1-gamma).*(2-gamma.^2).*t);

    cb = (0:0.01:0.5)'; 


    % P-AA 平台承担AI技术-代理模式
PN_PAA=(gamma^2*k - 2*t - 2*k + gamma^2*t + gamma*k*mu + gamma*mu*t)./(gamma^2 - 4);
PS_PAA=(gamma*k - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t + gamma^2*k*mu)./(gamma^2 - 4);

DN_PAA=qn-PN_PAA-gamma.*(qs-PS_PAA);
DS_PAA=qs-PS_PAA-gamma.*(qn-PN_PAA);
Pi_P_PAA=DN_PAA .* (alpha .* PN_PAA-t.*cb) + DS_PAA .* (alpha .* PS_PAA-t.*cb);
Pi_S_PAA= (1 - alpha) .* (DS_PAA .* PS_PAA);
Pi_N_PAA= (1 - alpha) .* (DN_PAA .* PN_PAA);


    % P-WA 平台承担AI技术-批发模式
    % 均衡解
w_PWA=(gamma*k + 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t + alpha*gamma*k + alpha*gamma*t - 2*cb*gamma*t - alpha*gamma^3*k - alpha*gamma^3*t - cb*gamma^2*t + cb*gamma^3*t + gamma^2*k*mu)./(2*gamma^2 - 4);
PN_PWA=(gamma^2*k - 2*t - gamma*w_PWA - 2*k + gamma^2*t - cb*gamma*t + gamma*k*mu + gamma*mu*t + cb*gamma^2*t)./(alpha*gamma^2 + gamma^2 - 4); 
PS_PWA=(gamma*k - 2*w_PWA - 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t - alpha*gamma*k - alpha*gamma*t + 2*cb*gamma*t + gamma^2*k*mu + alpha*gamma^2*k*mu + alpha*gamma^2*mu*t)./(alpha*gamma^2 + gamma^2 - 4);

DN_PWA=qn-PN_PWA-gamma.*(qs-PS_PWA);
DS_PWA=qs-PS_PWA-gamma.*(qn-PN_PWA);
Pi_P_PWA= DN_PWA .* (alpha .* PN_PWA-t.*cb ) + DS_PWA.* (PS_PWA-w_PWA-t.*cb);
Pi_S_PWA= w_PWA .* DS_PWA;
Pi_N_PWA= (1 - alpha) .* (DN_PWA .* PN_PWA);

    plot(cb,Pi_P_PAA+0*mu,'r'); hold on;
    plot(cb,Pi_P_PWA+0*mu,'b'); hold off;

    plot(cb,Pi_S_PAA+0*cb,'r'); hold on;
    plot(cb,Pi_S_PWA+0*cb,'b'); hold off;

    dif = abs(Pi_P_PWA - Pi_P_PAA);
    indcb = find(dif == min(dif));
    cbet = cb(indcb(1));
    mincb1 = [mincb1; cbet];

    dif = abs(Pi_S_PWA - Pi_S_PAA);
    indcb = find(dif == min(dif));
    cbet = cb(indcb(1));
    mincb2 = [mincb2; cbet];
    
    dif = abs(Pi_N_PWA - Pi_N_PAA);
    indcb = find(dif == min(dif));
    cbet = cb(indcb(1));
    mincb3 = [mincb3; cbet];

    i=i+1;
end

c_bmax1=((2-gamma.^2).*mus+gamma.*(-1 + alpha .*(-1 +gamma.^2) ))./((1-gamma).*(2-  gamma.^2));

plot(mus,mincb1,'-r','MarkerSize',6,'LineWidth',3);hold on;
plot(mus,mincb2,'--b','MarkerSize',6,'LineWidth',3);hold on;
%plot(mus,mincb3,'k','MarkerSize',6,'LineWidth',3);hold on;
%plot(mus,c_bmax1,'-k','MarkerSize',6,'LineWidth',3);hold on;

% 平滑处理

%Pi_p_c_b= smooth(mus, mincb1, 0.2, 'lowess'); % 使用局部加权回归 (lowess)
%Pi_m_c_b= smooth(mus, mincb2, 0.2, 'lowess');
% 绘图

%plot(ms, ms2min, 'b-', 'DisplayName', '原始数据'); hold on;
%plot(mus, Pi_p_c_b, 'r', 'LineWidth', 2, 'DisplayName', '平台'); hold on;
%plot(mus, Pi_m_c_b, 'b', 'LineWidth', 2, 'DisplayName', '制造商'); hold off;

hold off;

%xlabel('二手正品的质量{\it{\mu}}');
%ylabel('区块链的实施成本{\it{c_b}}');