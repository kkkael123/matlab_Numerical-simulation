%第三章PAA与PWA的对比
clc;
clear;
%gamma = 0.2; 
%alpha = 0.2;
k = 0.5; gamma = 0.41; 
alpha=0.71;
A = 2 - gamma^2;
B = 4 - gamma^2 - alpha .* gamma^2;
cb=0.3;
%cb=0.5;t=0.3;
mu_low1=gamma*(alpha*A+1-alpha)./A;
mu_low2=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alpha.*(1-gamma.^2)))./(2-gamma.^2);
mu_low=max(mu_low1,mu_low2);
mus = (0.6:0.01:1)';

qc=0;
q1 = 1;
minmu1=[];
minmu2=[];
i=1;
while i<=length(mus)
    mu= mus(i);
  t=(0:0.001:1)';  
    %c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));

 
   q2=mu.*q1;
    qn=(1-k-t).*qc+(k+t).*q1;
    qs=(1-k-t).*qc+(k+t).*q2;


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

    %plot(mu,Pi_P_NAA+0*mu,'r'); hold on;
    %plot(mu,Pi_P_NWA+0*mu,'b'); hold off;

    plot(t,Pi_P_PAA+0*mu,'r'); hold on;
    plot(t,Pi_P_PWA+0*mu,'b'); hold off;

    plot(t,Pi_S_PAA+0*mu,'r'); hold on;
    plot(t,Pi_S_PWA+0*mu,'b'); hold off;

    %plot(cb,Pi_S_PAA+0*cb,'r'); hold on;
    %plot(cb,Pi_S_PWA+0*cb,'b'); hold off;
  
    dif = abs(Pi_P_PWA - Pi_P_PAA);
    indmu = find(dif == min(dif));
    muet = t(indmu(1));
    minmu1 = [minmu1; muet];

    dif = abs(Pi_S_PWA - Pi_S_PAA);
    indmu = find(dif == min(dif));
    muet = t(indmu(1));
    minmu2 = [minmu2; muet];
    

    i=i+1;
end
mu_low11=gamma*(alpha*A+1-alpha)./A;
mu_low22=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alpha.*(1-gamma.^2)))./(2-gamma.^2);
mu_low33=max(mu_low11,mu_low22);
plot(mus,minmu1,'-r','MarkerSize',6,'LineWidth',3);hold on;
plot(mus,minmu2,'-b','MarkerSize',6,'LineWidth',3);hold on;
plot(mus,mu_low33,'-k','MarkerSize',6,'LineWidth',3);hold on;
hold off;

xlabel('mu');
ylabel('t');
%xlabel('二手正品的质量{\it{\mu}}');
%ylabel('区块链的实施成本{\it{c_b}}');