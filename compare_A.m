%第三章PAA与PWA的对比
clc;
clear;

k = 0.3; gamma = 0.2; 
alpha=0.5;
A = 2 - gamma^2;
B = 4 - gamma^2 - alpha .* gamma^2;
cb=0.3;
%cb=0.5;t=0.3;
mu_low1=gamma*(alpha*A+1-alpha)./A;
mu_low2=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alpha.*(1-gamma.^2)))./(2-gamma.^2);
mu_low=max(mu_low1,mu_low2);
mus = (0.4:0.01:1)';

qc=0;
q1 = 1;
minmu1=[];
minmu2=[];
i=1;
while i<=length(mus)
    mu= mus(i);
  t=(0.1:0.001:1-k)';  
    %c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));

 
   q2=mu.*q1;
qn = k * q1;
qs=  k*q2; 
    qn_A=(1-k-t).*qc+(k+t).*q1;
    qs_A=(1-k-t).*qc+(k+t).*q2;

    % N-AA 无AI技术-代理模式
    % 均衡解
    PS_NAA = (A .* qs - gamma .* qn) / (2 + A);
    PN_NAA = (A .* qn - gamma .* qs) / (2 + A);
    % 需求函数
    DS_NAA = (A .* qs - gamma .* qn) / (2 + A);
    DN_NAA = (A .* qn - gamma .* qs) / (2 + A);
    %利润函数
    Pi_P_NAA= DN_NAA .* (alpha .* PN_NAA ) + DS_NAA .* (alpha .* PS_NAA);
    Pi_S_NAA= (1 - alpha) .* (DS_NAA .* PS_NAA);
    Pi_N_NAA= (1 - alpha) .* (DN_NAA .* PN_NAA);

PN_PAA=(gamma^2*k - 2*t - 2*k + gamma^2*t + gamma*k*mu + gamma*mu*t)./(gamma^2 - 4);
PS_PAA=(gamma*k - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t + gamma^2*k*mu)./(gamma^2 - 4);

DN_PAA=qn_A-PN_PAA-gamma.*(qs_A-PS_PAA);
DS_PAA=qs_A-PS_PAA-gamma.*(qn_A-PN_PAA);
Pi_P_PAA=DN_PAA .* (alpha .* PN_PAA-t.*cb) + DS_PAA .* (alpha .* PS_PAA-t.*cb);
Pi_S_PAA= (1 - alpha) .* (DS_PAA .* PS_PAA);
Pi_N_PAA= (1 - alpha) .* (DN_PAA .* PN_PAA);

    %plot(mu,Pi_P_NAA+0*mu,'r'); hold on;
    %plot(mu,Pi_P_NWA+0*mu,'b'); hold off;

    plot(t,Pi_P_NAA+0.*t,'r'); hold on;
    plot(t,Pi_P_PAA,'b'); hold off;

    plot(t,Pi_S_NAA+0.*t,'r'); hold on;
    plot(t,Pi_S_PAA,'b'); hold off;

    %plot(cb,Pi_S_PAA+0*cb,'r'); hold on;
    %plot(cb,Pi_S_PWA+0*cb,'b'); hold off;
  
abs4=Pi_P_PAA - Pi_P_NAA;

    plot(t,abs4,'b'); hold off;


    dif = abs(abs4-0);
    indmu = find(dif == min(dif));
    muet = t(indmu(1));
    minmu1 = [minmu1; muet];




   % dif = abs(Pi_S_PAA - Pi_S_NAA);
    %indmu = find(dif == min(dif));
   % muet = t(indmu(1));
   % minmu2 = [minmu2; muet];
    

    i=i+1;
end
mu_low11=gamma*(alpha*A+1-alpha)./A;
mu_low22=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alpha.*(1-gamma.^2)))./(2-gamma.^2);
mu_low33=max(mu_low11,mu_low22);
plot(mus,minmu1,'-r','MarkerSize',6,'LineWidth',3);hold on;
%plot(mus,minmu2,'-b','MarkerSize',6,'LineWidth',3);hold on;


tmax11=k.*(2.*mus+gamma.*(-1+alpha.*(-1+gamma.^2)-gamma.*mus))./(cb.*(1-gamma).*(2-gamma.^2)-2.*mus+gamma.*(1+alpha.*(1-gamma.^2)+gamma.*mus));


%plot(mus,tmax11,'-k','MarkerSize',6,'LineWidth',3);hold on;
hold off;

xlabel('mu');
ylabel('t');
%xlabel('二手正品的质量{\it{\mu}}');
%ylabel('区块链的实施成本{\it{c_b}}');