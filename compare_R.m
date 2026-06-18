%第三章PAA与PWA的对比
clc;
clear;

k = 0.1; gamma = 0.7; 
alpha=0.1;
A = 2 - gamma^2;
B = 4 - gamma^2 - alpha .* gamma^2;
cb=0.1;
%cb=0.5;t=0.3;
mu_low1=gamma*(alpha*A+1-alpha)./A;
mu_low2=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alpha.*(1-gamma.^2)))./(2-gamma.^2);
mu_low=max(mu_low1,mu_low2);
mus = (0.6:0.01:0.9)';

qc=0;
q1 = 1;
minmu1=[];
minmu2=[];
i=1;
while i<=length(mus)
    mu= mus(i);
   tmax=k.*(2.*mu+gamma.*(-1+alpha.*(-1+gamma.^2)-gamma.*mu))./(cb.*(1-gamma).*(2-gamma.^2)-2.*mu+gamma.*(1+alpha.*(1-gamma.^2)+gamma.*mu));
 tmax2=-(8*k - 9*gamma^2*k + 2*gamma^4*k - 2*gamma*k*mu - alpha*gamma^2*k + alpha*gamma^4*k + gamma^3*k*mu)/(2*cb*gamma - 2*gamma*mu - alpha*gamma^2 + alpha*gamma^4 - 2*cb*gamma^2 - cb*gamma^3 + cb*gamma^4 + gamma^3*mu - 9*gamma^2 + 2*gamma^4 + 8);

  t=(0.1:0.01:1-k)';  
    %c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));

 
   q2=mu.*q1;
qn = k * q1;
qs=  k*q2; 
    qn_A=(1-k-t).*qc+(k+t).*q1;
    qs_A=(1-k-t).*qc+(k+t).*q2;

    % N-AA 无AI技术-代理模式
    w_NWA = (A .* qs + gamma .* (alpha - 1 - alpha .* A) .* qn) / (2 .* A);
    PS_NWA = (((alpha + 1) .* A.^2 + (1 - 2 .* alpha) .* A) * qs + gamma .* (alpha - 1 - A) * qn) ./ (A .* B);
    PN_NWA = (((alpha + 2) .* A.^2 + (1 - 3 .* alpha) .* A + 2.*(alpha-1)) .* qn + gamma .* (-A) .* qs)./ (2 .* A .* B);
    % 需求函数
    DS_NWA =((A.^2).*qs+(-alpha.*(A.^2)+(alpha-1).*A).*gamma.*qn)./(2.*A.*B);
    DN_NWA =((-A).*gamma.*qs+((alpha+2).*A.^2+(-3.*alpha+1).*A+2.*alpha-2).*qn)./(2.*A.*B);
    %利润函数
    Pi_P_NWA= DN_NWA .* (alpha .* PN_NWA ) + DS_NWA.* (PS_NWA-w_NWA);
    Pi_S_NWA= w_NWA .* DS_NWA;
    Pi_N_NWA= (1 - alpha) .* (DN_NWA .* PN_NWA);

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

    plot(t,Pi_P_NWA+0.*t,'r'); hold on;
    plot(t,Pi_P_PWA,'b'); hold off;

    plot(t,Pi_S_NWA+0.*t,'r'); hold on;
    plot(t,Pi_S_PWA,'b'); hold off;

    %plot(cb,Pi_S_PAA+0*cb,'r'); hold on;
    %plot(cb,Pi_S_PWA+0*cb,'b'); hold off;
  
abs4=Pi_P_PWA - Pi_P_NWA;

    plot(t,abs4,'b'); hold off;


    dif = abs(abs4-0);
    indmu = find(dif == min(dif));
    muet = t(indmu(1));
    minmu1 = [minmu1; muet];




    dif = abs(Pi_S_PWA - Pi_S_NWA);
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
%tmax1111=-(8.*k - 9.*gamma.^2.*k + 2.*gamma.^4.*k - 2.*gamma.*k.*mus - alpha.*gamma.^2.*k + alpha.*gamma.^4.*k + gamma.^3.*k.*mus)./(2.*cb.*gamma - 2.*gamma.*mus - alpha.*gamma.^2 + alpha.*gamma.^4 - 2.*cb.*gamma.^2 - cb.*gamma.^3 + cb.*gamma.^4 + gamma.^3.*mus - 9.*gamma.^2 + 2.*gamma.^4 + 8);

%plot(mus,tmax1111,'-k','MarkerSize',6,'LineWidth',3);hold on;
hold off;

xlabel('mu');
ylabel('t');
%xlabel('二手正品的质量{\it{\mu}}');
%ylabel('区块链的实施成本{\it{c_b}}');