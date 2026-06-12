%第三章PAA与PWA的对比
clc;
clear;
%gamma = 0.2; 
%alpha = 0.2;
k = 0.3; 
gamma = 0.2; 
cb=0.3;

q1 = 1; qn = k * q1;
alphas = (0.001:0.01:1)';
minmu1=[];
minmu2=[];
i=1;
while i<=length(alphas)
    alpha= alphas(i);
    
    %c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));
A = 2 - gamma^2;
B = 4 - gamma^2 - alpha .* gamma^2;
mu_low1=gamma*(alpha*A+1-alpha)./A;
mu_low2=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alpha.*(1-gamma.^2)))./(2-gamma.^2);
mu_low=max(mu_low1,mu_low2);
mu = (mu_low:0.001:1)';
q2 = mu * q1; qs = k * mu; 
    % P-AA 平台承担AI技术-代理模式
    % 均衡解
    PS_PAA = (A .* q2 - gamma .* q1) ./ (2 + A);
    PN_PAA = (A .* q1 - gamma .* q2) ./ (2 + A);
    % 需求函数
    DS_PAA = (A .* q2 - gamma .* q1) ./ (2 + A);
    DN_PAA = (A .* q1 - gamma .* q2) ./ (2 + A);
    %利润函数
    Pi_P_PAA=DN_PAA .* (alpha .* PN_PAA-cb) + DS_PAA .* (alpha .* PS_PAA-cb);
    Pi_S_PAA= (1 - alpha) .* (DS_PAA .* PS_PAA);
    Pi_N_PAA= (1 - alpha) .* (DN_PAA .* PN_PAA);


    % P-WA 平台承担AI技术-批发模式
    % 均衡解
    w_PWA = (A.*q2+(-alpha.*A+(alpha-1)).*gamma.*q1+(1-A).*(1-gamma).*cb)./(2.*(A-1));
    PS_PWA =(((alpha-2).*A+(1-alpha)).*gamma.*q1+(2.*(alpha+1).*A.^2-(6.*alpha+1).*A+4.*alpha).*q2+(1-gamma).*(A-1).*cb)./(2.*(A-1).*B);
    PN_PWA =(((alpha+2).*A.^2+(-3.*alpha-1).*A+2.*alpha-2).*q1+(2-A).*gamma.*q2+(1-gamma).*(A-1).*gamma.*cb)./(2.*(A-1).*B);
    % 需求函数
    DS_PWA =((A.*A-A).*q2+(-alpha.*A.*A+(2.*alpha-1).*A+(-alpha+1)).*gamma.*q1+(gamma.*gamma-1).*(1-gamma).*(A-1).*cb)./(2.*(A-1).*B);
    DN_PWA =((2.*A.^2+(-2).*A).*q1+(-2.*A+2).*gamma.*q2)./(2.*(A-1).*B);
    %利润函数
    Pi_P_PWA= DN_PWA .* (alpha .* PN_PWA-cb ) + DS_PWA.* (PS_PWA-w_PWA-cb);
    Pi_S_PWA= w_PWA .* DS_PWA;
    Pi_N_PWA= (1 - alpha) .* (DN_PWA .* PN_PWA);

    %plot(mu,Pi_P_NAA+0*mu,'r'); hold on;
    %plot(mu,Pi_P_NWA+0*mu,'b'); hold off;

    plot(mu,Pi_P_PAA+0*mu,'r'); hold on;
    plot(mu,Pi_P_PWA+0*mu,'b'); hold off;

    plot(mu,Pi_S_PAA+0*mu,'r'); hold on;
    plot(mu,Pi_S_PWA+0*mu,'b'); hold off;

    %plot(cb,Pi_S_PAA+0*cb,'r'); hold on;
    %plot(cb,Pi_S_PWA+0*cb,'b'); hold off;


    
    dif = abs(Pi_P_PWA - Pi_P_PAA);
    indmu = find(dif == min(dif));
    muet = mu(indmu(1));
    minmu1 = [minmu1; muet];

    dif = abs(Pi_S_PWA - Pi_S_PAA);
    indmu = find(dif == min(dif));
    muet = mu(indmu(1));
    minmu2 = [minmu2; muet];
    

    i=i+1;
end
mu_low11=gamma*(alphas*A+1-alphas)./A;
mu_low22=(cb.*(1-gamma).^2.*(1+gamma)+gamma.*(1+alphas.*(1-gamma.^2)))./(2-gamma.^2);
mu_low33=max(mu_low11,mu_low22);
plot(alphas,minmu1,'-r','MarkerSize',6,'LineWidth',3);hold on;
plot(alphas,minmu2,'-b','MarkerSize',6,'LineWidth',3);hold on;
plot(alphas,mu_low33,'-k','MarkerSize',6,'LineWidth',3);hold on;
hold off;



%xlabel('二手正品的质量{\it{\mu}}');
%ylabel('区块链的实施成本{\it{c_b}}');