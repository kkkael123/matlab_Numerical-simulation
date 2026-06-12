%第五章N情景两种模式的对比：平台和二手供应商利润
clc;
clear;

%k = 0.7; gamma=0.3;
k = 0.3; alpha=0.6;

q1 = 1; qn = k * q1;
gammas = (0:0.01:1)';
minmu1=[];
minmu2=[];
i=1;
while i<=length(gammas)
    gamma= gammas(i);
    A = 2 - gamma^2;
    B = 4 - gamma^2 - alpha * gamma^2;
    mu_low=gamma*(alpha*A+1-alpha)./A;
    mu = (mu_low:0.001:1)';
    q2 = mu * q1; qs = k * q2; 

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

    % N-WA 无AI技术-批发模式
    % 均衡解
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

     plot(mu,Pi_P_NAA+0*mu,'r'); hold on;
     plot(mu,Pi_P_NWA+0*mu,'b'); hold off;
       
    plot(mu,Pi_S_NAA+0*mu,'r'); hold on;
    plot(mu,Pi_S_NWA+0*mu,'b'); hold off;
    Pi_P_NWA - Pi_P_NAA;
    Pi_S_NWA - Pi_S_NAA;
    dif = abs(Pi_P_NWA - Pi_P_NAA);
    indmu = find(dif == min(dif));
    muet = mu(indmu(1));
    minmu1 = [minmu1; muet];

    dif = abs(Pi_S_NWA - Pi_S_NAA);
    indmu = find(dif == min(dif));
    muet = mu(indmu(1));
    minmu2 = [minmu2; muet];

    i=i+1;
end
mu_low2=gammas.*(alpha.*(2-gammas.^2)+1-alpha)./(2-gammas.^2);
plot(gammas,minmu1,'-r','MarkerSize',6,'LineWidth',3);hold on;
plot(gammas,minmu2,'-b','MarkerSize',6,'LineWidth',3);hold on;
plot(gammas,mu_low2,'-k','MarkerSize',6,'LineWidth',3);hold on;

hold off;

%xlabel('二手产品保值率{\it{\mu}}');
%ylabel('区块链的实施成本{\it{c_b}}');