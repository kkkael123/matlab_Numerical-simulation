%NA vs BA的对比仿真:平台利润和新品供应商的利润

clc;
clear;

%alpha = 0.8; gamma = 0.2; k = 0.7; 
alpha = 0.7; gamma = 0.2; k = 0.7; 
A = 2 - gamma^2; 
%A/gamma
%A*(1+k)/(2*gamma)
B = 4 - gamma^2 - alpha * gamma^2;
q1 = 1 ; qn = k * q1;
    mu_low=gamma*(alpha*A+1-alpha)./A;
    mus = (mu_low:0.01:1)';
mincb1=[];
mincb2=[];
i=1;
while i<=length(mus)
    mu= mus(i);
    q2 = mu * q1; qs = k * q2; 
    %cb_max1=((A.*A-A).*mu+(-alpha.*A.*A+(2.*alpha-1).*A+(-alpha+1)).*gamma)./((1-gamma.*gamma).*(1-gamma).*(A-1));
    %cb_max2=0;
    %cb_max=max(cb_max1,cb_max2);
    c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));
    cb = (0:0.001:c_bmax)'; 
    
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


    
    plot(cb,Pi_P_PAA+0*cb,'r'); hold on;
    plot(cb,Pi_P_NAA+0*cb,'b'); hold off;


    dif = abs(Pi_P_PAA-Pi_P_NAA);
    indcb = find(dif == min(dif));
    cbet = cb(indcb(1));
    mincb1 = [mincb1; cbet];



    i=i+1;
end
c_bmax1=(2.*mus+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mus))./((1-gamma).^2.*(1 + gamma));
c_bmax2=max(0,c_bmax1);
plot(mus,mincb1,'-r','MarkerSize',6,'LineWidth',3);hold on;
%plot(mus,mincb2,'--b','MarkerSize',6,'LineWidth',3);hold on;
plot(mus,c_bmax2,'k','MarkerSize',6,'LineWidth',3);hold on;
hold off;

%xlabel('产品保值率{\it{\mu}}');
%ylabel('AI技术的实施成本{\it{c_b}}');