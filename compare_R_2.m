%第四章n和p的对比仿真:平台利润和新品供应商的利润
clc;
clear;

%alpha = 0.8; gamma = 0.2; k = 0.7; 
alpha = 0.8; gamma = 0.4; k = 0.3; 
A = 2 - gamma^2; 
t=0.5;
qc=0;
%A/gamma
%A*(1+k)/(2*gamma)
B = 4 - gamma^2 - alpha * gamma^2;
q1 = 1 ; qn = k * q1;
    mu_low=gamma*(alpha*A+1-alpha)./A
    mus = (mu_low:0.01:1)';
mincb1=[];
mincb2=[];
i=1;
while i<=length(mus)
    mu= mus(i);

    q2=mu.*q1;
qn = k * q1;
qs=  k*q2; 
    qn_A=(1-k-t).*qc+(k+t).*q1;
    qs_A=(1-k-t).*qc+(k+t).*q2;
    %cb_max1=((A.*A-A).*mu+(-alpha.*A.*A+(2.*alpha-1).*A+(-alpha+1)).*gamma)./((1-gamma.*gamma).*(1-gamma).*(A-1));
    %cb_max2=0;
    %cb_max=max(cb_max1,cb_max2);
    c_bmax=(2.*mu+gamma.*(-1 + alpha .*(-1 +gamma.^2) - gamma.*mu))./((1-gamma).^2.*(1 + gamma));
    cb = (0:0.001:c_bmax)'; 
    

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

    % P-WA 平台承担AI技术-批发模式
w_PWA=(gamma*k + 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t + alpha*gamma*k + alpha*gamma*t - 2*cb*gamma*t - alpha*gamma^3*k - alpha*gamma^3*t - cb*gamma^2*t + cb*gamma^3*t + gamma^2*k*mu)./(2*gamma^2 - 4);
PN_PWA=(gamma^2*k - 2*t - gamma*w_PWA - 2*k + gamma^2*t - cb*gamma*t + gamma*k*mu + gamma*mu*t + cb*gamma^2*t)./(alpha*gamma^2 + gamma^2 - 4); 
PS_PWA=(gamma*k - 2*w_PWA - 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t - alpha*gamma*k - alpha*gamma*t + 2*cb*gamma*t + gamma^2*k*mu + alpha*gamma^2*k*mu + alpha*gamma^2*mu*t)./(alpha*gamma^2 + gamma^2 - 4);

DN_PWA=qn-PN_PWA-gamma.*(qs-PS_PWA);
DS_PWA=qs-PS_PWA-gamma.*(qn-PN_PWA);
Pi_P_PWA= DN_PWA .* (alpha .* PN_PWA-t.*cb ) + DS_PWA.* (PS_PWA-w_PWA-t.*cb);
Pi_S_PWA= w_PWA .* DS_PWA;
Pi_N_PWA= (1 - alpha) .* (DN_PWA .* PN_PWA);

    plot(cb,Pi_P_PWA+0*cb,'r'); hold on;
    plot(cb,Pi_P_NWA+0*cb,'b'); hold off;

    plot(cb,Pi_N_PWA+0*cb,'r'); hold on;
    plot(cb,Pi_N_NWA+0*cb,'b'); hold off;

    plot(cb,Pi_N_PWA-Pi_N_NWA);hold off;

    dif = abs(Pi_P_PWA-Pi_P_NWA);
    indcb = find(dif == min(dif));
    cbet = cb(indcb(1));
    mincb1 = [mincb1; cbet];

    dif = abs(Pi_N_PWA-Pi_N_NWA);
    indcb = find(dif == min(dif));
    cbet = cb(indcb(1));
    mincb2 = [mincb2; cbet];

    i=i+1;
end
plot(mus,mincb1,'-r','MarkerSize',6,'LineWidth',3);hold on;
plot(mus,mincb2,'--b','MarkerSize',6,'LineWidth',3);hold on;

hold off;

%xlabel('产品保值率{\it{\mu}}');
%ylabel('AI技术的实施成本{\it{c_b}}');