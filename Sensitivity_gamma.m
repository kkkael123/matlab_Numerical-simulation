%AI技术的成本对于二手产品价格的影响
clc;
clc;
clear;

alpha = 0.2;
gamma = (0:0.01:1)';
q1 = 1;
mu =0.7; % 请根据需要设置 mu 的值
q2 = mu * q1;
k = 0.7;
qs = k * mu;
qn = k * q1;
cb = 0.01;

A = 2 - gamma.^2;
B = 4 - gamma.^2 - alpha * gamma.^2;


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


plot(gamma, Pi_S_PAA, 'k', 'MarkerSize', 5, 'LineWidth', 2); hold on;
plot(gamma, Pi_S_PWA, 'k--', 'MarkerSize', 5, 'LineWidth', 2); hold on;
plot(gamma, Pi_P_PAA, 'b', 'MarkerSize', 5, 'LineWidth', 2.5); hold on;
plot(gamma, Pi_P_PWA, 'b--', 'MarkerSize', 5, 'LineWidth', 2.5); hold off;

%plot(k, Pi_S_PAA+0.*k, 'k', 'MarkerSize', 5, 'LineWidth', 2); hold on;
%plot(k, Pi_S_PWA+0.*k, 'k-^', 'MarkerSize', 5, 'LineWidth', 2); hold on;
%plot(k, Pi_P_PAA+0.*k, 'b', 'MarkerSize', 5, 'LineWidth', 2.5); hold on;
%plot(k, Pi_P_PWA+0.*k, 'b--', 'MarkerSize', 5, 'LineWidth', 2.5); hold off;


xlabel('{\it{\gamma}}');
ylabel('{\it{\Pi}}');
legend('\it{\Pi_S^{BA}}', '\it{\Pi_S^{BR}}', '\it{\Pi_P^{BA}}','\it{\Pi_P^{BR}}');
set(gca, 'LineWidth', 1.2);
set(gca, 'FontSize', 11);
xlim([0 0.8])
ylim([0 0.1])


