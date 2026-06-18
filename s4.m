clc;
clear;
syms k gamma alpha q1 qn ps pn PS_PWA PN_PWA q2 mu t cb w_PWA
qc=0;
q1=1;
q2=mu.*q1;
qn=(1-k-t).*qc+(k+t).*q1;
qs=(1-k-t).*qc+(k+t).*q2;
DN_PWA=qn-PN_PWA-gamma.*(qs-PS_PWA);
DS_PWA=qs-PS_PWA-gamma.*(qn-PN_PWA);
Pi_P_PWA= DN_PWA .* (alpha .* PN_PWA-t.*cb ) + DS_PWA.* (PS_PWA-w_PWA-t.*cb);
Pi_S_PWA= w_PWA .* DS_PWA;
Pi_N_PWA= (1 - alpha) .* (DN_PWA .* PN_PWA);

pnmax1= diff(Pi_N_PWA, PN_PWA);
psmax1= diff(Pi_P_PWA, PS_PWA);
[p_nmax,p_smax] = solve(pnmax1, psmax1, PN_PWA,PS_PWA)



clc;
clear;
syms k gamma alpha q1 qn ps pn PS_PWA PN_PWA q2 mu t cb w_PWA
PN_PWA=(gamma^2*k - 2*t - gamma*w_PWA - 2*k + gamma^2*t - cb*gamma*t + gamma*k*mu + gamma*mu*t + cb*gamma^2*t)/(alpha*gamma^2 + gamma^2 - 4); 
PS_PWA=(gamma*k - 2*w_PWA - 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t - alpha*gamma*k - alpha*gamma*t + 2*cb*gamma*t + gamma^2*k*mu + alpha*gamma^2*k*mu + alpha*gamma^2*mu*t)/(alpha*gamma^2 + gamma^2 - 4);
qc=0;
q1=1;
q2=mu.*q1;
qn=(1-k-t).*qc+(k+t).*q1;
qs=(1-k-t).*qc+(k+t).*q2;
DN_PWA=qn-PN_PWA-gamma.*(qs-PS_PWA);
DS_PWA=qs-PS_PWA-gamma.*(qn-PN_PWA);
Pi_P_PWA= DN_PWA .* (alpha .* PN_PWA-t.*cb ) + DS_PWA.* (PS_PWA-w_PWA-t.*cb);
Pi_S_PWA= w_PWA .* DS_PWA;
Pi_N_PWA= (1 - alpha) .* (DN_PWA .* PN_PWA);
wmax1=diff(Pi_S_PWA, w_PWA);
wmax=solve(wmax1,w_PWA)


clc;
clear;
syms k gamma alpha q1 qn ps pn PS_PWA PN_PWA q2 mu t cb w_PWA
w_PWA=(gamma*k + 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t + alpha*gamma*k + alpha*gamma*t - 2*cb*gamma*t - alpha*gamma^3*k - alpha*gamma^3*t - cb*gamma^2*t + cb*gamma^3*t + gamma^2*k*mu)/(2*gamma^2 - 4);
PN_PWA=(gamma^2*k - 2*t - gamma*w_PWA - 2*k + gamma^2*t - cb*gamma*t + gamma*k*mu + gamma*mu*t + cb*gamma^2*t)/(alpha*gamma^2 + gamma^2 - 4); 
PS_PWA=(gamma*k - 2*w_PWA - 2*cb*t - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t - alpha*gamma*k - alpha*gamma*t + 2*cb*gamma*t + gamma^2*k*mu + alpha*gamma^2*k*mu + alpha*gamma^2*mu*t)/(alpha*gamma^2 + gamma^2 - 4);
qc=0;
q1=1;
q2=mu.*q1;
qn=(1-k-t).*qc+(k+t).*q1;
qs=(1-k-t).*qc+(k+t).*q2;
DN_PWA=qn-PN_PWA-gamma.*(qs-PS_PWA)
DS_PWA=qs-PS_PWA-gamma.*(qn-PN_PWA);
Pi_P_PWA= DN_PWA .* (alpha .* PN_PWA-t.*cb ) + DS_PWA.* (PS_PWA-w_PWA-t.*cb);
Pi_S_PWA= w_PWA .* DS_PWA;
Pi_N_PWA= (1 - alpha) .* (DN_PWA .* PN_PWA);

solve(DN_PWA,t)