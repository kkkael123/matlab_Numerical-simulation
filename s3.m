clc;
clear;
syms k gamma alpha q1 qn ps pn PS_PAA PN_PAA q2 mu t cb
qc=0;
q1=1;
q2=mu.*q1;
qn=(1-k-t).*qc+(k+t).*q1;
qs=(1-k-t).*qc+(k+t).*q2;
DN_PAA=qn-PN_PAA-gamma.*(qs-PS_PAA);
DS_PAA=qs-PS_PAA-gamma.*(qn-PN_PAA);
Pi_P_PAA=DN_PAA .* (alpha .* PN_PAA-t.*cb) + DS_PAA .* (alpha .* PS_PAA-t.*cb);
Pi_S_PAA= (1 - alpha) .* (DS_PAA .* PS_PAA);
Pi_N_PAA= (1 - alpha) .* (DN_PAA .* PN_PAA);

pnmax1= diff(Pi_N_PAA, PN_PAA);
psmax1= diff(Pi_S_PAA, PS_PAA);
[p_nmax,p_smax] = solve(pnmax1, psmax1, PN_PAA,PS_PAA)

clc;
clear;
syms k gamma alpha q1 qn ps pn PS_PAA PN_PAA q2 mu t cb
PN_PAA=(gamma^2*k - 2*t - 2*k + gamma^2*t + gamma*k*mu + gamma*mu*t)/(gamma^2 - 4);
PS_PAA=(gamma*k - 2*k*mu + gamma*t - 2*mu*t + gamma^2*mu*t + gamma^2*k*mu)/(gamma^2 - 4);
qc=0;
q1=1;
q2=mu.*q1;
qn=(1-k-t).*qc+(k+t).*q1;
qs=(1-k-t).*qc+(k+t).*q2;
DN_PAA=qn-PN_PAA-gamma.*(qs-PS_PAA);
DS_PAA=qs-PS_PAA-gamma.*(qn-PN_PAA);
Pi_P_PAA=DN_PAA .* (alpha .* PN_PAA-t.*cb) + DS_PAA .* (alpha .* PS_PAA-t.*cb);
Pi_S_PAA= (1 - alpha) .* (DS_PAA .* PS_PAA);
Pi_N_PAA= (1 - alpha) .* (DN_PAA .* PN_PAA);