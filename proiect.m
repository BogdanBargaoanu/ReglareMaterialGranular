%% Calcule initiale

n=5.03;


%amplificator de putere
Kap=20+0.25*n;
tap=10^-2; %secunde,10msec

%motor de antrenare

K1=0.25+(10^-2)*n; %N.m/U
K2=5+0.1*n; % rad/sec/N.m
TM1=0.05+(2*10^-3)*n; % sec
TM2=0.5+(10^-2)*n; %sec
Hmot1 = tf(K1,[TM1 1]);
Hmot2 = tf(K2,[TM2 1]);

%Tahogeneratorul de masurare a turatiei TG
KT=0.1; % U/rad/sec
TT=10^-2; % sec
Ht = tf(KT,[TT 1]);


%Transportul melcat TM
TQ1 = 5;% sec
TTM=5;
KTM = 0.12+(10^-2)*n; %(kg/sec)/(rad/sec)
TB=60; %sec
Htm1=tf(0.05,[TB*TQ1 TB 0]);
Htm2=tf(KTM,[TTM 1]);

%Transportatorul cu cupe TC
K=0.9;
T=10; % sec
omega_punct_TC=1;
taum=60/omega_punct_TC;

%Doza gravimetrica cu adaptor
KG=0.16; %mA/kg/sec
Tg=2; %sec
delay = 1/(1+tap);

%Ventil pneumatic si convertorul electropneumatic
Kpg=10^-2;
KCE_KV = 0.025;
TV=4; %sec


%Cuptorul C
KC=200+2*n;
K_OmegaT=0.8;
k_OmegaZ=0.3;
T_OmegaZ=120;%sec
T_C=600+5*n; % sec
TauC=0.15*T_C; % sec
T_OmegaT=100+2*n; % sec
%TauT=0.05*T?

%Traductoarele de temperatura

%pentru pirometru
K_OmegaM=0.16; %mA/C
T_OmegaM=4; %sec

%pentru termorezistenta
K_OmegaC=0.1;%mA/C
T_OmegaC=16;% sec


% *
TMStar = TM2 + TT;
KMStar = K2*KT;

bloc1 = tf(Kap,[tap 1]);
bloc2 = tf(K1,[TM1 1]);
bloc3 = tf(K2,[TM2 1]);
bloc4 = tf(KT,[TT 1]);
bloc5 = tf([TT 1],KT);
Hf = tf(Kap*K1*KMStar,conv([tap 1],conv([TM1 1],[TMStar 1])));

%% Performante GT Clasic

sr = 0.10;
zeta = abs(log(sr))/sqrt(log(sr)^2 + pi^2); %zeta= 0.6 ts=4.746 pentru
ts= 5;
wn = 4/(zeta*ts);
cv = wn/(2*zeta);
lbanda = wn*sqrt(1 - 2*zeta^2 + sqrt(2 - 4*zeta^2 + 4*zeta^4));
Ho1 = tf(wn^2,[1 2*zeta*wn wn^2]);
Hr = minreal(zpk(1/Hf * Ho1/(1-Ho1)));
figure;
step((Hr*Hf/(1+Hr*Hf)));
t = 0:0.1:30;
figure;
lsim(Hr*Hf/(1+Hr*Hf),t,t);
Hd = tf(wn^2,[1 2*zeta*wn 0]);

%% simplificare 1 GT

Hr1 = minreal(zpk(0.00017543 * tf(conv([1 100],[1 16.65]),[1 0])));
Ho1 = Hr1*Hf/(1+Hr1*Hf);
figure;
step(Ho1);
t = 0:0.1:30;
figure;
lsim(Ho1,t,t);

%% simplificare 2 GT

t1 = 0.00017543 * 1/100;
t2 = 0.00017543 * 1/16.65;

% t2 >> t1

Hr2 = tf([t2+t1 1],[1 0]);
Ho2 = Hr2*Hf/(1+Hr2*Hf);
figure;
step(Ho2);
t = 0:0.1:30;
figure;
lsim(Ho2,t,t);


%% Corectie pol-dipol
srP = 0.12;
ts = 2;
zeta = abs(log(srP))/sqrt(log(srP)^2 + pi^2);
wn = 4/(zeta*ts);
cv = wn/(2*zeta);
lbanda = wn*sqrt(1 - 2*zeta^2 + sqrt(2 - 4*zeta^2 + 4*zeta^4));
Ho1 = tf(wn^2,[1 2*zeta*wn wn^2]);
Hr = minreal(zpk(1/Hf * Ho1/(1-Ho1)));
figure;
step((Hr*Hf/(1+Hr*Hf)));
t = 0:0.1:30;
figure;
lsim(Hr*Hf/(1+Hr*Hf),t,t);
Hd = tf(wn^2,[1 2*zeta*wn 0]);

cvStar = 1/4;
src = 0.01;
pc_zc = 1 + src;
pc = src/(1/cv - cvStar);
zc = pc/(1 + src);

Ho2 =tf(conv(pc,conv(wn^2,[1 zc])),conv(zc,conv([1 pc],[1 2*zeta*wn wn^2])));
Hr2 = minreal(zpk(1/Hf * Ho2/(1-Ho2)));
figure;
step((Hr2*Hf/(1+Hr2*Hf)));
t = 0:0.1:30;
figure;
lsim(Hr2*Hf/(1+Hr2*Hf),t,t);


%% simplificare 1 GT corectie

Hr2prim = minreal(zpk(Hr2 * tf([1 0.1259],[1 0.1573])));
figure;
step((Hr2prim*Hf/(1+Hr2prim*Hf)));
t = 0:0.1:30;
figure;
lsim(Hr2prim*Hf/(1+Hr2prim*Hf),t,t);

%% simplificare 2 GT corectie

t1 = 0.00124 * 1/100;
t2 = 0.00124 * 1/16.65;
t3 = 1/4.033;
t4 = 0.00124 * 1/1.785;

%t3 >> t1
Hr2secund = 0.00124 * tf(conv([1 16.65],[1 1.785]),[t3*(t3-t1) t3 0]); 
figure;
step((Hr2secund*Hf/(1+Hr2secund*Hf)));
t = 0:0.1:30;
figure;
lsim(Hr2secund*Hf/(1+Hr2secund*Hf),t,t);
%% Metoda frecventiala
t = 0:0.1:50;
TF=tap+TM1;
KF=Kap*K1;
HMF=tf(KF,[TF 1]);
TMF=TMStar/KMStar; %ts
HMF2=tf(1,[TMF 0]);
%Hfasd=HMF*HMF2;

%% 2.2 regulator P
srP = 0.1;
zetaP = abs(log(srP))/sqrt(log(srP)^2 + pi^2);
trP = 2;
wnP = 4/(zetaP*trP);
cvP = wnP/(2*zetaP);
lbandaP = wnP*sqrt(1 - 2*zetaP^2 + sqrt(2 - 4*zetaP^2 + 4*zetaP^4));
F = -26.4; %de pe bode, db la wf (aprox la wf)
A = 1/(4*zetaP^2*sqrt(2));
N = 20*log10(A);
FN = abs(F-N);
Kp = 10^(FN/20);
bode(Kp * Hf);
wt = 31.7;
HoP = Kp*Hf/(1 + Kp*Hf);
figure;
step(HoP);
figure;
lsim(HoP,t,t);

%% 2.3 regulator PI

srPI = 0.1;
zetaPI = abs(log(srPI))/sqrt(log(srPI)^2 + pi^2);
trPI = 2;
wnPI = 4/(zetaPI*trPI);
cvPI = wnPI/(2*zetaPI);
lbandaPI = wnPI*sqrt(1 - 2*zetaPI^2 + sqrt(2 - 4*zetaPI^2 + 4*zetaPI^4));
F = -26.4;
A = 1/(4*zetaPI^2*sqrt(2));
N = 20*log10(A);
FN = abs(F-N);
Kp = 10^(FN/20);
%bode(Kp*Hf);
cv2 = 12;
cv1 = 0.3;
wt = 31.7; %bode Kp*Hf
wz = 0.1*wt;
wp = wz*cv1/cv2;
Tz = 1/wz;
Tp = 1/wp;
Kpi = Kp*cv2/cv1;
figure;
Hpi = tf([Kpi*Tz Kpi],[Tp 1]);
bode(Hpi*Hf);
Ho = Hpi*Hf/(1 + Hpi*Hf);
figure;
step(Ho);
figure;
lsim(Ho,t,t);
figure;
bode(Ho);

%% 2.4 regulator PD

srPD = 0.1;
zetaPD = abs(log(srPD))/sqrt(log(srPD)^2 + pi^2);
trPD = 2;
wnPD = 4/(zetaPD*trPD);
cvPD = wnPD/(2*zetaPD);
lbandaPD = wnPD*sqrt(1 - 2*zetaPD^2 + sqrt(2 - 4*zetaPD^2 + 4*zetaPD^4));
F = -26;
A = 1/(4*zetaPD^2*sqrt(2));
N = 20*log10(A);
FN = abs(F-N);
Kp = 10^(FN/20);
bode(Kp * Hf);
wt = 29.7;
ts1 = 2/zetaPD^2/wt;
ts2 = 0.5;
wt2 = wt* ts1/ts2;
Td = 1.785 + 1/16.65 + 1/100; % T din Hf
Tn = ts2/ts1 * Td;
Kpd = Kp * wt2/wt;
%tsfinal = 2/(wt*zeta^2);
Hpd = tf([Kpd*Td Kpd],[Tn 1]);
Ho = Hpd*Hf/(1 + Hpd*Hf);
figure;
step(Ho);
figure;
lsim(Ho,t,t);

%% 2.5 regulator PID

srPID = 0.1;
zetaPID = abs(log(srPID))/sqrt(log(srPID)^2 + pi^2);
trPID = 2;
wnPID = 4/(zetaPID*trPID);
cvPID = wnPID/(2*zetaPID);
lbandaPID = wnPID*sqrt(1 - 2*zetaPID^2 + sqrt(2 - 4*zetaPID^2 + 4*zetaPID^4));
F = -26;
A = 1/(4*zetaPID^2*sqrt(2));
N = 20*log10(A);
FN = abs(F-N);
Kp = 10^(FN/20);
bode(Kp * Hf);
wt = 29.7;
ts1 = 2/zetaPID^2/wt;
ts2 = 0.5;
wt2 = wt* ts1/ts2;
Td = 1.785 + 1/16.65 + 1/100;
Tn = ts2/ts1 * Td;
Kpd = Kp * wt2/wt;
%tsfinal = 2/(wt*zeta^2);
Hpd = tf([Kpd*Td Kpd],[Tn 1]);
Hod = Hpd*Hf/(1 + Hpd*Hf);
%figure;
%step(Hod);
%figure;
%lsim(Hod,t,t);
cv1 = 1/0.2;
cv2 = 10;
Kpid = Kpd * wt2/wt * cv2/cv1;
TnNou = 0.5 * ts2/ts1;
wz = 0.1 * wt2;
Tz = 1/wz;
wp = wz*cv1/cv2;
Tp = 1/wp;
Hpid = tf(conv(Kpid,conv([Td 1],[Tz 1])),conv([TnNou 1],[Tp 1]));
Hopid = Hpid * Hf / (1 + Hpid * Hf);
figure;
step(Hopid);
figure;
lsim(Hopid,t,t);

%% 3 Margine impusa

Hf3=tf(KCE_KV*KC*K_OmegaC,conv(conv([TV 1],[T_C 1]),[T_OmegaC 1]),'IoDelay',TauC);

%% PI
faza_impusa = 50;
phi_i = -180 + 15 + faza_impusa;
w_i = 0.006;
Kf = -17.5;
Ti = 4/w_i;
Kp = 10^(-Kf/20);
Hr = tf([Kp*Ti Kp],[Ti 0]);
figure;
bode(Hr*Hf3)

%% PD 

phi_d = -180;
beta = 0.1;
w_d = 0.0142;
Kf = -25;
Td = 1/(w_d*sqrt(beta));
Kp = sqrt(beta)*10^(-Kf/20);
Hr = tf([Kp*Td Kp],[Td*beta 1]);
figure;
bode(Hr*Hf3);

%% PID
faza_impusa = 55;
phi_pid = -180 + 55;
w0_pid=0.00733;
wprim_pid=0.0147;
raport=w0_pid/wprim_pid;
Kf_w0 = -19.1;
Kf_wprim = -25.2;
Kpid = w0_pid/-Kf_w0; %pt raport din interval w0_pid = 0.023?
T0= 2*pi/w0_pid;
Ti_pid = 1.2*T0;
Td_pid = 0.5*T0;
beta_pid = 0.2;
Hr = Kpid*tf(conv([Td_pid 1],[Ti_pid 1]),conv([beta_pid*Td_pid 1],[Ti_pid 0]));
bode(Hr*Hf3);

%% 4.1 Metoda modulului

T = tap; % Tmin = 1/100;
t = 0:0.1:50;
Hdm = tf(1,[2*T*T 2*T 0]);
Hrm = minreal(zpk(Hdm/Hf));
Hom = Hdm/(1+Hdm);
figure;
subplot(211);
step(Hom);
subplot(212);
lsim(Hom,t,t);

%% 4.2 Metoda simetriei

T = tap; % Tmin = 1/100;
t = 0:0.1:50;
Hds = tf([4*T 1],[8*T^3 8*T*T 0 0]);
Hrs = minreal(zpk(Hds/Hf));
Hrs = series(Hrs,tf(1,conv([1 2.5],[1 1])));
Hos = Hds/(1+Hds);
figure;
subplot(211);
step(Hos);
subplot(212);
lsim(Hos,t,t);


%% 5.0
KTMStar = KTM/KT;
TTMStar = TTM - TT;
KF = KCE_KV * KC * K_OmegaC;
KStar = K * KG;
TStar = T + taum + Tg;
TF = TV + T_C + T_OmegaC;
K_OmegaTStar = K_OmegaT/K_OmegaC;
T_OmegaTStar = T_OmegaT-T_OmegaC;

%% 5.2
Hr5 = tf(conv([TM1 1],[TMStar 1]),2*tap*Kap*K1*KMStar);
%HOmega = tf(1,[2*tap 1]);
HOmega = tf(KTMStar,[TTMStar 1]); 
TTMStar2 = TTMStar + 2*tap;
%HRQ = tf([TStar*KTMStar*KMStar*TStar TStar*KTMStar*KMStar],[2*TTMStar2*TStar]);
HRQ = 6.73*tf(conv([TM1 1],[TMStar 1]),[1 0]);
us = tf([0.01121 0.001121],[0.0001141 0.01161 0.5904 1]);

%% 5.3
wStar = 128.6174*4;
T0 = 2*pi/wStar;
Vrc = 0.9*TauC/(TF*KF);
tauIC = 3.3*TauC;
Hoc = tf([Vrc*KF*2.3*TauC Vrc*KF],[tauIC*TF*TauC tauIC*(TF + TauC) tauIC + Vrc*KF*2.3*TauC Vrc*KF]);
[num,den] = tfdata(Hoc);

tauim = 0.6*T0;
taudm = 0.1*T0;
Vrm = 0.101;
P = tf(Vrm, 1);

I = tf(1, [tauim 0]);

D = tf([taudm 0], 1);

PID = P + I + D;

Hrm = P;

%% 6

timp_mort = 1/(1+tap);
T1_stea = T + Tg;
num = KTM*KStar;
den = [conv(conv([2*tap 1],[TTMStar 1]),[T1_stea 1])];
Hf = tf(num, den, 'IOdelay', timp_mort);
Hf_prim = tf(num, den);


%Tmin = 2*tap + TTMStar;
Tmin = T1_stea;
Tr = 2*Tmin;
TB = 4*Tmin;
B = tf(1 , [4*Tmin*Tmin 4*Tmin 0]);
HR1 = B/Hf_prim;

VR = T1_stea/(4*KTMStar*KStar*Tmin);
taui = T1_stea;

HR = tf([VR*taui VR],[taui 0]);

H0 = tf(1,[4*Tmin*Tmin 4*Tmin 1], 'IOdelay', timp_mort);

tr_prim = 6*Tmin;
tr = 6*Tmin + timp_mort;

figure
bode(Hf_prim)

w_secund = 6.1;
modul_Hf = 83.3;
VR_fin = 1/modul_Hf;
taui_fin = 4/w_secund;

Hr_fin = tf([VR_fin*taui_fin VR_fin],[taui_fin 0]);
figure;
step(Hr_fin*Hf/(1+Hr_fin*Hf));