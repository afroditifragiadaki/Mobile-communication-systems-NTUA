clear all; close all;
am=18042;
am=mod(am,10);
Gt = 10+am; % dBi
Gr = 30+am; % dBi
T = 290+am; % K
B = (20+am)*(10^6); % Hz
Pt = -(10+am); % dBW
fc = (10+am); % GHz
lat = 30+am; % degrees North
long = 20+am; % degrees West
%Kata tin ektelesi tou kwdika prepei stin arxi
%tou section me to wrwtima 2 na epilexthoun ta
%swsta fc gia to analogo am kai na mpoun se 
%sxolio ta alla fc

%% erwtima 1
ath_lat = 37.98; ath_long = 23.73;
ioa_lat = 39.66; ioa_long = 20.85;
p=0:0.0001:1;
[rr,p0] = itur_p837_5(p,lat,long);
[ath_rr,Ath_p0] = itur_p837_5(p,ath_lat,ath_long);
[ioa_rr,Ioa_p0] = itur_p837_5(p,ioa_lat,ioa_long);
figure(1);
plot(rr,p,"red");
hold on
plot(ath_rr,p,"blue");
hold on
plot(ioa_rr,p,"green");
hold off
grid on
lgd = sprintf('%.1f North %.2f West',lat, long);
legend(lgd,"Athens","Ioannina")
xlabel("Rainfall Rate exceeded for p% of an average year (mm/h)");
ylabel("Outage Prob");
title("Erwtima 1");

%% erwtima 2

d = [4000,8000,12000,15000]; % miki zeuksis

%am 18115
%kh=0.04481; ah=1.1233; kv=0.05008; av=1.0440; %fc=15 
%kh2=0.2403; ah2=0.9485; kv2=0.2291; av2=0.9129; %2fc=30
%kh25=0.3895; ah25=0.8853; kv25=0.37385; av25=0.85865; %2.5fc=37.5
%kh3=0.5521; ah3=0.8355; kv3=0.5375; av3=0.8123; %3fc=45

%am 18006
%kh=0.05282; ah=1.1086; kv=0.05899; av=1.0273; %fc=16
%kh2=0.2778; ah2=0.9302; kv2=0.2646; av2=0.8981; %2fc=32
%kh25=0.4431; ah25=0.8673; kv25=0.4274; av25=0.8421; %2.5fc=40
%kh3=0.6172; ah3=0.8187; kv3=0.6037; av3=0.7967 ; %3fc=48


%am 18042
kh=0.02386; ah=1.1825; kv=0.02455; av=1.1216; %fc=12
kh2=0.1425 ; ah2=1.0101; kv2=0.1404; av2=0.9561; %2fc=24
kh25=0.2403; ah25=0.9485; kv25=0.2291; av25=0.9129; %2.5fc=30
kh3=0.3580; ah3=0.8967; kv3=0.3427; av3=0.8690; %3fc=36


%We assume linear or circular polarization, so t=45 degrees
k=(kh+kv)/2;
a=(kh*ah+kv*av)/(2*k);
g_R=k*rr(2)^a; %rr(2) gia prob 0.0001 / 0.01%
r=1./(0.477*d.^0.633*rr(2)^(0.733*a)*(fc)^0.123 - 10.579*(1-exp(-0.024*d)));
C0=0.12+0.4*(log10((fc/10)^0.8));
C1=(0.07^C0)*0.12^(1-C0);
C2=0.855*C0+0.546*(1-C0);
C3=0.139*C0+0.043*(1-C0);
A_001=g_R.*d.*r;
Ap4=A_001(1)+10*log10(C1*p.^(-C2-C3*log10(p)));
Ap8=A_001(2)+10*log10(C1*p.^(-C2-C3*log10(p)));
Ap12=A_001(3)+10*log10(C1*p.^(-C2-C3*log10(p)));
Ap15=A_001(4)+10*log10(C1*p.^(-C2-C3*log10(p)));

figure(2);
semilogy(Ap4,p,'black');
hold on
semilogy(Ap8,p,'blue');
hold on
semilogy(Ap12,p,'red');
hold on
semilogy(Ap15,p,'green');
hold off
grid on
legend("d=4km","d=8km","d=12km","d=15km")
xlabel("Rain Attenuation(dB)");
ylabel("Exceedance Probability");
title("Erwtima 2 fc");

%Gia 2.5fc
fc25=2.5*fc;

k_25=(kh25+kv25)/2;
a_25=(kh25*ah25+kv25*av25)/(2*k_25);
g_R_25=k_25*rr(2)^a_25; %rr(2) gia prob 0.0001 / 0.01%
r_25=1./(0.477*d.^0.633*rr(2)^(0.733*a_25)*(fc25)^0.123 - 10.579*(1-exp(-0.024*d)));
C0_25=0.12+0.4*(log10((fc25/10)^0.8));
C1_25=(0.07^C0_25)*0.12^(1-C0_25);
C2_25=0.855*C0_25+0.546*(1-C0_25);
C3_25=0.139*C0_25+0.043*(1-C0_25);
A_001_25=g_R_25.*d.*r_25;
Ap4_25=A_001_25(1)+10*log10(C1_25*p.^(-C2_25-C3_25*log10(p)));
Ap8_25=A_001_25(2)+10*log10(C1_25*p.^(-C2_25-C3_25*log10(p)));
Ap12_25=A_001_25(3)+10*log10(C1_25*p.^(-C2_25-C3_25*log10(p)));
Ap15_25=A_001_25(4)+10*log10(C1_25*p.^(-C2_25-C3_25*log10(p)));

figure(3);
semilogy(Ap4_25,p,'black');
hold on
semilogy(Ap8_25,p,'blue');
hold on
semilogy(Ap12_25,p,'red');
hold on
semilogy(Ap15_25,p,'green');
hold off
grid on
legend("d=4km","d=8km","d=12km","d=15km")
xlabel("Rain Attenuation(dB)");
ylabel("Exceedance Probability");
title("Erwtima 2 2.5fc");


%% Erwtima 3


%Gia 2fc
fc2=2*fc;

k_2=(kh2+kv2)/2;
a_2=(kh2*ah2+kv2*av2)/(2*k_2);
g_R_2=k_2*rr(2)^a_2; %rr(2) gia prob 0.0001 / 0.01%
r_2=1./(0.477*d.^0.633*rr(2)^(0.733*a_2)*(fc2)^0.123 - 10.579*(1-exp(-0.024*d)));
C0_2=0.12+0.4*(log10((fc2/10)^0.8));
C1_2=(0.07^C0_2)*0.12^(1-C0_2);
C2_2=0.855*C0_2+0.546*(1-C0_2);
C3_2=0.139*C0_2+0.043*(1-C0_2);
A_001_2=g_R_2.*d.*r_2;
Ap4_2=A_001_2(1)+10*log10(C1_2*p.^(-C2_2-C3_2*log10(p)));
Ap8_2=A_001_2(2)+10*log10(C1_2*p.^(-C2_2-C3_2*log10(p)));
Ap12_2=A_001_2(3)+10*log10(C1_2*p.^(-C2_2-C3_2*log10(p)));
Ap15_2=A_001_2(4)+10*log10(C1_2*p.^(-C2_2-C3_2*log10(p)));



%Gia 3fc
fc3=3*fc;

k_3=(kh3+kv3)/2;
a_3=(kh3*ah3+kv3*av3)/(2*k_3);
g_R_3=k_3*rr(2)^a_3; %rr(2) gia prob 0.0001 / 0.01%
r_3=1./(0.477*d.^0.633*rr(2)^(0.733*a_3)*(fc3)^0.123 - 10.579*(1-exp(-0.024*d)));
C0_3=0.12+0.4*(log10((fc3/10)^0.8));
C1_3=(0.07^C0_3)*0.12^(1-C0_3);
C2_3=0.855*C0_3+0.546*(1-C0_3);
C3_3=0.139*C0_3+0.043*(1-C0_3);
A_001_3=g_R_3.*d.*r_3;
Ap4_3=A_001_3(1)+10*log10(C1_3*p.^(-C2_3-C3_3*log10(p)));
Ap8_3=A_001_3(2)+10*log10(C1_3*p.^(-C2_3-C3_3*log10(p)));
Ap12_3=A_001_3(3)+10*log10(C1_3*p.^(-C2_3-C3_3*log10(p)));
Ap15_3=A_001_3(4)+10*log10(C1_3*p.^(-C2_3-C3_3*log10(p)));



f = [1, 2, 3]*fc*10^9;
[F,D] = meshgrid(f,d);
lambda = 3*(10^8)./F; 
k_con=1.38*(10^(-23));
FSL = 20*log10(4*pi*D./lambda);
N = 10*log10(k_con*T*B);
SNRcs = Pt + Gt + Gr - FSL - N;

SNRth_4= SNRcs(1,1) - Ap4;
SNRth_8= SNRcs(2,1) - Ap8;
SNRth_12= SNRcs(3,1) - Ap12;
SNRth_15= SNRcs(4,1) - Ap15;
SNRth_4_2= SNRcs(1,2) - Ap4_2;
SNRth_8_2= SNRcs(2,2) - Ap8_2;
SNRth_12_2= SNRcs(3,2) - Ap12_2;
SNRth_15_2= SNRcs(4,2) - Ap15_2;
SNRth_4_3 = SNRcs(1,3) - Ap4_3;
SNRth_8_3= SNRcs(2,3) - Ap8_3;
SNRth_12_3= SNRcs(3,3) - Ap12_3;
SNRth_15_3= SNRcs(4,3) - Ap15_3;



figure(4);
semilogy(SNRth_4,p,'black');
hold on
semilogy(SNRth_8,p,'blue');
hold on
semilogy(SNRth_12,p,'red');
hold on
semilogy(SNRth_15,p,'green');
hold off
grid on
legend("d=4km","d=8km","d=12km","d=15km")
xlabel("SNRth (dB)");
ylabel("Outage Probability");
title("Erwtima 3 fc");

figure(5);
semilogy(SNRth_4_2,p,'black');
hold on
semilogy(SNRth_8_2,p,'blue');
hold on
semilogy(SNRth_12_2,p,'red');
hold on
semilogy(SNRth_15_2,p,'green');
hold off
grid on
legend("d=4km","d=8km","d=12km","d=15km")
xlabel("SNRth (dB)");
ylabel("Outage Probability");
title("Erwtima 3 2fc");

figure(6);
semilogy(SNRth_4_3,p,'black');
hold on
semilogy(SNRth_8_3,p,'blue');
hold on
semilogy(SNRth_12_3,p,'red');
hold on
semilogy(SNRth_15_3,p,'green');
hold off
grid on
legend("d=4km","d=8km","d=12km","d=15km")
xlabel("SNRth (dB)");
ylabel("Outage Probability");
title("Erwtima 3 3fc");