clear all; close all;
L_max = 142; 
h_MS = 1.42; % Mobile Station 
h_BS = 30:0.01:200 ; % Base Station
Pthr = -142; 
fc = 1.85*(10^9); 

%% (A)
a_urban = 3.2*(log10(11.75*h_MS))^2 - 4.97; 
a_rural = (1.1*log10(fc/(10^6)) - 0.7)*h_MS - (1.56*log10(fc/(10^6)) - 0.8);
A_urban = 46.3 + 33.9*log10(fc/(10^6)) - 13.82*log10(h_BS) - a_urban;
A_rural = 46.3 + 33.9*log10(fc/(10^6)) - 13.82*log10(h_BS) - a_rural;
B = 44.9 - 6.55*log10(h_BS);
C_urban = 3; 
C_rural = 0; 

d_urban = 10.^((L_max - A_urban - C_urban)./B);
d_rural = 10.^((L_max - A_rural - C_rural)./B);

figure(1);
plot(h_BS,d_urban,'red');
hold on
plot(h_BS,d_rural,'green');
hold off
grid on
legend("Urban","Rural");
xlabel("Height of Base Station (m)");
ylabel("Maximum Distance of Coverage (km)")
title("(a)");

%% (B)

lamda = 300./(fc+10^(-6)); % wavelength (m)
Dh = 10; % interdecile height variation (m)
if (Dh/lamda) < 3000
    sigma_rural = 6 + 0.69*sqrt(Dh/lamda) - 0.0063*(Dh/lamda); 
else
    sigma_rural = 25;
end
sigma_urban = 5.25 + 0.42*log10((fc/(10^6))/100) + 1.01*(log10((fc/(10^6))/100))^2; 

A_urban_30 = 46.3 + 33.9*log10(fc/(10^6)) - 13.82*log10(30) - a_urban;
A_urban_200 = 46.3 + 33.9*log10(fc/(10^6)) - 13.82*log10(200) - a_urban;
A_rural_30 = 46.3 + 33.9*log10(fc/(10^6)) - 13.82*log10(30) - a_rural;
A_rural_200 = 46.3 + 33.9*log10(fc/(10^6)) - 13.82*log10(200) - a_rural;
B_30 = 44.9 - 6.55*log10(30);
B_200 = 44.9 - 6.55*log10(200);

P_coverage = 0:0.001:1; 

L_urban = L_max - sqrt(2)*sigma_urban*erfinv(2*P_coverage - 1);
L_rural = L_max - sqrt(2)*sigma_rural*erfinv(2*P_coverage - 1);
d_urban_30 = 10.^((L_urban - A_urban_30 - C_urban)./B_30);
d_urban_200 = 10.^((L_urban - A_urban_200 - C_urban)./B_200);
d_rural_30 = 10.^((L_rural - A_rural_30 - C_rural)./B_30);
d_rural_200 = 10.^((L_rural - A_rural_200 - C_rural)./B_200);

figure(2);
plot(P_coverage,d_urban_30,'red');
hold on
plot(P_coverage,d_urban_200,'blue');
hold off
grid on
legend("Urban h_BS = 30m","Urban h_BS = 200m");
xlabel("Probability of Coverage Pr[L<Lmax]");
ylabel("Maximum Distance of Coverage (km)")
title("URBAN");

figure(3)
plot(P_coverage,d_rural_30);
hold on
plot(P_coverage,d_rural_200);
hold off
grid on
legend("Rural h_BS = 30m","Rural h_BS = 200m");
xlabel("Probability of Coverage Pr[L<Lmax]");
ylabel("Maximum Distance of Coverage (km)")
title("RURAL");

%% (Γ)
Gr = 2;  
Gt = 15; 
d= 1:0.001:20; 

L_urban_30 = A_urban_30 + B_30*log(d) + C_urban;
L_urban_200 = A_urban_200 + B_200*log(d) + C_urban;
L_rural_30 = A_rural_30 + B_30*log(d) + C_rural;
L_rural_200 = A_rural_200 + B_200*log(d) + C_rural;

%urban_30
Pt_urban_30_95 = Pthr - Gt - Gr + L_urban_30 + sqrt(2)*sigma_urban*erfinv(2*0.95 - 1);
Pt_urban_30_97 = Pthr - Gt - Gr + L_urban_30 + sqrt(2)*sigma_urban*erfinv(2*0.97 - 1);
Pt_urban_30_99 = Pthr - Gt - Gr + L_urban_30 + sqrt(2)*sigma_urban*erfinv(2*0.99 - 1);
Pt_urban_30_995 = Pthr - Gt - Gr + L_urban_30 + sqrt(2)*sigma_urban*erfinv(2*0.995 - 1);

% urban_200
Pt_urban_200_95 = Pthr - Gt - Gr + L_urban_200 + sqrt(2)*sigma_urban*erfinv(2*0.95 - 1);
Pt_urban_200_97 = Pthr - Gt - Gr + L_urban_200 + sqrt(2)*sigma_urban*erfinv(2*0.97 - 1);
Pt_urban_200_99 = Pthr - Gt - Gr + L_urban_200 + sqrt(2)*sigma_urban*erfinv(2*0.99 - 1);
Pt_urban_200_995 = Pthr - Gt - Gr + L_urban_200 + sqrt(2)*sigma_urban*erfinv(2*0.995 - 1);

% rural_30
Pt_rural_30_95 = Pthr - Gt - Gr + L_rural_30 + sqrt(2)*sigma_rural*erfinv(2*0.95 - 1);
Pt_rural_30_97 = Pthr - Gt - Gr + L_rural_30 + sqrt(2)*sigma_rural*erfinv(2*0.97 - 1);
Pt_rural_30_99 = Pthr - Gt - Gr + L_rural_30 + sqrt(2)*sigma_rural*erfinv(2*0.99 - 1);
Pt_rural_30_995 = Pthr - Gt - Gr + L_rural_30 + sqrt(2)*sigma_rural*erfinv(2*0.995 - 1);

% rural_200
Pt_rural_200_95 = Pthr - Gt - Gr + L_rural_200 + sqrt(2)*sigma_rural*erfinv(2*0.95 - 1);
Pt_rural_200_97 = Pthr - Gt - Gr + L_rural_200 + sqrt(2)*sigma_rural*erfinv(2*0.97 - 1);
Pt_rural_200_99 = Pthr - Gt - Gr + L_rural_200 + sqrt(2)*sigma_rural*erfinv(2*0.99 - 1);
Pt_rural_200_995 = Pthr - Gt - Gr + L_rural_200 + sqrt(2)*sigma_rural*erfinv(2*0.995 - 1);

figure(4);
plot(d,Pt_urban_30_95,"red");
hold on
plot(d,Pt_urban_30_97,"green");
hold on
plot(d,Pt_urban_30_99,"blue");
hold on
plot(d,Pt_urban_30_995,"magenta");
hold off
grid on
legend("Probability Coverage 95%","Probability Coverage 97%","Probability Coverage 99%","Probability Coverage 99.5%")
xlabel("Maximum Coverage Distance (km)");
ylabel("Transmit Power (dBm)")
title("Urban-H_BS=30m");

figure(5);
plot(d,Pt_urban_200_95,"red");
hold on
plot(d,Pt_urban_200_97,"green");
hold on
plot(d,Pt_urban_200_99,"blue");
hold on
plot(d,Pt_urban_200_995,"magenta");
hold off
grid on
legend("Probability Coverage 95%","Probability Coverage 97%","Probability Coverage 99%","Probability Coverage 99.5%")
xlabel("Maximum Coverage Distance (km)");
ylabel("Transmit Power (dBm)")
title("Urban-HBS=200m");

figure(6);
plot(d,Pt_rural_30_95,"red");
hold on
plot(d,Pt_rural_30_97,"green");
hold on
plot(d,Pt_rural_30_99,"blue");
hold on
plot(d,Pt_rural_30_995,"magenta");
hold off
grid on
legend("Probability Coverage 95%","Probability Coverage 97%","Probability Coverage 99%","Probability Coverage 99.5%")
xlabel("Maximum Coverage Distance (km)");
ylabel("Transmit Power (dBm)")
title("Rural-HBS=30m");

figure(7);
plot(d,Pt_rural_200_95,"red");
hold on
plot(d,Pt_rural_200_97,"green");
hold on
plot(d,Pt_rural_200_99,"blue");
hold on
plot(d,Pt_rural_200_995,"magenta");
hold off
grid on
legend("Probability Coverage 95%","Probability Coverage 97%","Probability Coverage 99%","Probability Coverage 99.5%")
xlabel("Maximum Coverage Distance (km)");
ylabel("Transmit Power (dBm)")
title("Rural-HBS=200m");
