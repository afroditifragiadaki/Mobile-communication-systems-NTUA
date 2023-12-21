
% Parameters
am=03118042;
am=mod(am,10);

mean_SNR = am; % Mean signal-to-noise ratio (dB)
N = -100 - am; % Background noise level (dBmW)

L = [1 , 2 , 3 , 4 , 5 , 6]; % Number of branches
samples = 1000000; % Number of Monte Carlo trials

%%%%%% Erwtima (a)
% Theoretical curve

ratio_SNR_dB = -40:0.01:20; % (dB)
ratio_SNR = 10.^(ratio_SNR_dB/10);

Pout_theory = zeros(length(L), length(ratio_SNR));
for ii = 1:length(L)
    Pout_theory(ii, :) = (1 - exp(-ratio_SNR)).^L(ii);
end

b = 10^(mean_SNR/20)*sqrt(2/pi);
gamma_branches = cell(length(L), 1);
gamma_SC = cell(length(L), 1);

for i = 1:length(L)
    gamma_branches{i} = raylrnd(b, samples, 1);
    gamma_SC{i} = max(cat(2, gamma_branches{1:i}), [], 2);
end

gamma_thr_dB = (ratio_SNR_dB + mean_SNR);
gamma_thr = 10.^(gamma_thr_dB/20);

Pout_SC_simulation = zeros(length(L), length(gamma_thr));

for i = 1:length(L)
    for j = 1:length(gamma_thr)
        Pout_SC_simulation(i, j) = sum(gamma_SC{i} < gamma_thr(j), 'all')/samples;
    end
end

% Plot results

figure(1);
L_values = [1 2 4 5];
L_colors = {'red', 'yellow', 'blue', 'green'};
for i = 1:length(L_values)
    semilogy(ratio_SNR_dB, Pout_theory(i, :), 'Color', L_colors{i});
    hold on;
end
for i = 1:length(L_values)
    semilogy(ratio_SNR_dB, Pout_SC_simulation(i, :), 'o', 'Color', L_colors{i});
    hold on;
end
hold off;
grid on;
axis([-20 20 10^(-5) 1]);
legend("Theoretical L=1", "Theoretical L=2", "Theoretical L=4","Theoretical L=5",...
    "Simulation L=1", "Simulation L=2", "Simulation L=4","Simulation L=5");
xlabel("SNR threshold /SNR_mean (dB)");
ylabel("Probability Outage");
axis([ -40 20 10^(-5) 1]);


%%%%Erwtima (b)

figure(2);
L_values = [1 3 4];
L_colors = {'red', 'yellow', 'blue'};
for i = 1:length(L_values)
    semilogy(gamma_thr_dB, Pout_theory(i, :), 'Color', L_colors{i});
    hold on;
end
for i = 1:length(L_values)
    semilogy(gamma_thr_dB, Pout_SC_simulation(i, :), 'o', 'Color', L_colors{i});
    hold on;
end
hold off;
grid on;
axis([-40 20 10^(-5) 1]);
legend("Theoretical L=1", "Theoretical L=3", "Theoretical L=4",...
    "Simulation L=1", "Simulation L=3", "Simulation L=4");
xlabel("SNR threshold (dB)");
ylabel("Probability Outage");


%%%%Erwtima (c)

%MRC
gamma_MRC = cell(length(L), 1);
gamma_MRC{1}=gamma_branches{1};
for i = 2:length(L)
     gamma_MRC{i}=gamma_MRC{i-1}+gamma_branches{i};
end
  
Pout_MRC_simulation = zeros(length(L), length(gamma_thr));
for i = 1:length(L)
    for j = 1:length(gamma_thr)
        Pout_MRC_simulation(i,j) = sum(gamma_MRC{i} < gamma_thr(j),'all')/samples;
    end
end

%EGC
gamma_EGC = cell(length(L), 1);
for i = 1:length(L)
    gamma_EGC{i} = sum(sqrt(gamma_branches{i}), 2).^2./L(i);
end
Pout_EGC_simulation = zeros(length(L), length(gamma_thr));
for i = 1:length(L)
    for j = 1:length(gamma_thr)
        Pout_EGC_simulation(i,j) = sum(gamma_EGC{i} < gamma_thr(j),'all')/samples;
    end
end

figure(3);
L_values = [1 2 4 6];
L_colors = {'red', 'yellow', 'blue', 'magenta'};

for i = 1:length(L_values)
    semilogy(gamma_thr_dB, Pout_SC_simulation(i, :), '-', 'Color', L_colors{i});
    hold on;
end
hold off;
grid on;
axis([-40 20 10^(-5) 1]);
legend("Simulation L=1", "Simulation L=2", "Simulation L=4","Simulation L=6");
xlabel("SNR threshold (dB)");
ylabel("Probability Outage");
title('Selection Combining SC')

figure(4);
for i = 1:length(L_values)
    semilogy(gamma_thr_dB, Pout_MRC_simulation(i, :), '-', 'Color', L_colors{i});
    hold on;
end
hold off;
grid on;
axis([-40 20 10^(-5) 1]);
legend("Simulation L=1", "Simulation L=2", "Simulation L=4","Simulation L=6");
xlabel("SNR threshold (dB)");
ylabel("Probability Outage");
title('Maximal Ratio Combining (MRC)')

figure(5);
for i = 1:length(L_values)
    semilogy(gamma_thr_dB, Pout_EGC_simulation(i, :), '-', 'Color', L_colors{i});
    hold on;
end
hold off;
grid on;
axis([-40 20 10^(-5) 1]);
legend("Simulation L=1", "Simulation L=2", "Simulation L=4","Simulation L=6");
xlabel("SNR threshold (dB)");
ylabel("Probability Outage");
title('Equal Gain Combining (EGC)')

figure(6);
semilogy(gamma_thr_dB,Pout_SC_simulation(1,:),"red");
hold on
semilogy(gamma_thr_dB,Pout_SC_simulation(3,:),"-O");
hold on
semilogy(gamma_thr_dB,Pout_MRC_simulation(3,:),"-+");
hold on
semilogy(gamma_thr_dB,Pout_EGC_simulation(3,:),"-x");
hold off
grid on
axis([-40 20 10^(-5) 1]);
legend("Simulation L=1","Simulation L=3 SC","Simulation L=3 MRC","Simulation L=3 ERC");
xlabel("SNR threshold (dB)");
ylabel("Probability Outage");

%%%%Erwtima (d)

L = 5;
L_vec = 2:L;
threshold_01_SC = zeros(L-1,1);
threshold_0001_SC = zeros(L-1,1);
gain_SC_01 = zeros(1,length(L_vec));
gain_SC_0001 = zeros(1,length(L_vec));

gain_MRC_01 = zeros(1,length(L_vec));
gain_MRC_0001 = zeros(1,length(L_vec));

gain_EGC_01 = zeros(1,length(L_vec));
gain_EGC_0001 = zeros(1,length(L_vec));

threshold_01_MRC = zeros(L-1,1);
threshold_0001_MRC = zeros(L-1,1);

threshold_01_EGC = zeros(L-1,1);
threshold_0001_EGC = zeros(L-1,1);

for ii = 1:length(L_vec)
    jj = L_vec(ii);
    [threshold_01_SC(ii),~] = ksdensity(gamma_SC{jj},10^(-1),'Function','icdf');
    [threshold_0001_SC(ii),~] = ksdensity(gamma_SC{jj},10^(-3),'Function','icdf');
    [threshold_01_MRC(ii),~] = ksdensity(gamma_MRC{jj},10^(-1),'Function','icdf');
    [threshold_0001_MRC(ii),~] = ksdensity(gamma_MRC{jj},10^(-3),'Function','icdf');
    [threshold_01_EGC(ii),~] = ksdensity(gamma_EGC{jj},10^(-1),'Function','icdf');
    [threshold_0001_EGC(ii),~] = ksdensity(gamma_EGC{jj},10^(-3),'Function','icdf');

    gain_SC_01(ii) = db(threshold_01_SC(ii),'power') - db(threshold_01_SC(1),'power');
    gain_SC_0001(ii) = db(threshold_0001_SC(ii),'power') - db(threshold_0001_SC(1),'power');
    
    gain_MRC_01(ii) = db(threshold_01_MRC(ii),'power') - db(threshold_01_MRC(1),'power');
    gain_MRC_0001(ii) = db(threshold_0001_MRC(ii),'power') - db(threshold_0001_MRC(1),'power');

    gain_EGC_01(ii) = db(threshold_01_EGC(ii),'power') - db(threshold_01_EGC(1),'power');
    gain_EGC_0001(ii) = db(threshold_0001_EGC(ii),'power') - db(threshold_0001_EGC(1),'power');

end
L = [2,3,4,5];
figure(7);
subplot(3,2,1);
stem(L,gain_SC_01); 
grid on; xlabel("Number of branches"); ylabel("Gain (dB)"); title(" SC - 0.1");
subplot(3,2,2);
stem(L,gain_SC_0001);
grid on; xlabel("Number of branches"); ylabel("Gain (dB)"); title("SC - 0.001");
subplot(3,2,3);
stem(L,gain_MRC_01);
grid on; xlabel("Number of branches"); ylabel("Gain (dB)"); title("MRC - 0.1");
subplot(3,2,4);
stem(L,gain_MRC_0001);
grid on; xlabel("Number of branches"); ylabel("Gain (dB)"); title(" MRC - 0.001");
subplot(3,2,5);
stem(L,gain_EGC_01);
grid on; xlabel("Number of branches"); ylabel("Gain (dB)"); title("EGC - 0.1");
subplot(3,2,6);
stem(L,gain_EGC_0001);
grid on; xlabel("Number of branches"); ylabel("Gain (dB)"); title("EGC - 0.001");