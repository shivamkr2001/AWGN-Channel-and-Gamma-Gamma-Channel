clc;
clear all;
close all;

N = 10^4; % Number of Samples
BER_w = []; % Array to store BER for each SNR for weak turbulence
BER_m = []; % Array to store BER for each SNR for moderate turbulence
BER_s = []; % Array to store BER for each SNR for strong turbulence
alpha_w = 11.6; % Alpha for weak turbulence
alpha_m = 4; % Alpha for moderate turbulence
alpha_s = 4.2; % Alpha for strong turbulence
beta_w = 10.1; % Beta for weak turbulence
beta_m = 1.9; % Beta for medium turbulence
beta_s = 1.4; % Beta for strong turbulence

for SNRdB = 0:1:15 % SNR in dB
     SNR = 10^(SNRdB/10); % SNR
     input = randi([0,1],1,N); % Randomly generating 0 and 1
     data = 2*input - 1; % Converting data to -1 and 1
     error_w = 0;
     error_m = 0;
     error_s = 0;
     A_w = gamrnd(alpha_w, 1/alpha_w, [1,N]);
     B_w = gamrnd(beta_w, 1/beta_w, [1,N]);
     H_w = A_w .* B_w;
     A_m = gamrnd(alpha_m, 1/alpha_m, [1,N]);
     B_m = gamrnd(beta_m, 1/beta_m, [1,N]);
     H_m = A_m .* B_m;
     A_s = gamrnd(alpha_s, 1/alpha_s, [1,N]);
     B_s = gamrnd(beta_s, 1/beta_s, [1,N]);
     H_s = A_s .* B_s;
     
     for i = 1:N
         n = sqrt(1/2)*randn(1,1); % Distribution of noise N(0,1/2)
         y_w = (H_w(i)*sqrt(SNR)*data(i)) + n; % Received Signal for weak turbulence
         y_m = (H_m(i)*sqrt(SNR)*data(i)) + n; % Received Signal for moderate turbulence
         y_s = (H_s(i)*sqrt(SNR)*data(i)) + n; % Received Signal for strong turbulence
         if(y_w > 0)
             received_w = 1;
         else
             received_w = -1;
         end
          if(y_m > 0)
             received_m = 1;
         else
             received_m = -1;
          end
          if(y_s > 0)
             received_s = 1;
         else
             received_s = -1;
         end
     error_w = error_w + nnz(received_w - data(i));
     error_m = error_m + nnz(received_m - data(i));
     error_s = error_s + nnz(received_s - data(i));
     end
     
     BER_w = [BER_w error_w/N];
     BER_m = [BER_m error_m/N];
     BER_s = [BER_s error_s/N];
end
  
SNRdB = 0:1:15;
SNR = 10.^(SNRdB/10); % SNR
A_w = [(2-alpha_w)/2 (1-alpha_w)/2 (2-beta_w)/2 (1-beta_w)/2];
A_m = [(2-alpha_m)/2 (1-alpha_m)/2 (2-beta_m)/2 (1-beta_m)/2];
A_s = [(2-alpha_s)/2 (1-alpha_s)/2 (2-beta_s)/2 (1-beta_s)/2];
B = [1];
C = [0 0.5];
D = [];
Z_w = 16*SNR/(alpha_w*beta_w)^2;
Z_m = 16*SNR/(alpha_m*beta_m)^2;
Z_s = 16*SNR/(alpha_s*beta_s)^2;
BER_tw = (2^(alpha_w+beta_w-3)/(pi^1.5 * gamma(alpha_w) * gamma(beta_w))).* meijerG(A_w,B,C,D,Z_w); % Theoretical BER for a Weak Turbulent Channel
BER_tm = (2^(alpha_m+beta_m-3)/(pi^1.5 * gamma(alpha_m) * gamma(beta_m))).* meijerG(A_m,B,C,D,Z_m); % Theoretical BER for a Medium Turbulent Channel
BER_ts = (2^(alpha_s+beta_s-3)/(pi^1.5 * gamma(alpha_s) * gamma(beta_s))).* meijerG(A_s,B,C,D,Z_s); % Theoretical BER for a Strong Turbulent Channel
semilogy(SNRdB,BER_w,'r');
hold on
semilogy(SNRdB,BER_m,'--');
hold on
semilogy(SNRdB,BER_s,':');
hold on
semilogy(SNRdB,BER_tw,'+');
hold on
semilogy(SNRdB,BER_tm,'o');
hold on
semilogy(SNRdB,BER_ts,'*');
legend('Weak','Moderate','Strong','Weak Theoretical','Moderate Theoretical','Strong Theoretical');
ylim([1e-4 1]);
xlabel('SNR in dB');
ylabel('Bit Error Rate');
