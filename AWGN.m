clc;
clear all;
close all;

N = 10^6; % Number of Samples
BER = []; % Array to store BER for each SNR
for SNRdB = 0:1:10 % SNR in dB
     SNR = 10^(SNRdB/10); % SNR
     input = randi([0,1],1,N); % Randomly generating 0 and 1
     data = 2*input - 1; % Converting data to -1 and 1
     error = 0;
     
     for i = 1:N
         n = sqrt(1/2)*randn(1,1); % Distribution of noise N(0,1/2)
         y = (sqrt(SNR)*data(i)) + n; % Received Signal
         if(y > 0)
             received = 1;
         else
             received = -1;
         end
     error = error + nnz(received - data(i));
     
     end
     
     BER = [BER error/N];
     
end
  
SNRdB = 0:1:10;
SNR = 10.^(SNRdB/10); % SNR
BER_theo = qfunc(sqrt(2*SNR)); % Theoretical BER
semilogy(SNRdB,BER);
hold on;
semilogy(SNRdB,BER_theo,'+');
legend('Experimental','Theoretical');
ylabel('Bit Error Rate');
xlabel('SNR in dB');