%rng('default')
clear all
clc
%% channel generation
Hb = [0.32 0.66;1.24 0.58]; 

H21 = [-0.58 -1.15;
    -0.37 -1.07]; % channel betweel Alice to 1st Eve

H22 = [0.17 0.73;
    -0.07 -0.54]; % channel betweel Alice to 2nd Eve

He = [H21;H22]; % aggregated channel matrix as two Eves collaborating


H31 = [-0.23 -0.16;-0.05 -0.71]; % channel from Alice to 1st PR

H32 = [-0.47 0.09]; % channel from Alice to 2bnd PR


W31 = H31'*H31; % to form interfernece contraint
W32 = H32'*H32; % to form interfernece contraint

maxIter = 10;
%% transmit power parameters
SNRdB = [10 20];
IPC = cell(2);


% set up interference power constraints
Pi1 = 10^(5/10); % iterference limit for 1st PR (5dB)
Pi2 = 10^(5/10); % iterference limit for 1st PR (5dB)
IPC{1,1} = W31; 
IPC{1,2} = Pi1;
IPC{2,1} = W32; 
IPC{2,2} = Pi2;


%% run iterative algorithms of comparison
q = 5;
nA = size(Hb,2);
figure 
hold on

for iSNR=1:length(SNRdB)
    %
    P0 = 10^(SNRdB(iSNR)/10);
    PAPC = ones(nA,1)*P0/nA*1.2;
    S0 = zeros(nA)*P0/2;

    [~,objseqADCA] = Algorithm1ADCA_SPC_PAPC_IPC(Hb,He,P0,PAPC,IPC,S0,maxIter,q);
    plot(objseqADCA,'r')
    [~,~,objseqPBRA] = Algorithm3PBRA_SPC_PAPC_IPC(Hb,He,P0,PAPC,IPC,maxIter);
    plot(objseqPBRA,'b')
end
legend('Algorithm 1','Algorithm 3')
legend('Location','southeast')
xlabel('Iteration index')
saveas(gcf,'../results/Fig1b.png')