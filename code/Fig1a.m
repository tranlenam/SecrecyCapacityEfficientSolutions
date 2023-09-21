clear
clc
%% system setup

% the channels considered in Fig. 1
Hb = [-0.397388134804813 + 0.564069068925790i	-0.0939337320157974+0.253217235750467i;
-0.0216281544185800+0.805107628931355i	-0.673432238998531+0.260482934998363i;
-1.19032291480353-0.393916696585856i	-0.972781640103893-0.446827692307842i;
0.201667162510618-0.689726745440122i	-0.945047869851278-0.730641339097362i];

He = [-0.201455080897275+0.312725917599147i -0.617804134199973-1.04797479362569i;
-0.0558694759682173-0.300019080432181i	-0.385799249500664-0.281697157494074i;
0.693547640268284+0.0558760342293484i	-0.506416527519376-0.144333099022521i];
maxIter = 100;
%% transmit power parameters
SNRdB = [5 15];

%% run iterative algorithms of comparison
q = 5;
nA = size(Hb,2);
figure 
hold on
for iSNR=1:length(SNRdB)
    %
    P0 = 10^(SNRdB(iSNR)/10);
    PAPC = ones(nA,1)*P0/nA*1.2;
    S0 = eye(nA)*P0/(nA*2);

    [~,objseqADCA] = Algorithm1ADCA_SPC_PAPC(Hb,He,P0,PAPC,S0,maxIter,q);
    plot(objseqADCA,'r')
    [~,~,objseqPBRA] = Algorithm3PBRA_SPC_PAPC(Hb,He,P0,PAPC,maxIter);
    plot(objseqPBRA,'b')
end
legend('Algorithm 1','Algorithm 3')
legend('Location','southeast')
xlabel('Iteration index')
saveas(gcf,'../results/Fig1a.png')