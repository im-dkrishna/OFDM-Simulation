MBPSK = 2;
MQPSK = 4;
MQAM = 16;
kBPSK = log2(MBPSK);
kQPSK = log2(MQPSK);
kQAM = log2(MQAM);

nSC = 4096;
nFFT = nSC;
nUsefulSC = 3072;
fSCspacing = 30;
nSymb = 14;

cpLen0 = 352;
cpLen1 = 288;
nPreamble = 64;
nSlots = 10;
nGB = (nSC - nUsefulSC)/2;
GB = zeros(nGB, nSymb);
endIndex = nUsefulSC + nGB;
zRow0 = zeros(cpLen0, 1);
zRow1 = zeros(cpLen1, 1);
pilotSymb1 = 2;
pilotSymb2 = 9;

chCoeff = [-0.014194636478986 0.022561094299837 -0.034046661636224 0.050235872796492 -0.074564210902286 0.116183398205243 -0.209922908839027 0.708391695270929 0.555511290188523 -0.192634871535184 0.109629537392724 -0.071008310551437 0.047959330035209 -0.032461758969592 0.021413701190345 -0.013355503650246];

maxBitErrors = 1e5;
maxNumBits = 1e6;

nDataCarriers = nUsefulSC*nSymb - nUsefulSC;
RGSize = [nUsefulSC nSymb];
bitSize_qpsk = [kQPSK*nDataCarriers 1];
nQpskBits = kQPSK*nDataCarriers;
frameSize_qpsk = [nQpskBits*nSlots 1];
bitSize_bpsk = [kBPSK*nDataCarriers 1];
nBpskBits = kBPSK*nDataCarriers;
frameSize_bpsk = [nBpskBits*nSlots 1];
bitSize_qam = [kQAM*nDataCarriers 1];
nQamBits = kQAM*nDataCarriers;
frameSize_qam = [nQamBits*nSlots 1];

EbNoVec = (0:12);
snrVec_qpsk = EbNoVec + 10*log10(kQPSK) + 10*log10(nDataCarriers/(nUsefulSC*nSymb));
berVec_qpsk = zeros(length(EbNoVec),3);
errorStats_qpsk = zeros(1,3);
snrVec_bpsk = EbNoVec + 10*log10(kBPSK) + 10*log10(nDataCarriers/(nUsefulSC*nSymb));
berVec_bpsk = zeros(length(EbNoVec),3);
errorStats_bpsk = zeros(1,3);
snrVec_qam = EbNoVec + 10*log10(kQAM) + 10*log10(nDataCarriers/(nUsefulSC*nSymb));

% QPSK Modulation
bpskMod = comm.BPSKModulator;
bpskDemod = comm.BPSKDemodulator;

% QPSK Modulation
qpskMod = comm.QPSKModulator('BitInput', true);
qpskDemod = comm.QPSKDemodulator('BitOutput',true);

% Channel
channel = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');
errorRate = comm.ErrorRate('ResetInputPort', true);

pilotValue = pilot1 + pilot2;

for m = 1:length(EbNoVec)
    snr = snrVec_qpsk(m);
    %snr = snrVec_bpsk(m);
    %snr = snrVec_qam(m);
    txFrame = [];
    txRG = [];
    rxRG = [];
    dataInQpskBits = randi([0,1], frameSize_qpsk);
    %dataInQpskBits = randi([0,1], frameSize_bpsk);
    %dataInQpskBits = randi([0,1], frameSize_qam);
    dataOutQpskBits = [];
    for j = 1:nSlots
        dataIn_qpsk = dataInQpskBits(1+(j-1)*nQpskBits:j*nQpskBits);        % Generate binary data
        %dataIn_qpsk = dataInQpskBits(1+(j-1)*nBpskBits:j*nBpskBits);
        %dataIn_qpsk = dataInQpskBits(1+(j-1)*nQamBits:j*nQamBits);
        qpskMap = qpskMod(dataIn_qpsk);                                     % Apply QPSK modulation
        %qpskMap = bpskMod(dataIn_qpsk);
        %qpskMap = qammod(dataIn_qpsk,MQAM,'InputType','bit','UnitAveragePower',true);
        
        qpskTx = zeros(RGSize);                                             % Bit stream to Resource Grid transformation
        qpskTx(:,1:2) = reshape(qpskMap(1:2*nUsefulSC),nUsefulSC,2);
        qpskTx(2:2:end,3) = qpskMap(2*nUsefulSC+1:2*nUsefulSC+nUsefulSC/2);
        qpskTx(:,4:9) = reshape(qpskMap(2*nUsefulSC+nUsefulSC/2+1:8*nUsefulSC+nUsefulSC/2),nUsefulSC, 6);
        qpskTx(2:2:end,10) = qpskMap(8*nUsefulSC+nUsefulSC/2+1:8*nUsefulSC+nUsefulSC/2+nUsefulSC/2);
        qpskTx(:,11:end) = reshape(qpskMap(8*nUsefulSC+nUsefulSC/2+nUsefulSC/2+1:end),nUsefulSC,4);
        qpskTx(1:2:end,pilotSymb1+1) = pilot1(1:2:end);                     % Adding pilot values to
        qpskTx(1:2:end,pilotSymb2+1) = pilot2(1:2:end);                     % symbol 2 and 9; alternative subcarriers
        
        qpskTxSig = [GB; qpskTx; GB];                                       % Padding with Guard Band on both sides
        txRG = [txRG qpskTxSig];
    end
    
    txRGSig = ifft(txRG);                                                   % IFFT
    
    for k = 1:nSlots
        txSig = txRGSig(:,1+(k-1)*nSymb:k*nSymb);
        cp1 = txSig(endIndex-cpLen0+1:endIndex,1);                          % Adding Cyclic Prefix to symbol 0
        txSigTemp = reshape(txSig, [], 1);
        txFrame = [txFrame; cp1; txSig(:,1)];
        for j = 2:14                                                        % Adding CP to rest of the symbols
            cp = txSig(endIndex-cpLen1+1:endIndex,j);
            txFrame = [txFrame; cp; txSig(:,j)];
        end
    end
    
    ofdmTxSig = reshape(txFrame, [], 1);                                    % Parallel to series
    powerDB = 10*log10(var(ofdmTxSig));                                     % Calculate Tx signal power
    noiseVar = 10.^(0.1*(powerDB-snr));                                     % Calculate the noise variance
    rxSig = channel(ofdmTxSig, noiseVar);                                   % Pass the signal through a noisy channel

    rxFrame = rxSig;
    index = 0;
    for j = 1:nSlots                                                        % Removing CP from rest of the symbols
        rxRG = [rxRG rxFrame(index+cpLen0+1:index+cpLen0+nSC)];
        index = index + cpLen0 + nSC;
        for i = 2:14
            rxRG = [rxRG rxFrame(index+cpLen1+1:index+cpLen1+nSC)];
            index = index + cpLen1 + nSC;
        end
    end
    
    qpskRx = fft(rxRG);                                                     % FFT
    
    qpskRx(1:nGB,:) = [];                                                   % Remove guard band
    qpskRx(nUsefulSC+1:end,:) = [];
    
    qpskRx = reshape(qpskRx, [], 1);
    for k = 1:nSlots
        qpskRxTemp = qpskRx(1+(k-1)*nUsefulSC*nSymb:k*nUsefulSC*nSymb);
        qpskRxTemp(2*nUsefulSC+1:2:3*nUsefulSC) = [];
        qpskRxTemp(9*nUsefulSC-nUsefulSC/2+1:2:9*nUsefulSC+nUsefulSC/2) = [];
        
        dataOut_qpsk = qpskDemod(qpskRxTemp);                               % Apply QPSK demodulation
        %dataOut_qpsk = bpskDemod(qpskRxTemp);
        %dataOut_qpsk = qamdemod(qpskRxTemp, MQAM,'OutputType','bit','UnitAveragePower',true);
        dataOutQpskBits = [dataOutQpskBits; dataOut_qpsk];
    end
    
    errorStats_qpsk = errorRate(dataInQpskBits,dataOutQpskBits,0);          % Collect error statistics 
    berVec_qpsk(m,:) = errorStats_qpsk;                                     % Save BER data
    errorStats_qpsk = errorRate(dataInQpskBits,dataOutQpskBits,1);          % Reset the error rate calculator
    
end

berqpsk = berawgn(EbNoVec,'psk',MQPSK,'nondiff');
berbpsk = berawgn(EbNoVec,'psk',MBPSK,'nondiff');
berqam = berawgn(EbNoVec,'qam',MQAM,'nondiff');
figure
semilogy(EbNoVec,berVec_qpsk(:,1),'*')
hold on
semilogy(EbNoVec,berbpsk)
semilogy(EbNoVec,berqpsk)
semilogy(EbNoVec,berqam)
legend('Simulation','BPSK','QPSK','QAM','Location','Best')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
grid on
hold off