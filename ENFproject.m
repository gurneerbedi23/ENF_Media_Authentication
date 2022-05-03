%% Examining the harmonics of the ENF

% Frequency = 60 Hz 
x = transpose(audioread('recording.wav'));
Fs = 44100;
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 60;

[y1, y2, y3] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[m,n] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
k1 = Frequency-1;
k2 = Frequency+1;
factor = (k2-k1)/n;
a = k1:factor:(k2-factor); %made length of vector for frequency divisions an element extra
b = 1:m;

surf(a,b,y1); title('ENF Harmonic at Frequency = 60 Hz'); ylabel('Block Number'); xlabel('Frequency (Hz)'); 
%% Frequency = 120 Hz 
x = transpose(audioread('recording.wav'));
Fs = 44100;
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120;

[y1, y2, y3] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[m,n] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
k1 = Frequency-1;
k2 = Frequency+1;
factor = (k2-k1)/n;
a = k1:factor:(k2-factor);
b = 1:m;

surf(a,b,y1); title('ENF Harmonic at Frequency = 120 Hz'); ylabel('Block Number'); xlabel('Frequency (Hz)');
%% Frequency = 180 Hz 
x = transpose(audioread('recording.wav'));
Fs = 44100;
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 180;

[y1, y2, y3] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[m,n] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
k1 = Frequency-1;
k2 = Frequency+1;
factor = (k2-k1)/n;
a = k1:factor:(k2-factor);
b = 1:m;

surf(a,b,y1); title('ENF Harmonic at Frequency = 180 Hz'); ylabel('Block Number'); xlabel('Frequency (Hz)'); 
%% Frequency = 240 Hz 
x = transpose(audioread('recording.wav'));
Fs = 44100;
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 240;

[y1, y2, y3] = enf(x, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[m,n] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
k1 = Frequency-1;
k2 = Frequency+1;
factor = (k2-k1)/n;
a = k1:factor:(k2-factor);
b = 1:m;

surf(a,b,y1); title('ENF Harmonic at Frequency = 240 Hz'); ylabel('Block Number'); xlabel('Frequency (Hz)'); 
%% Comparison of reference and test sound files without pre-processing

ref = transpose(audioread('ground truth.wav'));
test = transpose(audioread('recording.wav'));
Fs = 44100;
BlockSize = Fs*16;
ZeroPad = 0;
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120;

[y1, y2, y3] = enf(ref, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[x1, x2, x3] = enf(test, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);  
k1 = Frequency-1;
k2 = Frequency+1;

[refrow,refcol] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
factor = (k2-k1)/refcol;
a = k1:factor:(k2-factor);
b = 1:refrow;

[testrow,testcol] = size(x1); 
factor2 = (k2-k1)/testcol;
c = k1:factor2:(k2-factor2);
d = 1:testrow;

surf(a,b,y1); title('ENF of Reference Signal "Ground Truth"'); ylabel('Block Number'); xlabel('Frequency (Hz)'); figure;
surf(c,d,x1); title('ENF of Test Signal "Recording"'); ylabel('Block Number'); xlabel('Frequency (Hz)'); figure;

% Finding the mean and subtracting it from the signals
meanGround = mean(y3);
GroundTruth = y3 - meanGround;LPFref

meanRecord = mean(x3);
Recording = x3 - meanRecord;

[r,lags] = xcorr(GroundTruth, Recording);
plot(lags,r); title('Graph of Lags when Correlations are Computed'); xlabel('Lags'); ylabel('Cross-Correlation'); figure;
plot(y3); hold on; plot(x3); title('Weighted Magnitude Comparision Before Allignment'); xlabel('Block Number'); ylabel('Magnitude'); figure;
plot(y2); hold on; plot(x2); title('Maximum Magnitude Comparision Before Allignment'); xlabel('Block Number'); ylabel('Magnitude'); figure;

%Zero padding the shorter signal to match the plots, from xcorr found lag
%of 2 windows

newx2 = padarray(x2,2,meanRecord);
newx3 = padarray(x3,2,meanRecord); 
plot(y3); hold on; plot(newx3); title('Weighted Magnitude Comparision After Allignment'); xlabel('Block Number'); ylabel('Magnitude'); figure;
plot(y2); hold on; plot(newx2); title('Maximum Magnitude Comparision After Allignment'); xlabel('Block Number'); ylabel('Magnitude'); 
%% Comparison of reference and test sound files with pre-processing

load('SOS.mat');
load('G.mat');
ref = transpose(audioread('ground truth.wav'));
test = transpose(audioread('recording.wav'));
LPFref = filtfilt(SOS,G,ref);
LPFtest = filtfilt(SOS,G,test); 
LPFref = downsample(LPFref,100);
LPFtest = downsample(LPFtest,100);

Fs = 441;
BlockSize = Fs*16;
ZeroPad = (2^14) - (Fs*16);
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120;

[y1, y2, y3] = enf(LPFref, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[x1, x2, x3] = enf(LPFtest, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);  
k1 = Frequency-1;
k2 = Frequency+1;

[refrow,refcol] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
factor = (k2-k1)/refcol;
a = k1:factor:(k2-factor);
b = 1:refrow;

[testrow,testcol] = size(x1); 
factor2 = (k2-k1)/testcol;
c = k1:factor2:(k2-factor2);
d = 1:testrow;
surf(a,b,y1); title('ENF of Reference Signal "Ground Truth" After Pre-Processing'); ylabel('Block Number'); xlabel('Frequency (Hz)'); figure;
surf(c,d,x1); title('ENF of Test Signal "Recording" After Pre-Processing'); ylabel('Block Number'); xlabel('Frequency (Hz)'); figure; 

% Finding the mean and subtracting it from the signals
meanGround = mean(y3);
GroundTruth = y3 - meanGround;

meanRecord = mean(x3);
Recording = x3 - meanRecord;

[r,lags] = xcorr(GroundTruth, Recording);
plot(lags,r); title('Graph of Lag when Correlations are Computed'); xlabel('Lags'); ylabel('Cross-Correlation'); figure;
plot(y3); hold on; plot(x3); title('Weighted Magnitude Comparision Before Allignment'); xlabel('Block Number'); ylabel('Magnitude'); figure;

%Padding the array with the delay found earlier
newx3 = padarray(x3,2,meanRecord); 
plot(y3); hold on; plot(newx3); title('Weighted Magnitude Comparision After Allignment'); xlabel('Block Number'); ylabel('Magnitude');
%% Comparison of second set of sound files with pre-processing

load('SOS.mat');
load('G.mat');
ref = transpose(audioread('ground truth 2.wav'));
test = transpose(audioread('recording 2.wav'));
LPFref = filtfilt(SOS,G,ref);
LPFtest = filtfilt(SOS,G,test); 
LPFref = downsample(LPFref,100);
LPFtest = downsample(LPFtest,100);

Fs = 441;
BlockSize = Fs*16;
ZeroPad = (2^14) - (Fs*16);
Overlap = 0.5;
Window = hanning(BlockSize, 'periodic');
Frequency = 120;

[y1, y2, y3] = enf(LPFref, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);
[x1, x2, x3] = enf(LPFtest, Fs, BlockSize, ZeroPad, Overlap, Window, Frequency);  

[refrow,refcol] = size(y1); % [m,n] = size(y1), x = 1:n, y = 1:m 
factor = (k2-k1)/refcol;
a = k1:factor:(k2-factor);
b = 1:refrow;

[testrow,testcol] = size(x1); 
factor2 = (k2-k1)/testcol;
c = k1:factor2:(k2-factor2);
d = 1:testrow;
surf(a,b,y1); title('ENF of Reference Signal "Ground Truth 2" After Pre-Processing'); ylabel('Block Number'); xlabel('Frequency (Hz)'); figure;
surf(c,d,x1); title('ENF of Test Signal "Recording 2" After Pre-Processing'); ylabel('Block Number'); xlabel('Frequency (Hz)'); figure;

% Finding the mean and subtracting it from the signals
meanGround = mean(y3);
GroundTruth = y3 - meanGround;

meanRecord = mean(x3);
Recording = x3 - meanRecord;

[r,lags] = xcorr(GroundTruth, Recording);
plot(lags,r); title('Graph of Lag when Correlations are Computed'); xlabel('Lags'); ylabel('Cross-Correlation'); figure;
plot(y3); hold on; plot(x3); title('Weighted Magnitude Comparision Before Allignment'); xlabel('Block Number'); ylabel('Magnitude'); figure;

%Padding the array with the delay found earlier
newx3 = padarray(x3,2,meanRecord); 
plot(y3); hold on; plot(newx3); title('Weighted Magnitude Comparision After Allignment'); xlabel('Block Number'); ylabel('Magnitude');

