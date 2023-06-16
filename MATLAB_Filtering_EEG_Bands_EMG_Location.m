

EEG_data_C3_4 = csvread('C:\Users\Mahsa\Desktop\Salar\a_dtw\csveeg\csv eeg\B-eeg.csv');
EMG_data_Right_3 = csvread('C:\Users\Mahsa\Desktop\Salar\a_dtw\EMG_salar\B-emg.csv');
%Total Size 1.176.300 - each secound 56,014
eeg = EEG_data_C3_4(:,1); 
emg = EMG_data_Right_3(:,3); 
V = eeg(:,1);
fs = 500;

% E = 3; % bin edge to discretize and categorize and classify 25% 50% 75%
% eeg_disc00=discretize(eeg,E);
% emg_disc00=discretize(emg,E);

long_eeg = lowpass(eeg,250,fs);
fc=1;% cut off frequency
fn=250; %nyquivst frequency = sample frequency/2;  sample interval of 500 Hz
order = 5; %5th order filter, high pass
[b1 a1]=butter(order,(fc/fn),'high');
x_data=filtfilt(b1,a1,long_eeg);
V=x_data;
waveletFunction = 'db8';
                    [C,L] = wavedec(V,8,waveletFunction); % Denoise EOG from EEG                  
                    D1 = wrcoef('d',C,L,waveletFunction,1);
                    D2 = wrcoef('d',C,L,waveletFunction,2);
                    D3 = wrcoef('d',C,L,waveletFunction,3);
                    D4 = wrcoef('d',C,L,waveletFunction,4);
                    D5 = wrcoef('d',C,L,waveletFunction,5); %GAMMA
                    D6 = wrcoef('d',C,L,waveletFunction,6); %BETA
                    D7 = wrcoef('d',C,L,waveletFunction,7); %ALPHA
                    D8 = wrcoef('d',C,L,waveletFunction,8); %THETA
                    A8 = wrcoef('a',C,L,waveletFunction,8); %DELTA  
%PSD Bands:
Gamma = D5;
Beta = D6;
Alpha = D7;
Theta = D8;
Delta = A8;
%1. The EMG was also ﬁrst downsampled to 250 Hz (with antialiasing down-pass ﬁltering)
%Reference: https://www.mathworks.com/help/signal/ref/lowpass.html#d123e98532
long = lowpass(emg,250,fs);
% lowpass(emg,250,fs);

%2.  and it was then high-pass ﬁltered with a 17th order Butterworth ﬁlter of above 110 Hz
%Reference: https://www.mathworks.com/matlabcentral/answers/5026-high-pass-butterworth-filter
fc=110;% cut off frequency
fn=250; %nyquivst frequency = sample frequency/2;  sample interval of 500 Hz
order = 17; %17th order filter, high pass
[b14 a14]=butter(order,(fc/fn),'high');
xf14=filtfilt(b14,a14,long);
% xemg = fvtool(b14,a14);

%3.  To generate Figure 2C we computed the Hilbert envelope of the signal and used it 
%  to obtain decibels of power density referred to the mean power of the signal during the epoch.
%Reference: https://www.mathworks.com/help/signal/ug/envelope-extraction-using-the-analytic-signal.html
y = hilbert(xf14);
env = abs(y);
iemg_plan = env;
%4. Similarity https://www.mathworks.com/help/signal/ug/compare-the-frequency-content-of-two-signals.html
Fs = fs;
Gamma2=abs(Gamma);
Gamma_plan = normalize(Gamma2, 'range', [0 0.06]);

figure(1);
subplot(2,1,1);
plot(Gamma_plan,'k');
grid;
ylabel('Gamma');
title('Power Spectrum');

subplot(2,1,2);
plot(iemg_plan,'r');
grid;
ylabel('iEmg');
xlabel('Record Number (1.176.300 )');


Edge = 3; % bin edge to discretize and categorize and classify 0% 25% 50%
[eeg_disc,Eeeg] =discretize(Gamma_plan ,Edge);
[emg_disc,Eemg] =discretize(iemg_plan,Edge);
% disp(Eeeg)
% disp(Eemg)
len_eeg=length(eeg_disc(:,1));
s=10; %10 gripping priod = 21 second
interval = len_eeg/s;
eeg_mat=zeros(1,s);
emg_mat=zeros(1,s);
for i = 1:s
    if i==1
        avg_eeg=(sum(eeg_disc(1:interval,1)))/(interval);
        avg_emg=(sum(emg_disc(1:interval,1)))/(interval);
        disp(i+":  avg_eeg: "+(avg_eeg) +" avg_emg: " + (avg_emg));
        eeg_mat(1,i)=avg_eeg;
        emg_mat(1,i)=avg_emg;
    else
       
        avg_eeg=(sum(eeg_disc(((i-1)*interval):(i*interval),1)))/(interval);
        avg_emg=(sum(emg_disc(((i-1)*interval):(i*interval),1)))/(interval);
        disp(i+":  avg_eeg: "+(avg_eeg) +" avg_emg: " + (avg_emg));
        eeg_mat(1,i)=avg_eeg;
        emg_mat(1,i)=avg_emg;
    end
    
end
figure(2);

x = 1:s;
y1 = eeg_mat;
y2 = emg_mat;

bar(x,[y1;y2]); % Creating both bars at once, they are aware of one another and will not overlap
% Legend will show names for each color
legend('EEG','EMG');
