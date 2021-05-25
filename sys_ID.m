%% Dataset 1

raw_data = readtable('INTERP_2_7-5_E4pin_20mm_v9_0.2Strain_velocityprof.csv');

Ro = table2array(raw_data(:,1));
Ri = table2array(raw_data(:,2));
Resistivity = table2array(raw_data(:,3));
Strain = table2array(raw_data(:,4));
Stress = table2array(raw_data(:,5));
t_lin = table2array(raw_data(:,6));

Fs = 8.76;

y = Stress;
ydft = fftshift(fft(y));

u = Strain;
udft = fftshift(fft(u));


G = fftshift(fft(y)./fft(u)); % transfer function
G_fil = lowpass(G,8,8.7)

df_G = Fs/length(u);
freq_G = -Fs/2:df_G:Fs/2-df_G;

figure
% plot(freq,abs(udft))
% plot(freq,abs(ydft))
plot(freq_G,abs(G))
figure
% plot(t_lin, Stress)
% figure
plot(t_lin, Strain)

%% Test data (Dataset 2)

raw_data = readtable('INTERP_2_7-5_E4pin_20mm_v10_0.3Strain.csv');

Ro_test = table2array(raw_data(:,1));
Ri_test = table2array(raw_data(:,2));
Resistivity = table2array(raw_data(:,3));
strain_test = table2array(raw_data(:,4));
stress_test = table2array(raw_data(:,5));
t_test = table2array(raw_data(:,6));

% max_strain = 0.3
% t_test = linspace(0,3000, 30000)
% strain_test = 0.5*(square(2*pi*(1/800)*t_test)+1)*max_strain
figure
plot(t_test,strain_test)

u = strain_test;
udft = fftshift(fft(u));

df = Fs/length(u);
freq = -Fs/2:df:Fs/2-df;

figure
plot(freq,abs(udft));

ydft = ifft(G .* udft);
ydft_fil = ifft(G_fil .* udft);

t_dft = linspace(0,max(t_test),length(ydft))

figure;hold on
plot(t_dft,ydft)
plot(t_test,stress_test)
figure;hold on
plot(t_dft,ydft_fil)
plot(t_test,stress_test)



%% Dataset 2

raw_data = readtable('INTERP_2_7-5_E4pin_20mm_v10_0.3Strain.csv');

Ro = table2array(raw_data(:,1));
Ri = table2array(raw_data(:,2));
Resistivity = table2array(raw_data(:,3));
Strain = table2array(raw_data(:,4));
Stress = table2array(raw_data(:,5));
t_lin = table2array(raw_data(:,6));

Fs = 8.76;

y = Stress;
ydft = fftshift(fft(y));

u = Strain;
udft = fftshift(fft(u));


G = fftshift(fft(y)./fft(u)); % transfer function?

df = Fs/length(u);
freq = -Fs/2:df:Fs/2-df;

figure
% plot(freq,abs(udft))
% plot(freq,abs(ydft))
plot(freq,abs(G))
figure
% plot(t_lin, Stress)
% figure
plot(t_lin, Strain)

% %% Dataset 3
% 
% raw_data = readtable('INTERP_2_7-5_E4pin_20mm_v17_0.3Strain.csv');
% 
% Ro = table2array(raw_data(:,1));
% Ri = table2array(raw_data(:,2));
% Resistivity = table2array(raw_data(:,3));
% Strain = table2array(raw_data(:,4));
% Stress = table2array(raw_data(:,5));
% t_lin = table2array(raw_data(:,6));
% 
% Fs = 8.76;
% 
% y = Stress;
% ydft = fftshift(fft(y));
% 
% u = Strain;
% udft = fftshift(fft(u));
% 
% 
% G = fftshift(fft(y)./fft(u)); % transfer function?
% 
% df = Fs/length(u);
% freq = -Fs/2:df:Fs/2-df;
% 
% figure
% % plot(freq,abs(udft))
% % plot(freq,abs(ydft))
% plot(freq,abs(G))
% figure
% % plot(t_lin, Stress)
% % figure
% plot(t_lin, Strain)
% 
% %% Dataset 4
% 
% raw_data = readtable('INTERP_2_7-5_E4pin_20mm_v18_0.3Strain.csv');
% 
% Ro = table2array(raw_data(:,1));
% Ri = table2array(raw_data(:,2));
% Resistivity = table2array(raw_data(:,3));
% Strain = table2array(raw_data(:,4));
% Stress = table2array(raw_data(:,5));
% t_lin = table2array(raw_data(:,6));
% 
% Fs = 8.76;
% 
% y = Stress;
% ydft = fftshift(fft(y));
% 
% u = Strain;
% udft = fftshift(fft(u));
% 
% 
% G = fftshift(fft(y)./fft(u)); % transfer function?
% 
% df = Fs/length(u);
% freq = -Fs/2:df:Fs/2-df;
% 
% figure
% % plot(freq,abs(udft))
% % plot(freq,abs(ydft))
% plot(freq,abs(G))
% figure
% % plot(t_lin, Stress)
% % figure
% plot(t_lin, Strain)
% 
% %% Dataset 5
% 
% raw_data = readtable('INTERP_2_7-5_E4pin_20mm_v16_0.3Strain.csv');
% 
% Ro = table2array(raw_data(:,1));
% Ri = table2array(raw_data(:,2));
% Resistivity = table2array(raw_data(:,3));
% Strain = table2array(raw_data(:,4));
% Stress = table2array(raw_data(:,5));
% t_lin = table2array(raw_data(:,6));
% 
% Fs = 8.76;
% 
% y = Stress;
% ydft = fftshift(fft(y));
% 
% u = Strain;
% udft = fftshift(fft(u));
% 
% 
% G = fftshift(fft(y)./fft(u)); % transfer function?
% 
% df = Fs/length(u);
% freq = -Fs/2:df:Fs/2-df;
% 
% figure
% % plot(freq,abs(udft))
% % plot(freq,abs(ydft))
% plot(freq,abs(G))
% figure
% % plot(t_lin, Stress)
% % figure
% plot(t_lin, Strain)

%% Compress FFT data

udft_new = compress_dft(freq, freq_G, udft);

function u_A_new = compress_dft(f_A, f_B, u_A)
    % Compressive freq binning from size A to B of DFT data with magnitudes u and frequencies f
    if length(f_A) < length(f_B)
        error('Error: f_A is smaller than f_B vector')
    end
    in_bin = 0; % how many frequencies from f_A are being binned into one of f_B's bins
    bin_sum = 0;
    bindex = 1;
    
    % lower and upper limits of new bin
    bin_lo = f_B(bindex);
    bin_hi = f_B(bindex+1);
    
    % init new compressed fft data
    u_A_new = zeros(1,length(f_B));
    
    for i = 1:length(f_A)
       bin_lo = f_B(bindex)
       bin_hi = f_B(bindex+1)
       f_A(i)
       if ((f_A(i) >= bin_lo) && (f_A(i) < bin_hi))
           bin_sum = bin_sum + u_A(i)
           in_bin = in_bin + 1
       elseif (f_A(i) > bin_hi)
           u_A_new(bindex) = bin_sum/in_bin;
           bin_sum = 0
           in_bin = 0
           bindex = bindex + 1
           bin_sum = bin_sum + u_A(i)
           in_bin = in_bin + 1
       elseif (f_A(i) < bin_lo)
           error('Error: data out of order or bin missed')
       else
           error('what the fluff')
       end
    end

end
