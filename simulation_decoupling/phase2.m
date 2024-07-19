function [data_ft, frq_ax] = proc_fid(data,sw,zf,lb,phase0, phase1, basl, win, zp, fs)
%K. Aebischer, 10.06.20
%based on phase1.m by M. Ernst
%Function for basic processing of FID signal including zero-filling,
%baseline and phase correction and apodization
%[data_ft, frq_ax] = proc_fid(data,sw,zf,lb,phase0, phase1, basl, win, zp, fs)
%
%Input
    %data:      array with FID datapoints
    %sw:        spectral width in Hz
    %zf:        Zero-filling
    %lb:        line-broadening (Hz) for exponential multiplication
    %phase0:    deg, 0 order phase correction
    %phase1:    deg, 1st order phase correction
    %basl:      option for baseline correction (0: none, 1: on FID, 2: FID and spectrum)
    %win:       apodization window (0: exponential, else: cos^2)
    %zp:        zero-point for first order phase correction (point index)
    %fs:        Shift FID by fs points (protection delay)
%Output:
    %data_ft:   spectrum of FID
    %frq_ax:    frequency axis of spectrum in Hz

data=data(:).';         %ensures data is a row vector
dw = 1/sw;              %dwell time (time-res. of FID)
l=length(data);

%Check input arguments
    %set default values if no argument given
if nargin < 3
    zf=l;
end
if nargin < 4
    lb=0;
end
if nargin < 5
    phase0=0;
end
if nargin < 6
    phase1=0;
end
if nargin < 7
    basl=0;
end
if nargin < 8
    win=0;
end
if nargin < 9
    zp=zf/2;
end
if nargin < 10
    fs=0;
end

phase1 = phase1-fs*180;
if zf==0
  zf = l;
end

%Offset correction FID
    %takes the last 20% of the FID and corrects by mean value
if basl > 0
  offset = mean(data(round(0.8*l):l));
  data(fs+1:l) = data(fs+1:l) - offset;
end

%Apodization
if win == 0
    %exponential window
    apod = ones(1,l);
    apod(fs+1:l) = exp(-lb*(0:(l-(fs+1)))*dw*pi);
else
    %cos^2 window
    apod = ones(1,l);
    apod(fs+1:l) = cos((0:(l-(fs+1)))/(l-(fs+1))*pi/2).^2;
end
data = data .* apod;
data(1)=0.5*data(1);

%Fourier tranform FID
data_ft = fftshift(fft(data,zf));

%0 order phase correction
data_ft = data_ft * exp(-1i*pi/180*phase0);

%1st order phase correction
x = (((0:(zf-1))-zp+zf/2)/zf)-0.5;
frq_ax = x*sw;
data_ft = data_ft .* exp(-1i*pi/180*2*x*phase1);

%Offset correction for spectrum
if basl > 1
    %use edge of spectral width for baseline correction
    offset = mean(data_ft(round(0.1*l):round(0.2*l)));
    data_ft = data_ft - offset;
end

end%end of function
