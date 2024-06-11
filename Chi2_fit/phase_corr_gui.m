function [] = phase_corr_gui(data, sw, zf, lb,p0, p1, basl, win, zp, fs)
%Function for interactive 0th and 1st order phase correction
%   phase_corr_gui(data, sw, zf, lb, p0, p1, basl, win, zp, fs)
%K. Aebischer, 25.06.20
%Input:
    %data:      complex Signal
    %sw:        Spectral width (Hz)
    %zf:        zero-filling
    %lb:        Line-broadening (Hz)
    %basl:      Integer code for baseline correction (0: none, 1: in FID, 2: in FID and spectrum)
    %win:       Integer code for apodization (0: none, else: cos^2)
    %zp:        
    %fs:        

data=data(:).';         %ensures data is a row vector
dw = 1/sw;              %dwell time (time-res. of FID, determined by spectral width)
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
    p0=0;
end
if nargin < 6
    p1 = 0;
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

%1st order phase correction due to time-shift
p1_init = -fs*180;

if zf == 0
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
    % cos^2 window
    apod = ones(1,l);
    apod(fs+1:l) = cos((0:(l-(fs+1)))/(l-(fs+1))*pi/2).^2;
end
data1 = data .* apod;
data1(1)=0.5*data1(1);

%Fourier tranform FID
    %center around 0 frq
data_ft = fftshift(fft(data1,zf));

%1st order phase correction
    %due to time shift specified in input
x = (((0:(zf-1))-zp+zf/2)/zf)-0.5;
data_ft_corr = data_ft .* exp(-1i*pi/180*2*x*(p1_init));
xax = x*sw;

%Offset correction for spectrum
if basl > 1
    %use edge of spectral width for baseline correction
    offset = mean(data_ft_corr(round(0.1*l):round(0.2*l)));
    data_ft_corr = data_ft_corr - offset;
end

%Create figure and plot spectrum
FigH = figure();
axes('XLim', [min(xax) max(xax)], 'units','normalized', ...
     'position',[0.15 0.15 0.7 0.6], 'NextPlot', 'add');
xlabel('\nu / Hz')
%Plot spectrum with phase correction values given as input
LineH = plot(xax, real(data_ft_corr * exp(-1i*pi/180*p0).*exp(-1i*pi/180*2*x*p1)), 'LineWidth', 1);
grid on

%Create Sliders for p0 and p1
TextP0 = uicontrol('style','text','HorizontalAlignment','Center',...
                   'units','normalized', 'position',[0.2 0.85 0.275 0.05]);
TextP0.String = "0th order : " + num2str(p0) + "°";

TextP1 = uicontrol('style','text','HorizontalAlignment','Center',...
                   'units','normalized', 'position',[0.525 0.85 0.275 0.05]);
TextP1.String = "1st order : " + num2str(p1) + "°";               
               
SliderP0 = uicontrol('style','slider', 'units', 'normalized', ...
                     'value', p0, 'position',[0.2 0.8 0.275 0.05],'min', -30, 'max', 360);
addlistener(SliderP0, 'Value', 'PostSet', @callbackfnP0);

SliderP1 = uicontrol('style','slider','units', 'normalized',...
                       'position',[0.525 0.8 0.275 0.05],'min', -180, 'max', 180);
addlistener(SliderP1, 'Value', 'PostSet', @callbackfnP1);

movegui(FigH, 'center')

%Callback functions for updating spectrum
function callbackfnP0(~, eventdata)
    num          = get(eventdata.AffectedObject, 'Value');
    p0 = num;
    LineH.YData = real(data_ft_corr * exp(-1i*pi/180*p0).*exp(-1i*pi/180*2*x*p1));
    TextP0.String = "0th order : " + num2str(num) + "°";
end

function callbackfnP1(~, eventdata)
    num          = get(eventdata.AffectedObject, 'Value');
    p1 = num;
    LineH.YData  = real(data_ft_corr * exp(-1i*pi/180*p0).*exp(-1i*pi/180*2*x*p1));
    TextP1.String = "1st order : " + num2str(num) + "°";
end

end
