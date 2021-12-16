clear all
close all

fs=75000;
ts=1/fs;
t=0:ts:2;
kf=40*pi;
A=1;
wc=2000*pi;
N=100000;

%Señal de mensaje
figure(1)
m=(t).*(t>=0 & t<=1) + (-t + 2).*(t>1 & t<=1.9)  + 0.1*(t>1.9);
subplot(311)
plot(t,m)
axis([0 2 -1.5 1.5])
title('Señal m(t)')

%Integracion de m(t)
Int_m=cumsum(m)*ts;
figure(1)
subplot(312)
plot(t,Int_m)
title('Integral de m(t)')

%Modulacion de FM e 
yfm=A*cos(wc*t + kf*Int_m);
figure(1)
subplot(313)
plot(t,yfm)
title('Señal modulada FM')


%Espectro de la señal m(t)
M=fftshift(fft(m,N))*ts;
w=linspace(-fs/2, fs/2,N)*2*pi;
figure(2)
plot(w/(2*pi),abs(M))
axis([-100 100 0 1.1*max(abs(M))])
title('Espectro de magnitud de la señal moduladora m(t)')

%Ruido
long_m=length(m);
Py=var(yfm);
SNR=10;
snr=10^(SNR/10);
Pn=Py/snr;
sigma_n=sqrt(Pn);
n=sigma_n*randn(1,long_m);

%Grafica de señal de ruido garussiano
figure(3);
plot (t,n);
title('Ruido gaussiano')
Pnsim=var(m);

%Transformada de ruido
R=fftshift(fft(n,N))*ts;
f=linspace(-fs/2, fs/2, N);
figure(4);
plot(f, abs(R));
%axis([-500 500 0 0.06])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal de ruido');
grid

%Señal de ruido con la señal modulada
yfm_r = yfm + n;

%Modulacion de FM e 
figure(5)
plot(t,yfm_r)
title('Señal modulada FM')

%Espectro de la señal modulada
YFM_r=fftshift(fft(yfm_r,N))*ts;
w=linspace(-fs/2, fs/2,N)*2*pi;
figure(6)
plot(w/(2*pi),abs(YFM_r))
axis([-100 100 0 1.1*max(abs(YFM_r))])
title('Espectro de magnitud de la señal moduladora m(t)')

%Filtrado pasabanda en la entrada del receptor
th = -0.2:ts:0.2;
W=200*pi;
wc=2000*pi;
N=100000;
w=linspace(-fs/2, fs/2,N)*2*pi;
h1=2*W/pi*sinc(W*th/pi).*cos(wc*th);
H1=fftshift(fft(h1,N))*ts;

%Grafica de filtro pasabanda
figure(7)
subplot(211)
plot(th,h1,'m')
axis([-0.03 0.03 -340 1.1*max(abs(h1))])
title('Filtro pasabanda')

%Espectro de Filtro pasabandas
subplot(212)
plot(w/(2*pi),abs(H1),'m')
axis([-1300 1300 0 1.1*max(abs(H1))])
title('Espectro del Filtro pasabanda')
grid

%Convolucion de la señal y el filtro pasabanda
r1 = conv(yfm_r,h1,'same')*ts;

%Extrayendo la fase desde la transformada de Hilbert
z=hilbert(r1);

%Se obtiene la fase (Parte imaginaria)
theta=angle(z)-wc*t;
theta1=unwrap(theta);

%Señal integral de m(t)
figure(8)
plot(t,kf*Int_m)
title('kf*Integral de m(t)')

%Recuperacion de la fase de la tranformda de Hilbert
figure(9)
plot(t, theta1)
title('Recuperacion de la fase de yfm (kf*Integral de m(t))')

%Derivada de la señal de mensaje para obtencion de fase
m_recphase=[0 diff(theta1)/ts]/kf;

%Comparacion de señales recuperadas
figure(10)
plot(t,m,'b')
hold on
plot(t, m_recphase,'r')
axis([0 2 -1.5 1.5])
title('Comparación de la señal m(t) y la señal recuperada (con filtro)')

