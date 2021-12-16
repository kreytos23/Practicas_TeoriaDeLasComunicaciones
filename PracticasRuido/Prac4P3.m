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
plot(t,m)
axis([0 2 -1.5 1.5])
title('Señal m(t)')

%Integracion de m(t)
Int_m=cumsum(m)*ts;
figure(2)
plot(t,Int_m)
title('Integral de m(t)')

%Modulacion de FM e 
yfm=A*cos(wc*t + kf*Int_m);
figure(3)
plot(t,yfm)
title('Señal modulada FM')


%Espectro de la señal m(t)
M=fftshift(fft(m,N))*ts;
w=linspace(-fs/2, fs/2,N)*2*pi;
figure(4)
plot(w/(2*pi),abs(M))
axis([-100 100 0 1.1*max(abs(M))])
title('Espectro de magnitud de la señal moduladora m(t)')

%Ruido
long_m=length(m);
Py=var(yfm);
SNR=20;
snr=10^(SNR/10);
Pn=Py/snr;
sigma_n=sqrt(Pn);
n=sigma_n*randn(1,long_m);

%Grafica de señal de ruido garussiano
figure(5);
plot (t,n);
title('Ruido gaussiano')
Pnsim=var(m);

%Transformada de ruido
R=fftshift(fft(n,N))*ts;
f=linspace(-fs/2, fs/2, N);
figure(6);
plot(f, abs(R));
%axis([-500 500 0 0.06])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal de ruido');
grid

%Señal de ruido con la señal modulada
yfm_r = yfm + n;

%Modulacion de FM con ruido
figure(7)
plot(t,yfm_r)
title('Señal modulada FM')

%Espectro de la señal modulada
YFM_r=fftshift(fft(yfm_r,N))*ts;
w=linspace(-fs/2, fs/2,N)*2*pi;
figure(8)
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
figure(9)
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

%Convolucion de señal con ruido con filtro pasabandas
r1=conv(yfm_r,h1,'same')*ts;
figure(10)
plot(t,r1);
%axis([-0.08 0.08 -1.5 1.5])
title('Señal pasada por el filtro para quitar ruido')
grid

%Transformada de R1
R1=fftshift(fft(r1,N))*ts;
f=linspace(-fs/2, fs/2, N);
figure(11);
plot(f, abs(R1));
%axis([-700 700 0 .04])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal despues de pasabandas');
grid


%%%%%Demodulacion no coherente con extraccion de fase
%Extrayendo la fase desde la transformada de Hilbert
z=hilbert(r1);


%Se obtiene la fase (Parte imaginaria)
theta=angle(z)-wc*t;
theta1=unwrap(theta);

%Señal integral de m(t)
figure(12)
plot(t,kf*Int_m)
title('kf*Integral de m(t) (Desviacion de fase)')

%Recuperacion de la fase de la tranformda de Hilbert
figure(13)
plot(t, theta1)
title('Recuperacion de la fase de yfm (kf*Integral de m(t))')

%Derivada de la señal de mensaje para obtencion de fase
m_recphase=[0 diff(theta1)/ts]/kf;
figure(14)
plot(t, m_recphase,'r')
axis([0 2 -1.5 1.5])
title('Derivada de la Desviacion de fase obtenida')
grid

%Transformada de R1
M_recphase=fftshift(fft(m_recphase,N))*ts;
f=linspace(-fs/2, fs/2, N);
figure(15);
plot(f, abs(M_recphase),'r');
%axis([-700 700 0 .04])
xlabel('Frecuencia');
ylabel('Magnitud');
title('Espectro de la señal recuperada de la fase');
grid

%Comparacion de señales recuperadas
figure(16)
plot(t, m_recphase,'r')
hold on
plot(t,m,'b')
axis([0 2 -1.5 1.5])
legend('Señal recuperada', 'Señal original ')
title('Comparación de la señal m(t) y la señal recuperada (con filtro)')

%Comparacion de espectros de señal recuperada y original
figure(17)
plot(f, abs(M_recphase),'r');
hold on
plot(f,abs(M),'b')
axis([-100 100 0 1.1*max(abs(M))])
xlabel('Frecuencia');
ylabel('Magnitud');
title('Comparacion de Espectro de la señal recuperada y original');
legend('Señal recuperada', 'Señal original ')
grid
