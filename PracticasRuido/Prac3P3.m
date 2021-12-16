clear all
close all

fs=75000;
tm=1/fs;
t=-0.3:tm:0.3;
t0=0.15;
u=0.85;
w=500*pi;
A = 2 / u;

%Señal de mensaje m(t)
m=(1).*(t>=0 & t<=0.05) + (-2).*(t>0.05 & t<=.1);
figure(1);
plot(t,m);
title('Señal de mensaje m(t)');
axis([-0.01 .11 -2.2 1.5])
xlabel('t');
grid

%Transformada de m(t)
N=1000000;
M=fftshift(fft(m,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(2);
plot(f, abs(M));
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
axis([-800 800 0 .12])
title('Espectro de señal moduladora m(t)');
grid

%Portadora
c=cos(500*pi*t);

%Señal modulada en AM s(t)
yam=(A + m).*c;
figure(3);
plot(t,yam);
title('Señal modulada en AM s(t)');
axis([-0.03 .14 -4 4])
grid


%Transformada de s(t)
s=(A + m).*c;
S=fftshift(fft(s,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(4);
plot(f, abs(S));
axis([-600 600 0 .2])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal modulada en AM S(t)');

%Ruido Gaussiano
long_m=length(m);
Py=var(yam);
SNR=20;
snr=10^(SNR/10);
Pn=Py/snr;
sigma_n=sqrt(Pn);
n=sigma_n*randn(1,long_m);


%Grafica de señal de ruido
figure(5);
plot (t,n);
title('Ruido gaussiano')
Pnsim=var(m);

%Espectro de ruido
n_ruido=fftshift(fft(n,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(6);
plot(f, abs(n_ruido));
axis([-10000 10000 -0.001 0.006])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal de ruido');
grid

%Señal modulada con ruido
r=yam + n;
figure(7)
plot(t,r);
axis([-0.04 .15 -5 5.5])
title('Señal s(t) con ruido')
grid

%Transformada de señal con ruido
R=fftshift(fft(r,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(8);
plot(f, abs(R));
axis([-700 700 0 0.09])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal s(t) con ruido');
grid

%Filtrado pasabanda en la entrada del receptor
%Filtro pasa bandas
th=-0.2:tm:0.2;
W=200*pi;
wc=500*pi;
w=linspace(-fs/2, fs/2,N)*2*pi;


%Grafica de filtro pasabanda
h1=2*W/pi*sinc(W*th/pi).*cos(wc*th);
figure(9)
subplot(211)
plot(th,h1,'m')
axis([-0.03 0.03 -340 1.1*max(abs(h1))])
title('Filtro pasabanda')

%Espectro de filtro pasabandas
H1=fftshift(fft(h1,N))*tm;
subplot(212)
plot(w/(2*pi),abs(H1),'m')
axis([-600 600 0 1.1*max(abs(H1))])
title('Espectro de Filtro pasabanda')
grid

%Señal pasada por el pasabandas para eliminar ruido
r1=conv(r,h1,'same')*tm;
figure(10)
plot(t,r1);
axis([-0.04 .15 -5 5.5])
title('Señal pasada por el filtro para quitar ruido')
grid

%Transformada de la señal pasada por el filtro pasabandas
R1=fftshift(fft(r1,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(11);
plot(f, abs(R1));
axis([-700 700 0 0.09])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal sin ruido');
grid

%Demodualcion no coherente
%Hilbert
z=hilbert(r1);
v=abs(z);
m_rec = v - A;
figure(12)
plot(t,m_rec,'r')
axis([-0.01 .11 -2.5 1.5])
title('Señal demodulada con demodulacion no coherente')
grid

%Espectro de la señal recuperada 
M_rec=fftshift(fft(m_rec,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(13)
plot(f, abs(M_rec), 'r');
axis([-800 800 0 .12])
title('Espectro de la señal recuperada')
grid

%Comparacion de señales recuperada y original
figure(14)
plot(t,m_rec,'r')
hold on
plot(t,m,'b')
axis([-0.01 .11 -2.5 1.5])
title('Señal recuperada y señal original')
legend('Señal recuperada', 'Señal original ')
grid

%Espectro de la señal recuperada 
M_rec=fftshift(fft(m_rec,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(15)
plot(f, abs(M_rec), 'r');
hold on
plot(f, abs(M),'b');
axis([-150 150 0 .12])
title('Espectro de la señal recuperada')
legend('Señal recuperada', 'Señal original ')
grid

