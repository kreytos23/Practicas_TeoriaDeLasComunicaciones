clear all
close all

fs=75000;
tm=1/fs;
t0=0.1;
t=-0.2:tm:0.2;


%Funcion m(t)
m=sinc(100*t).*(abs(t)<=t0);
plot(t,m);
axis([-0.13 0.13 -0.3 1])
title('Señal m(t)');
xlabel('tiempo');
grid

%Transformada de m(t)
N=100000;
M=fftshift(fft(m,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(2);
plot(f, abs(M));
axis([-100 100 0 .014])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal m(t)');
grid

%Señal modulada en DSB-SC
%Portadora
c=cos(500*pi*t);

%Modulacion
ydsb_sc=m.*c;
figure(3)
plot(t,ydsb_sc)
title('Señal modulada s(t)');
axis([-0.04 0.04 -1.3 1.5])
hold on
plot(t,-m,'--')
plot(t,m,'--')
grid


%Transformada de s(t)
s=m.*c;
S=fftshift(fft(s,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(4);
plot(f, abs(S));
axis([-400 400 0 .006])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal modulada S(t)');
grid

%Ruido Gaussiano
long_m=length(m);
Py=var(ydsb_sc);
SNR=10;
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
axis([-2000 2000 -0.0002 0.0012])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal de ruido');
grid


%Grafica de la señal modulada con ruido
r=ydsb_sc + n;
figure(7)
plot(t,r);
axis([-0.04 0.04 -1.3 1.5])
title('Señal s(t) con ruido')

%Transformada de la señal modulada con ruido
R=fftshift(fft(r,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(8);
plot(f, abs(R));
axis([-400 400 0 0.006])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal modulada con ruido');
grid

%Filtrado pasabanda en la entrada del receptor
%Creacion de filtro pasabanda
th=-0.2:tm:0.2;
W=200*pi;
wc=500*pi;
N=100000;
w=linspace(-fs/2, fs/2,N)*2*pi;
h1=2*W/pi*sinc(W*th/pi).*cos(wc*th);

%Grafica de filtro pasabanda
figure(9)
subplot(211)
plot(t,h1,'m')
axis([-0.03 0.03 -300 1.1*max(abs(h1))])
title('Filtro pasabanda')


%Espectro de Filtro pasabandas
H1=fftshift(fft(h1,N))*tm;
subplot(212)
plot(w/(2*pi),abs(H1),'m')
axis([-500 500 0 1.1*max(abs(H1))])
title('Espectro de Filtro pasabanda')
grid


%Señal con pasabandas para quitar ruido
r1=conv(r,h1,'same')*tm;
figure(10)
plot(t,r1);
axis([-0.04 0.04 -1.3 1.5])
title('Señal pasada por el filtro para quitar ruido')
grid

%Transformada de R1
R1=fftshift(fft(r1,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(11);
plot(f, abs(R1));
axis([-400 400 0 0.006])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal sin ruido');
grid

%Demodulacion coherente
Dem_r=r1.*c;
figure(12)
plot(t,Dem_r);
axis([-0.04 0.04 -0.3 1.05])
title('Señal demodulada sin filtrado')
grid

%Transformada de la señal demodulada
DEM_R1=fftshift(fft(Dem_r,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(13);
plot(f, abs(DEM_R1));
axis([-800 800 0 0.055])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal sin ruido recuperada');
axis([-700 700 0 0.006])
grid

%Creacion de Filtro pasa bajas
h=300*sinc(300*t);
figure(14)
plot(t,h,'m')
axis([-.2 0.2 -30 101])
title('Filtro H(t)')
grid

%Espectro del filtro pasa bajas
H=fftshift(fft(h,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(15);
plot(f,abs(H),'m');
axis([-200 200 -.1 1.1]);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title('Espectro del filtro H(t)')
grid

%Convolucion de la señal demodulada con el Filtro
m_rec=conv(Dem_r,h, 'same')*tm;
figure(16)
plot(t,2*m_rec,'r')
hold on
plot(t,m,'b')
axis([-0.05 0.05 -0.4 1.3])
title('Señal Recuperada y Señal Original')
legend('Señal recuperada', 'Señal original ')
grid

%Espectro de señal recuperada
M_REC=fftshift(fft(m_rec,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(17);
plot(f,abs(M_REC),'r');
hold on;
plot(f, abs(M),'b');
axis([-80 80 0 0.013]);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title('Espectro de la señal recuperada')
legend('Señal recuperada', 'Señal original ')
grid
