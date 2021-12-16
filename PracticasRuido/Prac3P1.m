clear all
close all

fs=75000;
tm=1/fs;
t0=0.15;
t=-0.2:tm:0.2;

%Señal de mensaje m(t)
m=1*(0<=t & t<=t0/3) -2*(t>t0/3 & t<=2*t0/3);
figure(1);
plot(t,m);
title('Señal moduladora m(t)');
axis([-0.01 .11 -2.5 1.5])
xlabel('t');
grid

%Señal portadora c(t)
c=cos(500*pi*t);

%Transformada de m(t)
N=100000;
M=fftshift(fft(m,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(2);
plot(f, abs(M));
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
axis([-850 850 0 .115])
title('Espectro de señal moduladora M(t)');
grid

%Modulación DSB-SC
ydsb_sc=m.*c;
figure(3);
plot(t,ydsb_sc)
title('Señal modulada en DSB-SC s(t)');
axis([-0.01 .11 -2.5 2.5])
grid

%Transformada de s(t)
s=m.*c;
S=fftshift(fft(s,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(4);
plot(f, abs(S));
axis([-800 800 0 .06])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal modulada S(t)');
grid

%Ruido
long_m=length(m);
Py=var(ydsb_sc);
SNR=20;
snr=10^(SNR/10); %Relacion señal a ruido
Pn=Py/snr;
sigma_n=sqrt(Pn);
n=sigma_n*randn(1,long_m);
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

%Señal original sumada con ruido Gaussiano
r = ydsb_sc + n;
figure(7)
plot(t,r);
axis([-0.02 .12 -2.6 2.6])
title('Señal modulada DSB-SC con ruido Gaussiano')

%Transformada de señal modulada con ruido
R=fftshift(fft(r,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(8);
plot(f, abs(R));
axis([-1000 1000 0 0.06])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal modulada DSB-SC con ruido');
grid

%Filtrado pasabanda en la entrada del receptor
%Creacion del filtro pasabanda
th=-0.2:tm:0.2;
W=300*pi;
wc=500*pi;
N=100000;
w=linspace(-fs/2, fs/2,N)*2*pi;
h1=2*W/pi*sinc(W*th/pi).*cos(wc*th);

%Grafica de filtro pasabanda
figure(9)
subplot(211)
plot(t,h1,'m')
axis([-0.03 0.03 -340 1.1*max(abs(h1))])
title('Filtro pasabanda')

%Espectro del filtro pasabanda

H1=fftshift(fft(h1,N))*tm;
subplot(212)
plot(w/(2*pi),abs(H1),'m')
axis([-600 600 0 1.1*max(abs(H1))])
title('Espectro del Filtro pasabanda')
grid

%Convolucion de la señal con ruido y el filtro pasabanda
r1=conv(r,h1,'same')*tm;
figure(10)
plot(t,r1);
axis([-0.02 .12 -2.6 2.6])
title('Señal pasda por el filtro pasabandas')
grid

%Transformada de la señal pasada por el filtro pasabandas
R1=fftshift(fft(r1,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(11);
plot(f, abs(R1));
axis([-1000 1000 0 0.06])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal sin ruido');
grid

%Demodulacion coherente 
Dem_r=r1.*c;
figure(12)
plot(t,Dem_r);
axis([-0.02 0.12 -2.3 1.3])
title('Señal demodulada')
grid

%Transformada de la señal demodulada
DEM_R1=fftshift(fft(Dem_r,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(13);
plot(f, abs(DEM_R1));
axis([-800 800 0 0.058])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal demodulada');
grid

%Filtrado de la señal con filtro pasabajas
%Filtro pasabajas
h=300*sinc(300*t);
figure(14)
plot(t,h,'m')
axis([-.2 0.2 -30 101])
title('Filtro H(t)')
grid

%Espectro del filtro 
H=fftshift(fft(h,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(15);
plot(f,abs(H),'m');
axis([-250 250 -.1 1.3]);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title('Espectro del filtro H(t)')
grid

%Paso 6: Señal Recuperada
m_rec=conv(Dem_r,h, 'same')*tm;
figure(16)
plot(t,2*m_rec,'r')
hold on
plot(t,m,'b')
axis([-0.01 .11 -2.5 1.5])
title('Señal Recuperada y Señal Original')
legend('Señal recuperada', 'Señal original ')
grid

%Paso 7: Espectro de señal recuperada
M_REC=fftshift(fft(m_rec,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(17);
plot(f, abs(M),'r');
hold on;
plot(f,abs(M_REC),'b');
axis([-850 850 0 .115]);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title('Espectro de la señal recuperada')
legend('Señal original','Señal recuperada')
grid
