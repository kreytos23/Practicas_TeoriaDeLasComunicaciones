clear all
close all

to=.1;
fs=75000;
tm=1/fs;
t=-.4:tm:.4;
kf = 100*pi;
wc = 500*pi;

%señal modulada
m=(sinc(100*t)).*(t>=-to & t<=to);
figure(1)
plot(t,m);
title('Señal m(t)');
xlabel('tiempo (t)');
axis([-0.06 0.06 -1 1.5])
grid

%Espectro de magnitud de la señal de mensaje m(t)
N=100000;
M=fftshift(fft(m,N))*tm;
f=linspace(-fs/2, fs/2, N); 
figure(2)
plot(f,abs(M))
axis([-100 100 0 .014])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud m(t)')
grid

%Integral de m(t) y modulacion en FM
Int_m = cumsum(m)*tm;
yfm= cos(wc*t + kf*Int_m);


%Desviacion de fase es la Int de m(t) por kf
figure(3)
plot(t,kf*Int_m)
title('Desviacion de fase de señal modulada')
axis([-0.04 0.04 -0.5 4])
grid

%Desviacion de frecuencia es la señal  m(t) por kf, se muestra en HZ
figure(4)
kf_hz = kf/(2*pi);
plot(t,kf_hz*m)
title('Desviacion de frecuencia de la señal modulada');
axis([-0.06 0.06 -40 110])
grid


%señal m(t) modulada
figure(5)
plot(t,yfm)
title('Señal modulada en FM s(t)')
xlabel('tiempo (t)');
axis([-0.08 0.08 -1.5 1.5])
grid

%Espectro de magnitud de la señal modulada s(t)
S=fftshift(fft(yfm,N))*tm;
figure(6)
plot(f,abs(S))
axis([-700 700 0 .04])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud S(t)')
grid

%Espectro de magnitud de la señal moduladora s(t)aumentada
S=fftshift(fft(yfm,N))*tm;
figure(7)
plot(f,abs(S))
axis([-600 600 0 .04])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud s(t) ')
grid 

%Señal de ruido gaussiano
long_m=length(m);
Py=var(yfm);
SNR=20;
snr=10^(SNR/10); 
Pn=Py/snr;
sigma_n=sqrt(Pn);
n=sigma_n*randn(1,long_m);

%Grafica de señal de ruido garussiano
figure(8);
plot (t,n);
title('Ruido gaussiano')
Pnsim=var(m);

%Transformada de ruido
R=fftshift(fft(n,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(9);
plot(f, abs(R));
%axis([-500 500 0 0.06])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal de ruido');
grid

%Señal modulada mas ruido

r=yfm + n;
figure(10)
plot(t,r);
axis([-0.05 0.05 -2 2])
title('Señal s(t) con ruido')
grid

%Espectro de señal con ruido
R=fftshift(fft(r,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(11);
plot(f, abs(R));
axis([-700 700 0 .04])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal modulada con ruido R(t)');
grid

%Filtrado pasabanda en la entrada del receptor
th=-0.2:tm:0.2;
W=200*pi;
wc=500*pi;
N=100000;
w=linspace(-fs/2, fs/2,N)*2*pi;
h1=2*W/pi*sinc(W*th/pi).*cos(wc*th);
H1=fftshift(fft(h1,N))*tm;

%Grafica de filtro pasabanda
figure(12)
subplot(211)
plot(th,h1,'m')
axis([-0.03 0.03 -340 1.1*max(abs(h1))])
title('Filtro pasabanda')

%Espectro de Filtro pasabandas
subplot(212)
plot(w/(2*pi),abs(H1),'m')
axis([-500 500 0 1.1*max(abs(H1))])
title('Espectro del Filtro pasabanda')
grid

%Convolucion de señal con ruido con filtro pasabandas
r1=conv(r,h1,'same')*tm;
figure(13)
plot(t,r1);
axis([-0.08 0.08 -1.5 1.5])
title('Señal pasada por el filtro para quitar ruido')
grid

%Transformada de R1
R1=fftshift(fft(r1,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(14);
plot(f, abs(R1));
axis([-700 700 0 .04])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal despues de pasabandas');
grid

%%%%%Demodulacion no coherente con extraccion de fase
%Extrayendo la fase desde la transformada de Hilbert
z=hilbert(r1);
wc = 500*pi;

%Se obtiene la fase (Parte imaginaria)
theta=angle(z)-wc*t;
theta1=unwrap(theta);

%Señal integral de m(t)
figure(15)
plot(t,kf*Int_m)
title('kf*Integral de m(t) (Desviacion de fase)')

%Recuperacion de la fase de la tranformda de Hilbert
figure(16)
plot(t, theta1)
title('Recuperacion de la fase de yfm (kf*Integral de m(t))')

%Derivada de la señal de mensaje para obtencion de fase
m_recphase=[0 diff(theta1)/tm]/kf;
figure(17)
plot(t, m_recphase,'r')
axis([-0.06 0.06 -1 1.5])
title('Derivada de la Desviacion de fase obtenida')
grid

%Transformada de R1
M_recphase=fftshift(fft(m_recphase,N))*tm;
f=linspace(-fs/2, fs/2, N);
figure(18);
plot(f, abs(M_recphase),'r');
%axis([-700 700 0 .04])
xlabel('Frecuencia');
ylabel('Magnitud');
title('Espectro de la señal recuperada de la fase');
grid

%Comparacion de señales recuperadas
figure(19)
plot(t, m_recphase,'r')
hold on
plot(t,m,'b')
axis([-0.06 0.06 -1 1.5])
legend('Señal recuperada', 'Señal original ')
title('Comparación de la señal m(t) y la señal recuperada (con filtro)')

%Comparacion de espectros de señal recuperada y original
figure(20)
plot(f, abs(M_recphase),'r');
hold on
plot(f,abs(M),'b')
axis([-100 100 0 .014])
xlabel('Frecuencia');
ylabel('Magnitud');
title('Comparacion de Espectro de la señal recuperada y original');
legend('Señal recuperada', 'Señal original ')
grid