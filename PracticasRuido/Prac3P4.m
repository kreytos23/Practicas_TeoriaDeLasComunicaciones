clear all
close all

to=.15;
fm=75000;
tm=1/fm;
t=-.4:tm:.4;
n=100000;

%señal modulada
m=(1).*(t>=0 & t<=0.05) + (-2).*(t>0.05 & t<=.1);
figure(1)
plot(t,m)
title('Señal de mensaje m(t)')
axis([-0.01 .11 -2.5 1.5])
grid

%portadora
c=cos(500*pi*t);
figure(2)
plot(t,c)
title('Portadora c(t)')
grid
axis([-0.2 0.2 -2 2])


%MODULACIÓN SINGLE SIDE BAND
%Transformada de Hilbert
mh=imag(hilbert(m));
ssb=m.*c+mh.*sin(500*pi*t)
ssb2=m.*c-mh.*sin(500*pi*t);


%Ruido Gaussiano
long_m=length(m);
Py=var(ssb);
SNR=20;
snr=10^(SNR/10);
Pn=Py/snr;
sigma_n=sqrt(Pn);
N=sigma_n*randn(1,long_m);


%Grafica de señal de ruido
figure(3);
plot (t,N);
title('Ruido gaussiano')
Pnsim=var(m);

%Espectro de ruido
n_ruido=fftshift(fft(N,n))*tm;
f=linspace(-fm/2, fm/2, n);
figure(4);
plot(f, abs(n_ruido));
%axis([-0.000001 0.000001 0 0.0016])
xlabel('Frecuencia[Hz]');
ylabel('Magnitud');
title('Espectro de señal de ruido');
grid

%Transformada de hilbert
figure(5)
plot(t,mh)
title('Transformada de Hilbert mh(t)')
axis([-0.2 .3 -6.5 9.5])
grid

%Grafica de señal modulada en SSB
figure(6)
plot(t,ssb)
title('Señal moduladada s(t) en SSB')
axis([-0.2 .2 -5 5])
grid

%Espectro de magnitud y fase de señal m(t)
T=to;
w0=2*pi/T;  %frecuencia fundamental
Cn1=0;
n0=-20:20;  %Se calculan y grafican solo 41 coeficientes
k=0;
for tt=t
    k=k+1;
    Cn1=Cn1 + 1/T*m(k)*exp(-j*n0*w0*tt)*tm;
end
figure(7)
subplot(211)
stem(n0*w0, abs(Cn1))
axis([-1000 1000 0 0.9])
title('Espectro de Magnitud de m(t)')
grid
subplot(212)
stem(n0*w0,angle(Cn1))
axis([-1000 1000 -4 4])
title('Espectro de Fase de m(t)')
grid



%Espectro y fase de señal mh(t)
T=to;
w0=2*pi/T;  %frecuencia fundamental
Cn2=0;
n0=-20:20;  %Se calculan y grafican solo 41 coeficientes
k=0;
for tt=t
    k=k+1;
    Cn2=Cn2 + 1/T*mh(k)*exp(-j*n0*w0*tt)*tm;
end
figure(8)
subplot(211)
stem(n0*w0, abs(Cn2))
axis([-1000 1000 0 0.9])
title('Espectro de Magnitud de mh(t)')
grid
subplot(212)
stem(n0*w0,angle(Cn2))
axis([-1000 1000 -4 4])
title('Espectro de Fase de mh(t)')
grid

%Transformada de Hilbert con ruido Gaussiano
mh_r= mh + N
ssb_r= ssb + N
ssb2_r= ssb2 + N

%Espectro y fase de señal s(t) LSB con ruido 
T=to;
w0=2*pi/T;  %frecuencia fundamental
Cn3=0;
n0=-20:20;  %Se calculan y grafican solo 41 coeficientes
k=0;
for tt=t
    k=k+1;
    Cn3=Cn3 + 1/T*ssb_r(k)*exp(-j*n0*w0*tt)*tm;
end
figure(9)
subplot(211)
stem(n0*w0, abs(Cn3))
axis([-1000 1000 0 0.035])
title('Espectro de Magnitud de s(t) LSB con ruido')
grid
subplot(212)
stem(n0*w0,angle(Cn3))
axis([-1000 1000 -4 4])
title('Espectro de Fase de s(t) LSB con ruido')
grid



%Espectro y fase de señal s(t) USB
T=to;
w0=2*pi/T;  %frecuencia fundamental
Cn4=0;
n0=-20:20;  %Se calculan y grafican solo 41 coeficientes
k=0;
for tt=t
    k=k+1;
    Cn4=Cn4 + 1/T*ssb2_r(k)*exp(-j*n0*w0*tt)*tm;
end
figure(10)
subplot(211)
stem(n0*w0, abs(Cn4))
axis([-900 900 0 0.008])
title('Espectro de Magnitud de s(t) USB con ruido')
grid
subplot(212)
stem(n0*w0,angle(Cn4))
axis([-1000 1000 -5 5])
title('Espectro de Fase de s(t) USB con ruido')
grid

%Filtrado pasabanda en la entrada del receptor
%Filtro pasa bandas
th=-0.2:tm:0.2;
W=300*pi;
wc=500*pi;
w=linspace(-fm/2, fm/2,n)*2*pi;


%Grafica de filtro pasabanda
h1=2*W/pi*sinc(W*th/pi).*cos(wc*th);
figure(11)
subplot(211)
plot(th,h1,'m')
axis([-0.03 0.03 -340 1.1*max(abs(h1))])
title('Filtro pasabanda')

%Espectro de filtro pasabandas
H1=fftshift(fft(h1,n))*tm;
subplot(212)
plot(w/(2*pi),abs(H1),'m')
axis([-600 600 0 1.1*max(abs(H1))])
title('Espectro de Filtro pasabanda')
grid

%Convolucion de la señal con ruido y el filtro pasabanda
ssb_sn=conv(ssb_r,h1,'same')*tm;

%demodulacion de la señal s(t)
r=ssb_sn.*c;
figure(12)
plot(t,r)
title('r(t)=s(t)*c(t)')
axis([-0.1 0.25 -3 2.5])
grid

%Espectro de la demodulacion
n=100000;
R=fftshift(fft(r,n))*tm;
f=linspace(-fm/2, fm/2,n);
figure(13)
plot(f,abs(R))
axis([-750 750 0 1.1*max(abs(R))])
title('Espectro de magnitud de r(t)')
grid

%Filtro pasabajas
%h=200*sinc(200*t);
h=(200*pi)*sinc(200*pi*t);
figure(14)
plot(t,h,'m')
axis([-0.15 0.15 -120 650])
title('Filtro h(t)')
grid

%Espectro del filtro
H=fftshift(fft(h,n))*tm;
figure(15)
plot(f,abs(H),'m')
axis([-450 450 0 1.2])
title('Espectro de magnitud H')
grid

%Grafica de señal recuperada
m_rec=conv(r,h, 'same')*tm;
figure(16)
plot(t,2*m_rec,'r')
title('Señal recuperada')
axis([-0.01 .11 -2.5 1.5])
grid

%Espectro de señal recuperada
M_rec=fftshift(fft(m_rec,n))*tm;
M=fftshift(fft(m,n))*tm;
figure(17)
plot(f,abs(M_rec),'r')
axis([-800 800 0 .12])
title('Espectro de magnitud H')
grid


%Grafica de señal recuperada
figure(18)
plot(t,2*m_rec,'r')
hold on
plot(t,m,'b');
title('Comparacion de Señales recuperada y Original')
legend('Señal recuperada', 'Señal original ')
axis([-0.01 .11 -2.5 1.5])
grid

%Espectro de señal recuperada
figure(19)
plot(f,abs(M_rec),'r')
hold on
plot(f,abs(M),'b')
axis([-800 800 0 .12])
title('Comparacion de Espectro de señales recuperada y original')
legend('Señal recuperada', 'Señal original ')
grid