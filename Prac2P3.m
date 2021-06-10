clear all
close all

fm=75000;
tm=1/fm;
t=0:tm:2;
kf = 50*pi;
wc = 2000*pi;

%señal modulada
m=(t).*(t>=0 & t<=1) + (-t + 2).*(t>1 & t<=1.9)  + 0.1*(t>1.9);
figure(1)
plot(t,m)
title('Señal de mensaje m(t)')
axis([-1.1 3 -0.1 1.1])
grid


%Modulacion de señal m(t)
Int_m = cumsum(m)*tm;
yfm= cos(wc*t + kf*Int_m);
figure(2)
plot(t,yfm)
title('Señal modulada en FM s(t)')
xlabel('tiempo (t)');
axis([0.992 1.008 -1.2 1.2])
grid

%Desviacion de fase es la Int de m(t) por kf
figure(3)
plot(t,kf*Int_m)
title('Desviacion de fase')
axis([-0.5 2.5 -10 180])
grid

%Desviacion de frecuencia es la señal  m(t) por kf, se muestra en HZ
figure(4)
kf_hz = kf/(2*pi);
plot(t,kf_hz*m)
title('Desviacion de frecuencia');
axis([-0.5 2.5 -3 27])
grid

%Espectro de magnitud de la señal de mensaje m(t)
n=100000;
M=fftshift(fft(m,n))*tm;
f=linspace(-fm/2, fm/2, n);
figure(5)
plot(f,abs(M))
axis([-100 100 0 .07])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud m(t)')
grid

%Espectro de magnitud de la señal modulada s(t)
S=fftshift(fft(yfm,n))*tm;
figure(6)
plot(f,abs(S))
axis([-1500 1500 0 .2])
xlabel('Frecuencia') 
ylabel('Magnitud')
title('Espectro de magnitud s(t)')
grid

%Espectro de magnitud de la señal modulada s(t)
S=fftshift(fft(yfm,n))*tm;
figure(7)
plot(f,abs(S))
axis([750 1250 0 .02])
xlabel('Frecuencia') 
ylabel('Magnitud')
title('Espectro de magnitud s(t)')
grid


