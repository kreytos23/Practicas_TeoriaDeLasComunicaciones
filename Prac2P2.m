clear all
close all

to=.1;
fm=75000;
tm=1/fm;
t=-.4:tm:.4;
kf = 200*pi;
wc = 500*pi;

%señal modulada
m=(sinc(100*t)).*(t>=-to & t<=to);
figure(1)
plot(t,m)
title('Señal de mensaje m(t)')
axis([-0.06 0.06 -1 1.5])
grid


%Modulacion de señal m(t)
Int_m = cumsum(m)*tm;
yfm= cos(wc*t + kf*Int_m);
figure(2)
plot(t,yfm)
title('Señal modulada en FM s(t)')
xlabel('tiempo (t)');
axis([-0.08 0.08 -1.5 1.5])
grid

%Desviacion de fase es la Int de m(t) por kf
figure(3)
plot(t,kf*Int_m)
title('Desviacion de fase')
axis([-0.1 0.1 -2 8])
grid

%Desviacion de frecuencia es la señal  m(t) por kf, se muestra en HZ
figure(4)
kf_hz = kf/(2*pi);
plot(t,kf_hz*m)
title('Desviacion de frecuencia');
axis([-0.06 0.06 -40 110])
grid

%Espectro de magnitud de la señal de mensaje m(t)
n=100000; %Armonicos
M=fftshift(fft(m,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(5)
plot(f,abs(M))
axis([-100 100 0 .014])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud m(t)')
grid

%Espectro de magnitud de la señal modulada s(t)
S=fftshift(fft(yfm,n))*tm;
figure(6)
plot(f,abs(S))
axis([-700 700 0 .04])
xlabel('Frecuencia') 
ylabel('Magnitud')
title('Espectro de magnitud s(t)')
grid

%Espectro de magnitud de la señal moduladora s(t)aumentada
S=fftshift(fft(yfm,n))*tm;
figure(7)
plot(f,abs(S))
axis([-700 100 0 .03])
xlabel('Frecuencia') 
ylabel('Magnitud')
title('Espectro de magnitud s(t) de cerca')
grid

%%%%%Demodulacion no coherente

%Derivando la señal para obtener la señal de mensaje
yfmdif = [0 diff(yfm)/tm];
figure(8)
plot(t,yfmdif)
title('Señal derivada de yfm')
axis([-0.08 0.08 -2500 2500])
grid

%Espectro de la señal derivada
YFMdif=fftshift(fft(yfmdif,n))*tm;
figure(9)
plot(f,abs(YFMdif))
axis([-400 400 0 100])
title('Espectro de magnitud de la señal derivada')
grid

%Sacando la envolvente con la transformada de Hilbert y quitando las
%amplitudes y kf

yfm_h = hilbert(yfmdif);
yfm_abs = abs(yfm_h);
yfm_rec = (yfm_abs/1-wc)/kf;
figure(10)
plot(t,yfm_rec,'r');
axis([-0.06 0.06 -1 1.5]);
title('Señal recuperada de la envolvente')
hold on
grid

%Espectro de la señal recuperada de la envolvente
YFM_rec=fftshift(fft(yfm_rec,n))*tm;
figure(11)
plot(f,abs(YFM_rec),'r')
axis([-100 100 0 .014])
title('Espectro de magnitud de la señal recuperada de la envolvente')
grid

%Comparacion entre m(t) y m^(t)
figure(12)
plot(t,yfm_rec);
hold on
plot(t,m);
axis([-0.06 0.06 -1 1.5]);
title('Comparacion entre m(t) y m¨(t)')
legend('Señal recuperada', 'Señal original ')
grid

%Comparacion entre espectros de m(t) y m^(t)
figure(13)
plot(f,abs(M));
hold on
plot(f,abs(YFM_rec),'r');
axis([-100 100 0 .014]);
title('Comparacion entre espectros de m(t) y m¨(t)')
legend('Señal original', 'Señal recuperada')
grid


