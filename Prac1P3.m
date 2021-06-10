clear all
close all

to=.15;
fm=90000;
tm=1/fm;
t=-.4:tm:.4;
u = 0.85;
A = 2 / u;

%señal modulada
m=(1).*(t>=0 & t<=0.05) + (-2).*(t>0.05 & t<=.1);
figure(1)
plot(t,m)
title('Señal de mensaje m(t)')
axis([-0.01 .11 -2.2 1.5])
grid

%Espectro de magnitud de la señal de mensaje m(t)
n=100000; %Armonicos
M=fftshift(fft(m,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(2)
plot(f,abs(M))
axis([-800 800 0 .12])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud m(t)')
grid

%Modular la señal multiplicando por portadora y sumando una amplitud
c=cos(500*pi*t);
yAM= (m + A).*c;
figure(3)
plot(t,yAM)
title('Señal modulada en AM s(t)')
axis([-0.03 .14 -4 4])


%Espectro de magnitud de la señal modulada Yam = s(t)
s=(m + A).*c;
S=fftshift(fft(s,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(4)
plot(f,abs(S))
axis([-600 600 0 .2])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud s(t)')
grid

%Demodulación no coherente con tranformada de Hilbert
r=hilbert(yAM);
v = abs(r);
m_rec = v - A;
figure(5)
plot(t,m_rec,'r')
title('Señal demodulada con demodulacion no coherente')
axis([-0.01 .11 -2.2 1.5])
grid

%Espectro de magnitud de la señal demodulada
esp_r = m_rec;
R_esp=fftshift(fft(esp_r,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(6)
plot(f,abs(R_esp),'r')
hold on
plot(f,abs(M),'b');
axis([-800 800 0 .12])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud de señal demodulada')
grid

%Paso 6: Señal demodulada vs señal original
figure(7)
plot(t,m_rec,'r')
hold on
plot(t,m,'b')
title('Señal Recuperada y Señal Original')
legend('Señal recuperada', 'Señal original ')
axis([-0.01 .11 -2.2 1.5])
grid
