clear all
close all

to=.1;
fm=75000;
tm=1/fm;
t=-.4:tm:.4;


%señal modulada
m=(sinc(100*t)).*(t>=-to & t<=to);
figure(1)
plot(t,m)
title('Señal de mensaje m(t)')
axis([-0.04 0.04 -1 1.5])
grid

%portadora
c=cos(500*pi*t);
ybsb_sc= m.*c;
figure(2)
plot(t,ybsb_sc)
title('Señal modulada s(t)')
hold on
plot(t,-m,'--')
plot(t,m,'--')
axis([-0.04 0.04 -1.3 1.5])
grid

%Espectro de magnitud de la señal moduladora m(t)
n=100000; %Armonicos
M=fftshift(fft(m,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(3)
plot(f,abs(M))
axis([-100 100 0 .02])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud m(t)')
grid

%Espectro de magnitud de la señal moduladora ybsb_sc = s(t)
s=m.*c;
S=fftshift(fft(s,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(4)
plot(f,abs(S))
axis([-400 400 0 .006])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud s(t)')
grid

%Demodulación
%Paso1: r=s(t).*c(t)
r=ybsb_sc.*c;
figure(5)
plot(t,r)
title('r(t)=s(t)*c(t)')
axis([-0.035 .035 -0.25 1.1])
grid

%Paso 2: Obtener el espectro de la señal r(t)
R=fftshift(fft(r,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(6)
plot(f,abs(R))
axis([-800 800 0 .006])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud r(t)')
grid

%Paso 3: Diseñar el filtro

h=600*sinc(600*t);
figure(7)
plot(t,h,'m')
axis([-.15 0.15 -30 101])
title('Filtro h(t)')
grid

%Paso 4: Espectro del filtro
H=fftshift(fft(h,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la función linspace para poder hacer el vector del mismo tamaño que el M
figure(8);
plot(f,abs(H),'m');
axis([-700 700 -0.1 1.1]);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title('Espectro del filtro h(t)')
grid

%Paso 5: Multiplico el espectro del filtro * el espectro de la señal r(t) Fig 6
%ESPRECTRO DE LA SEÑAL RECUPERADA
m_rec=R.*H;
figure(9)
plot(f,abs(m_rec),'r')
hold on
plot(f, abs(M),'b')
title('Espectro r(t) * Espectro h(t)')
legend('Señal Recuperada', 'Señal original')
axis([-85 85 0 .015])
grid

%Paso 6: Señal recuperada haciendo la convolución
m_rec=conv(r,h, 'same')*tm;
figure(10)
plot(t,2*m_rec,'r')
hold on
plot(t,m,'b')
title('Señal Recuperada y Señal Original')
legend('Señal recuperada', 'Señal original ')
axis([-0.04 .04 -0.25 1.1])
grid
