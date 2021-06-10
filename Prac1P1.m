clear all
close all

to=.15;
fm=75000;
tm=1/fm;
t=-.4:tm:.4;


%se�al modulada
m=(1).*(t>=0 & t<=0.05) + (-2).*(t>0.05 & t<=.1);
figure(1)
plot(t,m)
title('Se�al de mensaje m(t)')
axis([-0.01 .11 -2.01 2.01])
grid

%portadora
c=cos(500*pi*t);
ybsb_sc= m.*c;
figure(2)
plot(t,ybsb_sc)
title('Se�al modulada s(t)')
hold on
plot(t,-m,'--')
plot(t,m,'--')
axis([-0.01 .11 -2.01 2.01])

%Espectro de magnitud de la se�al moduladora m(t)
n=100000; %Armonicos
M=fftshift(fft(m,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la funci�n linspace para poder hacer el vector del mismo tama�o que el M
figure(3)
plot(f,abs(M))
axis([-800 800 0 .12])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud m(t)')
grid

%Espectro de magnitud de la se�al moduladora ybsb_sc = s(t)
s=m.*c;
S=fftshift(fft(s,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la funci�n linspace para poder hacer el vector del mismo tama�o que el M
figure(4)
plot(f,abs(S))
axis([-600 600 0 .06])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud s(t)')
grid


%Demodulaci�n
%Paso1: r=s(t).*c(t)
r=ybsb_sc.*c;
figure(5)
plot(t,r)
title('r(t)=s(t)*c(t)')
axis([-0.01 .11 -2.01 2.01])
grid

%Paso 2: Obtener el espectro de la se�al r(t)
R=fftshift(fft(r,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la funci�n linspace para poder hacer el vector del mismo tama�o que el M
figure(6)
plot(f,abs(R))
axis([-900 900 0 .06])
xlabel('Frecuencia [Hz]') 
ylabel('Magnitud')
title('Espectro de magnitud r(t)')
grid

%Paso 3: Dise�ar el filtro

h=300*sinc(300*t);
figure(7)
plot(t,h,'m')
axis([-.25 0.25 -30 101])
title('Filtro h(t)')
grid

%Paso 4: Espectro del filtro
H=fftshift(fft(h,n))*tm;
f=linspace(-fm/2, fm/2, n); %se usa la "n" en la funci�n linspace para poder hacer el vector del mismo tama�o que el M
figure(8);
plot(f,abs(H),'m');
axis([-200 200 -.1 1.1]);
xlabel('Frecuencia [Hz]') ;
ylabel('Magnitud')
title('Espectro del filtro h(t)')
grid

%Paso 5: Multiplico el espectro del filtro * el espectro de la se�al r(t) Fig 6
%ESPRECTRO DE LA SE�AL RECUPERADA
m_rec=R.*H;
figure(9)
plot(f,abs(m_rec),'r')
hold on
plot(f, abs(M),'b')
title('Espectro r(t) * Espectro h(t)')
legend('Se�al Recuperada', 'Se�al original')
axis([-100 100 0 .12])
grid

%Paso 6: Se�al recuperada haciendo la convoluci�n
m_rec=conv(r,h, 'same')*tm;
figure(10)
plot(t,2*m_rec,'r')
hold on
plot(t,m,'b')
title('Se�al Recuperada y Se�al Original')
legend('Se�al recuperada', 'Se�al original ')
axis([-0.1 .12 -2.1 2.1])
grid

