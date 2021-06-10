clear all
close all

fs=75000;
ts=1/fs;
%Señal m(t)
t0=0.15;
t=0:ts:t0;
m=1*(0<=t & t<=t0/3) -2*(t>t0/3 & t<=2*t0/3);


%Espectro de magnitud y fase de la señal m(t), periodica
%Usando los coeficientes de la serie de Fourier
T=t0;
w0=2*pi/T;  %frecuencia fundamental
Cn1=0;
n=-20:20;  %Se calculan y grafican solo 41 coeficientes
k=0;
for tt=t
    k=k+1;
    Cn1=Cn1 + 1/T*m(k)*exp(-j*n*w0*tt)*ts;
end
figure(1)
subplot(311)
plot(t,m,'linewidth',2)
axis([0,t0,1.1*min(m), 1.1*max(m)])
grid
subplot(312)
stem(n*w0, abs(Cn1))
grid
subplot(313)
stem(n*w0,angle(Cn1))
grid

%Transformada de hilber de m(t), Espectros de magnitud y fase
mh=imag(hilbert(m))
k=0;
Cn2=0
for tt=t
    k=k+1;
    Cn2=Cn2 + 1/T*mh(k)*exp(-j*n*w0*tt)*ts;
end
figure(2)
subplot(311)
plot(t,mh,'linewidth',2)
axis([0,t0,1.1*min(mh), 1.1*max(mh)])
grid
subplot(312)
stem(n*w0, abs(Cn2))
grid
subplot(313)
stem(n*w0,angle(Cn2))
grid
