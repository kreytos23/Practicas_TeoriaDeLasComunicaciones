
f = [1 5 10 15 20 25 30 35 40 45 50 55 60 65 73.45 75 80 85 90 95 100 200 300 400 500 600];

v = [0.165168 0.67248 1.665 2.419 3.353 4.235 5.089 5.983 6.815 7.65 8.397 8.974 9.428 9.762 9.97 9.981 9.827 9.653 9.412 9.106 9.095 4.551 2.991 2.228 1.78 1.482];

plot(f,v,'--ro','linewidth',3,'markersize',7,'markeredgecolor','b')
title("Circuito resonante en serie")
axis([0 650 0 11])
xlabel("frecuencia (kHz)")
ylabel("volts (V)")