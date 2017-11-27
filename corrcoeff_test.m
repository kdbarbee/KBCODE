close all
clear all
clc

x1 = linspace(1,100);
m = 2;
b = 0;
y1 = m*x1 + b;
iterations = 5;
summary = zeros(iterations,4);
noise_range_summary = zeros(iterations,2);
r2_summary = zeros(iterations,1);
plot(x1, y1, '-r')
hold on

MM = {'o','*','+','x','.','s','d','p','h'}
CC = {'k','b','r','g','y','c',[.5 .6 .7],[.8 .2 .6]} % Cell array of colros.


for i=1:iterations
    
    noise_band = 0.2*i;
    noise_upr = 1 + noise_band;
    noise_lwr = 1 - noise_band;
	noise_range = [noise_lwr, noise_upr];

    noise = (noise_upr - noise_lwr)*rand(size(y1)) + noise_lwr;
    x1_noisy = x1.*x1;
    y1_noisy = y1.*noise;
    r1 = corrcoef(x1, y1);
    r1n = corrcoef(x1, y1_noisy);
    r1xnyn = corrcoef(x1_noisy, y1_noisy);

    [r2, rmse] = rsquare(y1, y1_noisy)
    noise_range_summary(i,:) = noise_range;
    r2_summary(i,:) = r2
    summary(i,1:2) = noise_range;
    summary(i,3) = noise_band;
    summary(i,4) = r2;
    y1_summary(i,:) = y1_noisy;
    plot(x1, y1_noisy, 'Marker', MM{i},'color', CC{i}, 'LineStyle', 'none')
    
    hold on
end