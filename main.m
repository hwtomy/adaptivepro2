clc;
clear all;
clc;

[s,fs]=audioread("EQ2401Project2data2024.wav");
% plot(periodogram(s,hanning(length(s)), fs));

D = 0;
delay = dsp.Delay(D);
x = delay(s);

%%lsm
mu = 0.01; 
order = 512;
w = lms(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
soundsc(sv,fs);

%%rls
mu = 0.6; 
order = 512;
w = rls(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
soundsc(sv,fs);

%%nlsm
mu = 0.999; 
order = 512;
w = nlms(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
soundsc(sv,fs);





