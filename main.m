clc;
clear all;
clc;

[s,fs]=audioread("EQ2401Project2data2024.wav");

D = 0;
delay = dsp.Delay(D);
x = delay(s);

%%lsm
mu = 0.01; %learning rate
order = 256;
w = lms(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
audiowrite('lms1.wav', sv,fs);

%%rls
order = 512;
lambda = 0.999999;%decay rate
w = rls(s,x,order, lambda);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
audiowrite('rls.wav', sv,fs);

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
audiowrite('nlms.wav', sv,fs);


%%momentumlsm
mu = 0.01; 
order = 256;
ga = 0.8;%decay rate
w = mlms(s,x,mu,order,ga);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
audiowrite('momentumlms.wav', sv,fs);

%%adagradelsm
mu = 0.05; 
order = 256;
w = adag(s,x,mu,order);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
audiowrite('adagradelms.wav', sv,fs);

%%adamlsm
mu = 0.8; 
b1 = 0.2;%decay rate of the first moment
b2 = 0.0000001;%decay rate of the second moment
order = 256;
w = adam(s,x,mu,order,b1,b2);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
audiowrite('adamlms.wav', sv,fs);





