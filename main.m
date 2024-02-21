clc;
clear all;
clc;

[s,fs]=audioread("EQ2401Project2data2024.wav");
% plot(periodogram(s,hanning(length(s)), fs));

D = 0;
delay = dsp.Delay(D);
x = delay(s);

%%lsm
% mu = 0.01; 
% order = 256;
% ga = 0.1;
% w = lms(s,x,mu,order);
% sn=transpose(s(1:order-1));
% for i=order:length(s)
%     sn1 = s(i:-1:i-order+1)'*w(:,i-1);
%     sn=[sn,sn1];
% end
% sn = sn';
% sv = s-sn;
% soundsc(sv,fs);

% %%rls
% mu = 0.6; 
% order = 512;
% w = rls(s,x,mu,order);
% sn=transpose(s(1:order-1));
% for i=order:length(s)
%     sn1 = s(i:-1:i-order+1)'*w(:,i-1);
%     sn=[sn,sn1];
% end
% sn = sn';
% sv = s-sn;
% soundsc(sv,fs);
% 
% %%nlsm
% mu = 0.999; 
% order = 512;
% w = nlms(s,x,mu,order);
% sn=transpose(s(1:order-1));
% for i=order:length(s)
%     sn1 = s(i:-1:i-order+1)'*w(:,i-1);
%     sn=[sn,sn1];
% end
% sn = sn';
% sv = s-sn;
% soundsc(sv,fs);

%%adagradelsm
% mu = 0.05; 
% order = 256;
% w = adag(s,x,mu,order);
% sn=transpose(s(1:order-1));
% for i=order:length(s)
%     sn1 = s(i:-1:i-order+1)'*w(:,i-1);
%     sn=[sn,sn1];
% end
% sn = sn';
% sv = s-sn;
% soundsc(sv,fs);

%%adamlsm
mu = 0.8; 
b1 = 0.2;
b2 = 0.0000001;
order = 256;
w = adam(s,x,mu,order,b1,b2);
sn=transpose(s(1:order-1));
for i=order:length(s)
    sn1 = s(i:-1:i-order+1)'*w(:,i-1);
    sn=[sn,sn1];
end
sn = sn';
sv = s-sn;
soundsc(sv,fs);





