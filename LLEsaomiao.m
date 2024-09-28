%****************************************************************************************************
%在原来的LLE方程上的线性作用算符上加入由于热效应的作用使得相位发生改变的作用
%加入由于热效应作用后微环大小L的变化
%程序加入腔内功率
%****************************************************************************************************

clc; clear all; close all; clf;
cputime=0;
tic;
%--------------------------------------参数设置------------------------------------------------
%%------------------------------注：tao为腔内慢变时间-----------------------------------------


ln=1;
i=sqrt(-1);
alpha=0.009; % 一圈损耗
gamma=1; %非线性系数（ /W/m）
kappa=0.009;  %腔内耦合系数
beta2=-4.711e-26; %2nd order disp. (s2/m)      
pi=3.1415926535;
c=2.99792458e8;           %光速 m/s
FSR=226e9; %自由光谱范围（Hz），自由光谱范围是频率f的间距
tR=1/FSR; %循环一圈的时间
Q=3e5;%微环品质因数

L=0.628e-3;   %微环周长
lambda_pump=1550e-9;
w_pump=2*pi*c/lambda_pump;

k0=2*pi/lambda_pump;
neff=Q*kappa*lambda_pump/(pi*L);%有效折射率2.1212
beta=k0*neff;%传播常数 8.5987e+06




%%-----------------------------------------------------------------------------------------

k=2^9;% step size采样点数
tao=((-k/2):(k/2-1))*(tR/k);
w=FSR*2*pi*((-k/2):(k/2-1)); %傅里叶变换用的是圆频率
lambda=2*pi*c./(w+w_pump);
lambda=fftshift(lambda);
Po=0.755; Ao=sqrt(Po);                     % input pwr in watts// %Amplitude
w0=0;           % 泵浦光波长/频率
% Ein=Ao*exp(tao*w0*i);
[A]=load('matlab.mat','Ein');
Ein = A.Ein;
%Ein=sqrt(Ao)*(sech(T));
% Ein = awgn(Ein,140,'measured');                       %加入高斯白噪声
%加入的高斯噪声大小影响梳线的线宽，影响梳线的效果

%--------------------------------------------------------------------------------------- 

w=fftshift(w);
spectrum=fft(Ein); %Pulse spectrum   *2/k—求信号的真实幅值
length(w)
length(spectrum)


%------------------------------------------------------------------------------------
%---------------------------每运行1000次，保存相关数据-------------------------

Nstep =10;
iplot = 1;
zplot(iplot) = 0;  
Uplot(:,iplot) =abs(Ein).^2 ;
Splot(:,iplot) =10*log10((abs(spectrum)*2/k).^2);



%---------------------------------------腔内作用循环--------------------------------------------
%**************************************************************************************************



%**************************************************************************************************

%-------------------------------------------------------------------------------------


zeta0=-0.0045000000001; %初失谐量，根据循环圈数可调 -0.0045—0.0653
%20231128估计是扫频的问题，试着调一下扫频范围，应该能出结果。
s_zeta=2e-3;     %失谐速度 单位4e-3 ns-1
d_zeta=s_zeta*tR*1e9;  % 每圈失谐值d_zeta = 1.7699e-05


%--------------------------------------------------------------------------------------------------
%加入热效应，利用欧拉方法计算微腔内的温度变化，微腔温度的计算公式为微环内热动态方程

Ih=Po*kappa*(gamma ^0.5);
K= 2.78e-8;%e-6;%热导e-5
epsilon=2.7e-6;   %微环热膨胀系数，单位：/℃
toindex=2.45e-5; %热光系数（折射率温度系数）de/dT(单位：RIU/℃0)
a=epsilon+toindex/neff;   %  1.4250e-05
Cp= 8.81e-10;%711.76; %热容J/℃   8.81e-10
delta_lamb=lambda_pump/Q;%梳线宽度


nums=35000;  %****总仿真圈数*******************************************


delta_t=tR; %tR =4.4248e-12
t=delta_t*(0:1:nums-1); %循环一圈的时间
delta_T=zeros(1,length(t));
dT=zeros(1,length(t));

power=zeros(1,nums/10);

 jplot = 1;  
 Eplot(1:k) =0 ;



for m=2:1:nums
   
    %-------------------------------------调节失谐量------------------------------------------
    if m<=7893 %可调到最好状态 在不加热效应作用时最佳调节失谐圈数为3946*16*2.5  
        
        zeta=zeta0+m*d_zeta; %初失谐量，根据循环圈数可调
        s_pumpw=zeta*FSR;  %根据每圈设置的失谐量算得6.3662e+05
        s_pumpl=-2*pi*c*s_pumpw/(w_pump^2);%每个循环圈波长变化   8.1198e-16
        lamb=lambda_pump+s_pumpl;   %调节泵浦波长，即失谐量的调节，lamb即为调节后的变量表示
        
    end
    %---------------------------------------------------------------------------------------
    
    mm=(lamb-lambda_pump*(1+a*(delta_T(m-1))))./(delta_lamb/2);
    nn=(mm^2 )+1;
    dd= K*delta_T(m-1);
    delta_T(m)=delta_T(m-1)+delta_t*(Ih./nn-dd)/Cp ;
    dT(m-1)=delta_T(m)-delta_T(m-1);
    
    L=L+L*epsilon*dT(m-1);
    
    D=exp(FSR*((i*L*beta2*(w.^2)/2)*(tR/2)));  %线性作用算符
    
    spectrum=spectrum.*D ;                     
    f=ifft(spectrum); 
    %N=exp(FSR*(-kappa/2-alpha/2-i*zeta+i*L*gamma*((abs(f)).^2)+kappa^0.5*Ein./f)*tR);   %非线性作用算符
    %加入热光作用项
    N=exp(FSR*(-kappa/2-alpha/2-i*zeta+i*L*gamma*((abs(f)).^2)+kappa^0.5*Ein./f+i*L*toindex*beta*dT(m-1))*tR); 
    f=f.*N;
    spectrum=fft(f);
    spectrum=spectrum.*D ;
    %---------------------------------------------------------------------------------------
    f=ifft(spectrum);
   % op_pulse(ln,:)=abs(f);%saving output pulse at all intervals
    spectrum=fft(f);
    ln=ln+1;
    %----------------------------------------------------------------------------------------
    if (rem(m,Nstep) == 0)%每间隔Nstep保存相关的数据
        iplot=iplot+1;
        Uplot(:,iplot)=abs(f).^2;
        Splot(:,iplot) =10*log10((abs(spectrum)*2/k).^2);
        zplot(iplot) = m;   
        %腔内功率计算（每十圈计算一次）
        power(iplot)=sum(abs(f).^2);       
    end
     if (m>8000)%每间隔Nstep保存相关的数据
        Eplot(jplot:1:jplot+k-1) =f ;
        jplot =jplot+k;
    end
    %n=m;
end

toc


    %-------------------------------------画图------------------------------------------

% figure(1)%微环温度随时间的变化
% plot(t,delta_T,'linewidth',2);%axis([0 5e-6 0 5.5]);
% xlabel('t(s)');ylabel('ΔT(℃)');
% 
% figure(2)%每圈温度的变化
% plot(t,dT,'linewidth',2);
% xlabel('t(s)');ylabel('dT(℃)');


figure(7)
subplot(421);
plot(tao,Uplot(:,1),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(422);
% plot(lambda,Splot(:,35),'b');
fig=stem(lambda,Splot(:,1),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(423);
plot(tao,Uplot(:,160),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(424);
% plot(lambda,Splot(:,150),'b');
fig=stem(lambda,Splot(:,160),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(425);
plot(tao,Uplot(:,200),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(426);
fig=stem(lambda,Splot(:,200),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
% plot(lambda,Splot(:,240),'b');
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(427);
plot(tao,Uplot(:,300),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(428);
fig=stem(lambda,Splot(:,300),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
% axis([1.194e-6 2.212e-6 -80 max(Splot(:,600))+5]);
xlabel('wavelength(m)');ylabel('spectrum(dB)');
% plot(lambda,Splot(:,600),'b');
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;

figure(8)
subplot(421);
plot(tao,Uplot(:,350),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(422);
% plot(lambda,Splot(:,35),'b');
fig=stem(lambda,Splot(:,350),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(8)
subplot(423);
plot(tao,Uplot(:,500),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(424);
% plot(lambda,Splot(:,150),'b');
fig=stem(lambda,Splot(:,500),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(8)
subplot(425);
plot(tao,Uplot(:,1000),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(426);
fig=stem(lambda,Splot(:,1000),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
% plot(lambda,Splot(:,240),'b');
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(8)
subplot(427);
plot(tao,Uplot(:,1500),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(428);
fig=stem(lambda,Splot(:,1500),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
% axis([1.194e-6 2.212e-6 -80 max(Splot(:,600))+5]);
xlabel('wavelength(m)');ylabel('spectrum(dB)');
% plot(lambda,Splot(:,600),'b');
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;







figure(3);%时域演变图
pcolor(zplot,tao,Uplot);
Umin=min(min(Uplot));
Umax=max(max(Uplot));
caxis([Umin,Umax])
colorbar
shading interp
%set(gcf,'color','w')
title('Time Domain Evolution');
% xlabel('s','FontSize',18);
% ylabel('\xi','FontSize',18);
% zlabel('|U|^2','FontSize',14,'Rotation',0,'Position',0);

figure(4);%频域演变图
H=surface(zplot,lambda,Splot);
Smin=min(min(Splot));
Smax=max(max(Splot));
caxis([Smin,Smax])
colorbar
shading interp
axis([0 m 1.194e-6 2.212e-6]);
title('Frequency Domain Evolution');
% xlabel('s','FontSize',18);
% ylabel('\xi','FontSize',18);
% zlabel('|U|^2','FontSize',14,'Rotation',0,'Position',0);



figure(5);
subplot(2,1,1);
plot(tao,abs(f).^2,'r');
xlabel('Time');ylabel('intension');title('Output');
subplot(2,1,2);
spec_out=10*log10((abs(spectrum)*2/k).^2);   %将功率转换成dB表示

%stem有一个属性叫做BaseValue。你可以先用k = stem(...)画出图，调整好坐标后，再调用set(k, 'BaseValue', 10e-2);

fig=stem(lambda,spec_out,'fill','b','Marker','none');set(fig, 'BaseValue', -80);
axis([1.194e-6 2.212e-6 -80 max(spec_out)+5]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
toc;

%腔内功率变化
figure(6);
plot(zplot,power/(max(power)),'r');
