%****************************************************************************************************
%��ԭ����LLE�����ϵ�������������ϼ���������ЧӦ������ʹ����λ�����ı������
%����������ЧӦ���ú�΢����СL�ı仯
%�������ǻ�ڹ���
%****************************************************************************************************

clc; clear all; close all; clf;
cputime=0;
tic;
%--------------------------------------��������------------------------------------------------
%%------------------------------ע��taoΪǻ������ʱ��-----------------------------------------


ln=1;
i=sqrt(-1);
alpha=0.009; % һȦ���
gamma=1; %������ϵ���� /W/m��
kappa=0.009;  %ǻ�����ϵ��
beta2=-4.711e-26; %2nd order disp. (s2/m)      
pi=3.1415926535;
c=2.99792458e8;           %���� m/s
FSR=226e9; %���ɹ��׷�Χ��Hz�������ɹ��׷�Χ��Ƶ��f�ļ��
tR=1/FSR; %ѭ��һȦ��ʱ��
Q=3e5;%΢��Ʒ������

L=0.628e-3;   %΢���ܳ�
lambda_pump=1550e-9;
w_pump=2*pi*c/lambda_pump;

k0=2*pi/lambda_pump;
neff=Q*kappa*lambda_pump/(pi*L);%��Ч������2.1212
beta=k0*neff;%�������� 8.5987e+06




%%-----------------------------------------------------------------------------------------

k=2^9;% step size��������
tao=((-k/2):(k/2-1))*(tR/k);
w=FSR*2*pi*((-k/2):(k/2-1)); %����Ҷ�任�õ���ԲƵ��
lambda=2*pi*c./(w+w_pump);
lambda=fftshift(lambda);
Po=0.755; Ao=sqrt(Po);                     % input pwr in watts// %Amplitude
w0=0;           % ���ֹⲨ��/Ƶ��
% Ein=Ao*exp(tao*w0*i);
[A]=load('matlab.mat','Ein');
Ein = A.Ein;
% Ein=sqrt(Ao)*(sech(tao));
% Ein = awgn(Ein,140,'measured');                       %�����˹������
%����ĸ�˹������СӰ�����ߵ��߿�Ӱ�����ߵ�Ч��

%--------------------------------------------------------------------------------------- 

w=fftshift(w);
spectrum=fft(Ein); %Pulse spectrum   *2/k�����źŵ���ʵ��ֵ
length(w)
length(spectrum)


%------------------------------------------------------------------------------------
%---------------------------ÿ����1000�Σ������������-------------------------

Nstep =10;
iplot = 1;
zplot(iplot) = 0;  
Uplot(:,iplot) =abs(Ein).^2 ;
Splot(:,iplot) =10*log10((abs(spectrum)*2/k).^2);



%---------------------------------------ǻ������ѭ��--------------------------------------------
%**************************************************************************************************



%**************************************************************************************************

%-------------------------------------------------------------------------------------


zeta0=0.02955506065; %��ʧг��������ѭ��Ȧ���ɵ� -0.0045��0.0653
s_zeta=0;     %ʧг�ٶ� ��λ4e-3ns-1
d_zeta=s_zeta*tR*1e9;  % ÿȦʧгֵd_zeta = 1.7699e-05


%--------------------------------------------------------------------------------------------------
%������ЧӦ������ŷ����������΢ǻ�ڵ��¶ȱ仯��΢ǻ�¶ȵļ��㹫ʽΪ΢�����ȶ�̬����

Ih=Po*kappa*(gamma ^0.5);
K= 2.78e-8;%e-6;%�ȵ�e-5
epsilon=2.7e-6;   %΢��������ϵ������λ��/��
toindex=2.45e-5; %�ȹ�ϵ�����������¶�ϵ����de/dT(��λ��RIU/��0)
a=epsilon+toindex/neff;   %  1.4250e-05
Cp= 8.81e-10;%711.76; %����J/��   8.81e-10
delta_lamb=lambda_pump/Q;%���߿��


nums=12000;  %****�ܷ���Ȧ��*******************************************


delta_t=tR; %tR =4.4248e-12
t=delta_t*(0:1:nums-1); %ѭ��һȦ��ʱ��
delta_T=zeros(1,length(t));
dT=zeros(1,length(t));

power=zeros(1,nums/10);

 jplot = 1;  
 Eplot(1:k) =0 ;



for m=2:1:nums
   
    %-------------------------------------����ʧг��------------------------------------------
    if m<=7893 %�ɵ������״̬ �ڲ�����ЧӦ����ʱ��ѵ���ʧгȦ��Ϊ3946*16*2.5  
        
        zeta=zeta0+m*d_zeta; %��ʧг��������ѭ��Ȧ���ɵ�
        s_pumpw=zeta*FSR;  %����ÿȦ���õ�ʧг�����6.3662e+05
        s_pumpl=-2*pi*c*s_pumpw/(w_pump^2);%ÿ��ѭ��Ȧ�����仯   8.1198e-16
        lamb=lambda_pump+s_pumpl;   %���ڱ��ֲ�������ʧг���ĵ��ڣ�lamb��Ϊ���ں�ı�����ʾ
        
    end
    %---------------------------------------------------------------------------------------
    
    mm=(lamb-lambda_pump*(1+a*(delta_T(m-1))))./(delta_lamb/2);
    nn=(mm^2 )+1;
    dd= K*delta_T(m-1);
    delta_T(m)=delta_T(m-1)+delta_t*(Ih./nn-dd)/Cp ;
    dT(m-1)=delta_T(m)-delta_T(m-1);
    
    L=L+L*epsilon*dT(m-1);
    
    D=exp(FSR*((i*L*beta2*(w.^2)/2)*(tR/2)));  %�����������
    
    spectrum=spectrum.*D ;                     
    f=ifft(spectrum); 
    %N=exp(FSR*(-kappa/2-alpha/2-i*zeta+i*L*gamma*((abs(f)).^2)+kappa^0.5*Ein./f)*tR);   %�������������
    %�����ȹ�������
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
    if (rem(m,Nstep) == 0)%ÿ���Nstep������ص�����
        iplot=iplot+1;
        Uplot(:,iplot)=abs(f).^2;
        Splot(:,iplot) =10*log10((abs(spectrum)*2/k).^2);
        zplot(iplot) = m;   
        %ǻ�ڹ��ʼ��㣨ÿʮȦ����һ�Σ�
        power(iplot)=sum(abs(f).^2);       
    end
     if (m>8000)%ÿ���Nstep������ص�����
        Eplot(jplot:1:jplot+k-1) =f ;
        jplot =jplot+k;
    end
    %n=m;
end

toc


    %-------------------------------------��ͼ------------------------------------------

% figure(1)%΢���¶���ʱ��ı仯
% plot(t,delta_T,'linewidth',2);%axis([0 5e-6 0 5.5]);
% xlabel('t(s)');ylabel('��T(��)');
% 
% figure(2)%ÿȦ�¶ȵı仯
% plot(t,dT,'linewidth',2);
% xlabel('t(s)');ylabel('dT(��)');
colormap(1-jet)

figure(7)
subplot(521);
plot(tao,Uplot(:,1),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(522);
% plot(lambda,Splot(:,35),'b');
fig=stem(lambda,Splot(:,1),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(523);
plot(tao,Uplot(:,35),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(524);
% plot(lambda,Splot(:,35),'b');
fig=stem(lambda,Splot(:,35),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(525);
plot(tao,Uplot(:,150),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(526);
% plot(lambda,Splot(:,150),'b');
fig=stem(lambda,Splot(:,150),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(527);
plot(tao,Uplot(:,250),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(528);
fig=stem(lambda,Splot(:,250),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
% plot(lambda,Splot(:,240),'b');
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
figure(7)
subplot(529);
plot(tao,Uplot(:,600),'r');title('Input'); xlabel('Time'); ylabel('intension');
xlim([-2.2e-12,2.2e-12]);
subplot(5,2,10);
fig=stem(lambda,Splot(:,500),'fill','b','Marker','none');set(fig, 'BaseValue', -150);
% axis([1.194e-6 2.212e-6 -80 max(Splot(:,600))+5]);
xlabel('wavelength(m)');ylabel('spectrum(dB)');
% plot(lambda,Splot(:,600),'b');
xlim([1.2e-6,2.2e-6]);
ylim([-150,20]);
hold on;





figure(3);%ʱ���ݱ�ͼ
pcolor(zplot,tao,Uplot);
Umin=min(min(Uplot));
Umax=max(max(Uplot));
caxis([Umin,Umax])
% colormap(flipud(jet)) % ����ɫӳ�䷽����Ϊ�Ӻ�ɫ����ɫ
% colormap(my_colormap);
colorbar
shading interp
%set(gcf,'color','w')
title('Time Domain Evolution');
% xlabel('s','FontSize',18);
% ylabel('\xi','FontSize',18);
% zlabel('|U|^2','FontSize',14,'Rotation',0,'Position',0);

figure(4);%Ƶ���ݱ�ͼ
H=surface(zplot,lambda,Splot);
Smin=min(min(Splot));
Smax=max(max(Splot));
caxis([Smin,Smax])
% colormap(my_colormap);
% colormap(flipud(jet)) % ����ɫӳ�䷽����Ϊ�Ӻ�ɫ����ɫ
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
spec_out=10*log10((abs(spectrum)*2/k).^2);   %������ת����dB��ʾ

%stem��һ�����Խ���BaseValue�����������k = stem(...)����ͼ��������������ٵ���set(k, 'BaseValue', 10e-2);

fig=stem(lambda,spec_out,'fill','b','Marker','none');set(fig, 'BaseValue', -80);
axis([1.194e-6 2.212e-6 -80 max(spec_out)+5]);
hold on;
xlabel('wavelength(m)');ylabel('spectrum(dB)');
toc;

%ǻ�ڹ��ʱ仯
figure(6);
plot(zplot,power/(max(power)),'r');
