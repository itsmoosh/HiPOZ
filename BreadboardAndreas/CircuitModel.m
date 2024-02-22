clear, clf ,clc 
lwdth=2; %linewidth for plot

C1=47e-12; %47pF  These are the values which we used in the last experiment:
C2=10e-9; %10nF
R0=1e3; %1kΩ                    
R1=100e3; %100kΩ                        
R2=10e3; %10kΩ


%                   --R1--    --R2--
%           o--R0--|      |--|      |--o
%                   --C1--    --C2--



f=logspace(2,6,1000); %create frequency sweep according to the setup specifications. Only number was increased to receive a smoother plot
Z=R0+1./(1/R1+1i*2*pi*f*C1)+1./(1/R2+1i*2*pi*f*C2); %sum of impedances according to the circuit above



subplot(1,3,1) %Bode Magnitude
plot(f,abs(Z),'LineWidth',lwdth)
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ylim([1e3 1e5]) 
box on
xlabel('f in Hz ');
ylabel('|Z| in \Omega');
title('','Bode Magnitude')

subplot(1,3,2) %Bode phase plot
plot(f,-angle(Z)/pi*180,'LineWidth',lwdth)
set(gca, 'XScale', 'log')
box on
xlabel('f in Hz ');
ylabel('-phase in \circ');
ylim([0 90])
title(strcat('R0=',num2str(R0/1000),'k\Omega, R1=',num2str(R1/1000),'k\Omega, R2=',num2str(R2/1000),'k\Omega, C1=',num2str(C1*1e12),'pF, C2=',num2str(C2*1e9),'nF'),'Bode Phase'); %generate title according to inserted values. comment: this is very unflexible, however, I chose this simple approach, since we have only two capacitor and it is not planned to do this experiment whith more capacitor values. 

subplot(1,3,3) %Nyquist Plot
plot(real(Z)*1e-3,-imag(Z)*1e-3,'LineWidth',lwdth)
grid on
box on
xlabel('real(Z) in k\Omega');
ylabel('-imag(Z) in k\Omega');
title('','Nyquist Plot')


%export Figure
name=strcat('R0_',num2str(R0/1000),'kOhm_R1_',num2str(R1/1000),'kOhm_R2_',num2str(R2/1000),'kOhm_C1_',num2str(C1*1e12),'pF_C2_',num2str(C2*1e9),'nF');
fig = figure(1);
fig.Units = 'centimeters';
fig.Position(3:4) = [40 40/3];
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 12)
exportgraphics(gcf, strcat(name, '.png'),'ContentType', 'vector','BackgroundColor', 'none');