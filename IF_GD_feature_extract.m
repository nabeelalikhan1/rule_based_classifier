clear all;
close all;
addpath('E:\Qatar data 30 july 2015\D\QU_EEG DATA\dataset');
addpath('E:\TFSA7\TFSA7')
filt_leng=5;
w=fir1(30,0.05,'high');
maxiter = 5;
residual_energy = 0.1;
%load('IF_GD_new');

for i=1:36
    load([num2str(i) '_data.mat']);
    eegData=data.';
   clear eeg;
       mask=mask(1:8:end);
length(mask)
    for c=1:20
       
        eeg(c,:)=filter(w,1,resample(eegData(c,:),32,256));
        
        eeg(c,:)=filter([1 -1],1,eeg(c,:));
    end
    %eeg=sum(eeg,1);
i    
    mask(eegState(:,3)==5)=1;
    mask(eegState(:,3)==2)=1;
    mask(eegState(:,3)==3)=0;
    mask(eegState(:,3)==4)=0;
    mask(eegState(:,3)==1)=1;
    
    
    %mask=mask(1:end/2);
    %eeg=eeg(:,1:end/2);
    %eeg=filter([1 -1],1,eeg);
    iii=1;
    ST=8;
    %T=2*ones(1,20);
        T=4;
m_out=zeros(size(mask));
Thresh=median(std(eeg(:,1:1+256-1).'));
for j=1:1*256:length(eeg)-256
        for c=1:20
            signal=eeg(c,j:j+256-1);
            %Thresh=std(signal);
            %Thresh=median(std(eeg(:,j:j+256-1).'));
                        Thresh=10*0.1*mean(std(eeg(:,j:j+256-1).'))+0*0.9*Thresh;

            %Thresh=mean(std(eeg(:,j:j+255)));
            %        [~,signal] = mp_Gabor(signal,residual_energy,maxiter);
            %         signal=filter([1 -1],1,signal);
            IF_GD(i,iii,c,8)=Thresh;
            signal1=signal;
            [IF,~,~]= FAST_IF_EEG1(fft(signal),195, 4, 4,3,0.0,0);
            IF=IF(:,16:end-16);
            IF_GD(i,iii,c,1)=mean(std(IF.'));
            IF_GD(i,iii,c,2)=min(diff(sort(mean(IF.'))));
            
            ind=256-round(mean(IF.')*256);
            for i2=1:4
                if and(ind(i2)>3,ind(i2)<254)
                    if max(abs(signal(ind(i2)-3:ind(i2)+3)))>3*IF_GD(i,iii,c,8)%std((signal));%*median(abs(signal))
                        spike(i2)=1;
                    else
                        spike(i2)=0;
                    end
                else
                    %                        if abs(signal(ind(i)))>3*std(eeg(c,:))%(signal)
                    %if abs(signal(ind(i)))>3*std(eeg(c,:))%(signal)
                    if max(abs(signal(ind(i2))))>3*IF_GD(i,iii,c,8)%std((signal))
                        spike(i2)=1;
                    else
                        spike(i2)=0;
                    end
                end
            end
            IF_GD(i,iii,c,3)=sum(spike);
            IF_GD(i,iii,c,4)=mean(abs(signal));
            
            [IF,S,R]= FAST_IF_EEG1(hilbert(signal),155, 2, 4,3,0.1,0);
                        IF=IF(:,16:end-16);

            %IF_GD_new(i,iii,c,5)=TFHH(c);
            IF_GD(i,iii,c,5)=mean(mean(IF));
            IF_GD(i,iii,c,6)=mean(std(IF.'));% std(diff(IF).')
             IF_GD(i,iii,c,7)=std(diff(IF).');
            
        end
        %iii
        %mask(iii)==m_out(iii)
        iii=iii+1;
    end
  
end
save('IF_GD','IF_GD');  