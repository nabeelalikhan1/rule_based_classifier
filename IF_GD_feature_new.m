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
            IF_GD_new(i,iii,c,9)=Thresh;
            signal1=signal;
            [IF,~,~]= FAST_IF_EEG1(fft(signal),195, 4, 4,3,0.0,0);
            IF=IF(:,16:end-16);
            IF_GD_new(i,iii,c,1)=mean(std(IF.'));
            IF_GD_new(i,iii,c,2)=min(diff(sort(mean(IF.'))));
            
            ind=256-round(mean(IF.')*256);
            for i2=1:4
                if and(ind(i2)>3,ind(i2)<254)
                    if max(abs(signal(ind(i2)-3:ind(i2)+3)))>3*IF_GD_new(i,iii,c,9)%std((signal));%*median(abs(signal))
                        spike(i2)=1;
                    else
                        spike(i2)=0;
                    end
                else
                    %                        if abs(signal(ind(i)))>3*std(eeg(c,:))%(signal)
                    %if abs(signal(ind(i)))>3*std(eeg(c,:))%(signal)
                    if max(abs(signal(ind(i2))))>3*IF_GD_new(i,iii,c,9)%std((signal))
                        spike(i2)=1;
                    else
                        spike(i2)=0;
                    end
                end
            end
            IF_GD_new(i,iii,c,3)=sum(spike);
            IF_GD_new(i,iii,c,4)=mean(abs(signal));
            
            if and(IF_GD_new(i,iii,c,3)>3,and(IF_GD_new(i,iii,c,4)>T,and(IF_GD_new(i,iii,c,1)<0.01,IF_GD_new(i,iii,c,2)>0.1)))
                %   if  and(sum(spike)==4,and(mean(abs(signal))>2,and(TF(c)<0.01,TFG(c)>0.1)))
                m_out(iii)=1;
                if mask(iii)==1
                %display('TP');
                end
                %break;
            end
            [IF,S,R]= FAST_IF_EEG1(hilbert(signal),155, 2, 4,3,0.1,0);
                        IF=IF(:,16:end-16);

            %IF_GD_new(i,iii,c,5)=TFHH(c);
            IF_GD_new(i,iii,c,5)=mean(mean(IF));
            IF_GD_new(i,iii,c,6)=mean(std(IF.'));% std(diff(IF).')
             IF_GD_new(i,iii,c,7)=std(diff(IF).');
            [x,~]=size(IF);
            if x>1
                IF_GD_new(i,iii,c,8)=diff(mean(IF.'));
            else
                IF_GD_new(i,iii,c,8)=0.1;
            end
            if and(IF_GD_new(i,iii,c,8)>0.05,and(and(and(IF_GD_new(i,iii,c,5)<0.3,IF_GD_new(i,iii,c,6)<0.01), IF_GD_new(i,iii,c,4)>T),IF_GD_new(i,iii,c,7)<0.01))
                    m_out(iii)=1;
                   if mask(iii)==0
                display('FP');
                   end
              %  break;
            end
        end
        %iii
        %mask(iii)==m_out(iii)
        iii=iii+1;
    end
    
   mm=filter([1 1 1 1 1],1,m_out);
    mm(mm<=2)=0;
    mm(mm>0)=1;
    ind=find(mm==1);
    for i4=1:length(ind)
    if ind(i4)>6
        mm(ind(i4)-4:ind(i4))=1;
    end
        if ind(i4)<length(mask)-2 
        mm(ind(i4):ind(i4)+2)=1;
        end
    end
   
    figure;
    plot(mask);
    hold on;
    stem(mm*0.8,'r');
    sum(mask==mm)/length(mask)
    TA(i)=sum(mask==mm)/length(mask);
    indd=find(mask==1);
    TP(i)=sum(mm(indd)==1);

    sen(i)=sum(mm(indd)==1)/sum(mask==1);
    indd=find(mask==0);
    spec(i)=sum(mm(indd)==0)/sum(mask==0);
    TN(i)=sum(mm(indd)==0);
    N(i)=sum(mask==0);
    P(i)=sum(mask==1);
    Acc(i)=sum(mask==mm);
end
sensitivity=sum(TP)/sum(P)
specificity=sum(TN)/sum(N)
Total_accuracy=sum(Acc)/(sum(P)+sum(N))
save('IF_GD_new','IF_GD_new');  