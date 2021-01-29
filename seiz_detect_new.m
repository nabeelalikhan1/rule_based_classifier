clear all;
close all;
addpath('E:\Qatar data 30 july 2015\D\QU_EEG DATA\dataset');
addpath('E:\TFSA7\TFSA7')
filt_leng=3;
%load('TT');
load('IF_GD');
IF_GD_new=IF_GD;
IF_GD=IF_GD_new;
load('masked');
TPP=0;
FPP=0;
        T=5 ;  % 84.5% total accuracy 78% specificity 78% sensitivity
       T=4;% 83% total accuracy, specifity==83%, Sensitivity=83%
for i=1:36
    mask=masked{i};
    
    
    
    iii=1;
    ST=8;
 
   m_out=zeros(size(mask));
    
    for j=1:1*256:length(mask)*256-256
        for c=1:20
            
            if and(IF_GD_new(i,iii,c,3)>3,and(IF_GD_new(i,iii,c,4)>T,and(IF_GD_new(i,iii,c,1)<0.015,IF_GD_new(i,iii,c,2)>0.1)))
                %   if  and(sum(spike)==4,and(mean(abs(signal))>2,and(TF(c)<0.01,TFG(c)>0.1)))
                m_out(iii)=m_out(iii)+1*1;
                break;
            end
            
            
            if and(and(and(IF_GD_new(i,iii,c,5)<0.14,IF_GD_new(i,iii,c,6)<0.01), IF_GD_new(i,iii,c,4)>T),IF_GD_new(i,iii,c,7)<0.01)
                    m_out(iii)=m_out(iii)+1*1;
                    if mask(iii)==1
                        %display('IF_TP')
                   %     break;
                   TPP=TPP+1;
                    end
                    if mask(iii)==0
                        FPP=FPP+1;
                       % display('IF_FP')
                   %     break;
                    end
                    
              %      display('IF');
                    break;
            end
        end
       % iii
       if m_out(iii)>=1
           m_out(iii)=1;
       else
           m_out(iii)=0;
       end
      %  m_out(iii)==mask(iii)
        iii=iii+1;
    end
   
    mm=filter([1 1 1  1 ],1,m_out);
    
    mm(mm<3)=0;

     % mm(mm<2)=0;

    mm(mm>0)=1;
    %mm(mm>2 )=1;
    ind=find(mm==1);
       %mm=m_out;

    for i4=1:length(ind)
    if ind(i4)>6
        mm(ind(i4)-6:ind(i4))=1;
    end
        if ind(i4)<length(mask)-4 
        mm(ind(i4):ind(i4)+4)=1;
        end
    end

% for i4=1:length(ind)
%     if ind(i4)>5
%         mm(ind(i4)-5:ind(i4))=1;
%     end
%         if ind(i4)<length(mask)-4 
%         mm(ind(i4):ind(i4)+3)=1;
%         end
%     end
    
    
    %figure;
    if i==	1
    plot(mask);
    hold on;
    stem(mm*0.8,'r');
    end
    %sum(mask==mm)/length(mask)
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
Total_accuracy=(sum(TP)+sum(TN))/(sum(P)+sum(N))
