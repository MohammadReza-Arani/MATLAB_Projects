
function [peaks,T,F]=voiceprint(data,fs1)
meaned_clip=mean(data,2);
[sample , channel]=size(data);
DC=mean(meaned_clip(:))*ones(sample,1);
DCless_data=(meaned_clip)-(DC);
fs2=8000;
resampled_data=resample(DCless_data,fs2,fs1);
window=64*10^(-3);
nfft=window;
noverlap=32*10^(-3);
[S,F,T]=spectrogram(resampled_data,window*fs2,noverlap*fs2,nfft*fs2,fs2);
%1
% figure();
% imagesc(T,F,abs(S));
 loged_s=log(abs(S));
% xlabel('Time(s)');
% ylabel('frequency');
% legend on;
%2
% figure();
% imagesc(T,F,loged_s);
% xlabel('Time(s)');
% ylabel('frequency');

CSa1=circshift(loged_s , [0,-1]);
CSa2=circshift(loged_s , [0,1]);
CSa3=circshift(loged_s , [-1,-1]);
CSa4=circshift(loged_s , [1,1]);
CSa5=circshift(loged_s , [1,-1]);
CSa6=circshift(loged_s , [1,0]);
CSa7=circshift(loged_s , [-1,0]);
CSa8=circshift(loged_s , [-1,1]);

Piks=((loged_s-CSa1)>0)&((loged_s-CSa2)>0)&((loged_s-CSa3)>0)&((loged_s-CSa4)>0)...
    &((loged_s-CSa5)>0)&((loged_s-CSa6)>0)&((loged_s-CSa7)>0)&((loged_s-CSa8)>0);
%3
% figure();
% imagesc(T,F,Piks);
% colormap(1-gray);
% xlabel('Time(s)');
% ylabel('frequency');


%A=findpeaks(abs(S),fs2,'MinPeakDistance',0.033); this could work if it
%wasn't because of erors

OUR_peaks=Piks.*loged_s;
[m,n]=size(OUR_peaks);
limitation=zeros(m*n,1);

for i=1:n
    limitation(m*(i-1)+1:m*(i))=OUR_peaks(:,i);
end

    limitation=sort(limitation,'descend');
    
    limitness=zeros(m,n);
    
    for j=1:ceil(30*T(end))
        indices=find(OUR_peaks==limitation(j));
        limitness(indices)=1;
    end
    peaks=limitness;
end























