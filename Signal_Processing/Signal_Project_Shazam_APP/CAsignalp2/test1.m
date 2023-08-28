clc;
clear;
song_list = get_mp3_list('C:\Users\Mohammad Reza\Desktop\CAsignalp2\Songs\Clips');
L=length(song_list(:,1));

for j=1:L
    Clip = song_list{j};


% Clip = 'Khob_Shod-4.mp3';

SNR=-15:3:15;


for i=1:size(SNR,2)
 disp(i);   
match_clip(Clip,SNR(1,i));

end

end


