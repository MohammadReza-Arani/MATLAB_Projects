clc;
clear;

song_list = get_mp3_list('C:\Users\Mohammad Reza\Desktop\CAsignalp2\Songs\Train');
for i=1:length(song_list(:,1))
songname=song_list(i,1);
hash_producer(songname{1},num2str(i));
end

%%hash_producer('07_Homayoun_Shajarian_Sohrab_Pournazeri-Norouz.mp3','noroz');

%%hash_producer('11_Homayoun_Shajarian_Sohrab_Pournazeri-Irane_Man.mp3','Iraneman');
