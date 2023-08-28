function  chance = searcher(Clip,hashing,SNR)
    [data,fs] = audioread(Clip);
    data=awgn(data , SNR );
    peaks = voiceprint(data,fs);
    pairs = peak_to_pair(peaks);
    rightous_pairs=matched(pairs);
    L1 = length(rightous_pairs);
    chance = 0;
    
    
    for i = 1 : L1
         hash1 = hash_func(rightous_pairs(i,1),rightous_pairs(i,2),rightous_pairs(i,4));
         p=0;
        for p=1:length(hashing)
             if(hashing(p,1)==hash1)
                 chance = chance+1;
             end
         end 
     end
            