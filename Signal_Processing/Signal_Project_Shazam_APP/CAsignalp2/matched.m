function rightous_pairs=matched(pairs)

rightous_pairs(:,1)=pairs(:,3);
rightous_pairs(:,2)=pairs(:,4);
rightous_pairs(:,3)=pairs(:,1);
rightous_pairs(:,4)=pairs(:,2)-pairs(:,1);
end
