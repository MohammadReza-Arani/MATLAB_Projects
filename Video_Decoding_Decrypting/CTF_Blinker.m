clear; clc;
Video = VideoReader("blinker.mp4");
%% 
nFrames   = Video.NumFrames ;
Frame_Rate = Video.FrameRate;

tFrame    = (1:nFrames) / Frame_Rate ;
ghostCom  = zeros( nFrames, 1 ) ;

for fId = 1 : nFrames
  grayImage = rgb2gray( read( Video, fId )) ;    
  ghostCom(fId) = sum( grayImage(:) ) ;
end

%%
figure() ;  clf ;
set( gcf, 'Color', 'White', 'Units', 'Normalized', ...
    'OuterPosition', [0, 0.1, 1, 0.6] ) ;
plot( tFrame, ghostCom/max(ghostCom), 'b' ) ;
set( gca, 'YTick', [0, 1] ) ;
xlabel( 'Time [s]' ) ;
%%
message = zeros(1,length(tFrame));
Temp_UP = [find(ghostCom/max(ghostCom)>0.9)];
Temp_Down = [find(ghostCom/max(ghostCom) < 0.9)];

for i=1:length(tFrame)
    if(sum(ismember(Temp_UP,i)))
        message(i) =1;
    else
        message(i) = 0;
    end
end

%%
[r]         =risetime(ghostCom);
FREQUENCY   =numel(r)/60;



%%

% % create a sample vector
% v = [1 2 2 3 3 3 4 5 5 5 5 2 2];
v = message;

% find the difference between consecutive elements
d = diff(v);

% identify the start of each sequence of repeated elements
starts = [1 find(d~=0)+1];

% identify the length of each sequence
lengths = diff([starts, numel(v)+1]);

% extract the repeated elements
repeated_nums = v(starts);

% display the results
disp([repeated_nums' lengths'])

%% 

Sequence_1_0 = [repeated_nums' lengths'];
Sequence_1_0(:,2) = Sequence_1_0(:,2)/10;

Code = [];
for j=1:length(Sequence_1_0)

    for p=1:Sequence_1_0(j,2)
        Spec = Sequence_1_0(j,1);
        Code = [Code,num2str(Spec)];
    end

end

