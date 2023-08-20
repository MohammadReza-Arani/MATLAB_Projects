
function text=TIME_CALENDER()
clc;

 date = datestr(now, 'dd/mm/yy-HH:MM');
%  set(handles.edit,'String',date);
TIME=datestr(now,'dd/mm/yy-HH:MM');

miladi_date(1)= 2000+str2num(TIME(7:8)) ;    %input('enter the year of miladi date:');
miladi_date(2)=  str2num(TIME(4:5));    %input('enter the month of miladi date:');

while miladi_date(2)>=13||miladi_date(2)<=0%for making sure that the input is in the range%
    disp('input is out of range')
    disp('press any key to continue')
    pause
    miladi_date(2)=input('enter the month of miladi date:');
end

miladi_date(3)=   str2double(TIME(1:2))    ;%input('enter the day of miladi date ');
disp("Miladi date is :");
datestr(now,'dd/mm/yy-HH:MM')
while miladi_date(3)>=32||miladi_date(3)<=0
    disp('input is out of range')
    disp('press any key to continue')
    pause
    miladi_date(3)=input('enter the day of miladi date');
end

% date equivalent for program 1385/10/11=2007/1/1
% This program works for miladi year more than 2007
% This program doesn't work for leap years

if miladi_date(2)==2%for making sure that the input isn't leap year%
    while miladi_date(3)>=29||miladi_date(3)<=0
       disp('this program does''nt work for leap years and max day for Febroery is 28')
       disp('press any key to continue')
       pause
       miladi_date(3)=input('enter the day for miladi date');
    end
end
extra_day=(miladi_date(1)-2007)*365;%for computing how many years have pased since the 2007/1/1%
switch miladi_date(2)%for adding the day less than a year that have passed from 2007/1/1 to extra variable%
    case 1
       extra_day=extra_day+miladi_date(3);
    case 2
       extra_day=extra_day+31+miladi_date(3);%31 for january%
    case 3
       extra_day=extra_day+31+28+miladi_date(3);%31+28 for january and february%
    case 4
       extra_day=extra_day+31+28+31+miladi_date(3);%31+28+31 for january and february and march%
    case 5
       extra_day=extra_day+31+28+31+30+miladi_date(3);%adding the month in order in the following cases%
    case 6
       extra_day=extra_day+31+28+31+30+31+miladi_date(3);
    case 7
       extra_day=extra_day+31+28+31+30+31+30+miladi_date(3);
    case 8
       extra_day=extra_day+31+28+31+30+31+30+31+miladi_date(3);
    case 9
       extra_day=extra_day+31+28+31+30+31+30+31+31+miladi_date(3);
    case 10
       extra_day=extra_day+31+28+31+30+31+30+31+31+30+miladi_date(3);
    case 11
       extra_day=extra_day+31+28+31+30+31+30+31+31+30+31+miladi_date(3);
    case 12
       extra_day=extra_day+31+28+31+30+31+30+31+31+30+31+30+miladi_date(3);
end
shamsi_date(1)=1385;
while extra_day>=365
    shamsi_date(1)=shamsi_date(1)+1;
    extra_day=extra_day-365;
end
if extra_day<=19
    shamsi_date(2)=10;
    shamsi_date(3)=extra_day;
elseif extra_day<=(19+30)
    shamsi_date(2)=11;
    shamsi_date(3)=extra_day-19;
elseif extra_day<=(19+30+29)
    shamsi_date(2)=12;
    shamsi_date(3)=extra_day-(19+30);
elseif extra_day<=(19+30+29+31)
    shamsi_date(2)=1;
    shamsi_date(3)=extra_day-(19+30+29);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+2*31)
    shamsi_date(2)=2;
    shamsi_date(3)=extra_day-(19+30+29+31);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+3*31)
    shamsi_date(2)=3;
    shamsi_date(3)=extra_day-(19+30+29+2*31);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+4*31)
    shamsi_date(2)=4;
    shamsi_date(3)=extra_day-(19+30+29+3*31);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+5*31)
    shamsi_date(2)=5;
    shamsi_date(3)=extra_day-(19+30+29+4*31);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+6*31)
    shamsi_date(2)=6;
    shamsi_date(3)=extra_day-(19+30+29+5*31);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+6*31+30)
    shamsi_date(2)=7;
    shamsi_date(3)=extra_day-(19+30+29+6*31);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+6*31+2*30)
    shamsi_date(2)=8;
    shamsi_date(3)=extra_day-(19+30+29+6*31+30);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+6*31+3*30)
    shamsi_date(2)=9;
    shamsi_date(3)=extra_day-(19+30+29+6*31+2*30);
    shamsi_date(1)=shamsi_date(1)+1;
elseif extra_day<=(19+30+29+6*31+3*30+10)
    shamsi_date(2)=10;
    shamsi_date(3)=extra_day-(19+30+29+6*31+3*30);
    shamsi_date(1)=shamsi_date(1)+1;
end
shamsi_date(3)=shamsi_date(3)-1;
disp('the shamsi date  is: ')
disp([num2str(shamsi_date)])
text=shamsi_date;
datetime(now,'ConvertFrom','datenum')
end