function z=Fourierfunc(a,x,fs,p)
if a==1
    b=0;
else if a==2
        b=-1;
    end
end

   name=p;     
L=length(x);
NFFT=2.^nextpow2(L);
X=fftshift(fft(x,NFFT)/(16*L));
z=X;

if a==1
f=(fs/2)*linspace(b,1,NFFT/2+1);
figure
plot(f,2*abs(X(1:NFFT/2+1)),'DisplayName',name);
title('single sided amplitude');
else if a==2
f=(fs/2)*linspace(b,1,NFFT);
figure
plot(f,abs(X(1:NFFT)),'DisplayName',name);
title('double sided amplitude');
    end
end

xlabel('Frequency(Hz)');
ylabel('|X(f)|');
grid on
legend('show');


