function [g,filter,ghat] = apply_filter(R,os,filter,mu,d,npols)
% Ntheta=numel(th);
[Ns, Ntheta]=size(R);
Nse=os*(Ns);
Re=zeros(Nse,Ntheta);
Re(Nse/2+1+(-floor(Ns/2):ceil(Ns/2)-1),:)=R;

Rhate=fftshift(fft(ifftshift(Re)));


t=(0:(Nse/2))/Nse;
switch filter
    case 'ramp'
        wfa=Nse*0.5*wint(12,t)';%.*(t/(2*d)<=1);%compute the weigths
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        wfa = Nse*0.5*wint(12,t)'.*sinc(t/(2*d)).*(t/d<=2);
    case 'cosine'
        wfa = Nse*0.5*wint(12,t)'.*cos(pi*t./(2*d)).*(t/d<=1);
    case 'cosine2'
        wfa = Nse*0.5*wint(12,t)'.*(cos(pi*t./(2*d))).^2.*(t/d<=1); 
    case 'hamming'
        wfa = Nse*0.5*wint(12,t)'.*(.54 + .46 * cos(pi*t./d)).*(t/d<=1);
    case 'hann'
        wfa=Nse*0.5*wint(12,t)'.*(1+cos(pi*t./d)) / 2.*(t/d<=1);
    otherwise
        error(message('images:iradon:invalidFilter'))
end
wfa=wfa.*(wfa>=0);%.*(t/d<=1);
wfamid=2*wfa(1);%2*
wfa=[fliplr(wfa(2:end)) wfamid wfa(2:end)]';
wfa=wfa(1:end-1);

g=fftshift(ifft(ifftshift(diag(wfa)*Rhate)));
g=g(Nse/2+1+(-floor(Ns/2):ceil(Ns/2)-1),:);
g=real(g');
filter=wfa;
t=[fliplr(-t(2:end)) 0 t(2:end-1)]';

% [Ns, Ntheta]=size(R);
% Nse=os*(Ns);
% Re=zeros(Nse,Ntheta);
% Re(Nse/2+1+(-floor(Ns/2):ceil(Ns/2)-1),:)=R;
% 
% Rhate=fftshift(fft(ifftshift(Re)));
% 
% t1=ceil(Nse*mu/(2*pi))+1;%start index
% t=(t1-1:Nse/2)/Nse;
% % keyboard
% wfa=zeros(1,Nse/2+1);
% switch filter
%     case 'ram-lak'
%         wfa(t1:end)=Nse*0.5*wint(12,t)';%.*(t/(2*d)<=1);%compute the weigths
%         % Do nothing
%     case 'shepp-logan'
%         % be careful not to divide by 0:
%         wfa(t1:end) = Nse*0.5*wint(12,t)'.*sinc(t/(2*d)).*(t/d<=2);
%     case 'cosine'
%         wfa(t1:end) = Nse*0.5*wint(12,t)'.*cos(pi*t./(2*d)).*(t/d<=1);
%     case 'cosine2'
%         wfa(t1:end) = Nse*0.5*wint(12,t)'.*(cos(pi*t./(2*d))).^2.*(t/d<=1); 
%     case 'hamming'
%         wfa(t1:end) = Nse*0.5*wint(12,t)'.*(.54 + .46 * cos(pi*t./d)).*(t/d<=1);
%     case 'hann'
%         wfa(t1:end)=Nse*0.5*wint(12,t)'.*(1+cos(pi*t./d)) / 2.*(t/d<=1);
%     otherwise
%         error(message('images:iradon:invalidFilter'))
% end
% wfamid=2*wfa(1);
% 
% % wfamid
% wfa=[fliplr(wfa(2:end)) wfamid wfa(2:end)]';
% wfa=wfa(1:end-1);
% 
% % keyboard%  figure(1);plot(wfa0(end/2+(-end/32+1+10:end/32+1-10)));axis([0 128-20 -0.005 0.013]);  set(gca, 'XTick', []);  set(gca, 'YTick', []);name=['../bilder_weights/w0', '.png'];set(gcf,'Color','w');img = getframe(gcf);imwrite(img.cdata,name);RemoveWhiteSpace([],0,'file',name,'output',name);
% % 
% %  figure(2);plot(wfa(end/2+(-end/32+1+10:end/32+1-10)));axis([0 128-20 -0.005 0.013]);  set(gca, 'XTick', []);set(gca, 'YTick', []);name=['../bilder_weights/w1', '.png'];set(gcf,'Color','w');img = getframe(gcf);imwrite(img.cdata,name);RemoveWhiteSpace([],0,'file',name,'output',name);
% %  figure(3);plot(wfa(200:end-200));axis([0 2*1024-400 -0.02 0.22]);  set(gca, 'XTick', []);set(gca, 'YTick', []);name=['../bilder_weights/w', '.png'];set(gcf,'Color','w');img = getframe(gcf);imwrite(img.cdata,name);RemoveWhiteSpace([],0,'file',name,'output',name);
% ghat=diag(wfa)*Rhate;
% g=fftshift(ifft(ifftshift(diag(wfa)*Rhate)));
% 
% 
% se=(-Nse/2:Nse/2-1)';
% sigma_short=linspace((mu*Nse/(2*pi)),ceil(mu*Nse/(2*pi)),4)'/length(se);
% sigma_short=[-flipud(sigma_short);sigma_short];
% T=exp(-2*pi*i*sigma_short*se');
% g_compensation=real(T'*diag(1/8*[1 3 3 1 1 3 3 1]'*length(se)*(sigma_short(8)-sigma_short(5)).*abs(sigma_short)/2)*T*Re/(size(Re,1)));
% g=g+g_compensation;
% 
% g=g(Nse/2+1+(-floor(Ns/2):ceil(Ns/2)-1),:);filter=wfa;
end