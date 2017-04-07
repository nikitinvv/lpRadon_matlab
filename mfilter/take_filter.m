function [filter] = take_filter(Ns,filter)
os=4;d=1/2;
Nse=os*Ns;
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
filter=wfa;
end