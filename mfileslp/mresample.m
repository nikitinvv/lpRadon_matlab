function  y = mresample( x, p, q, N )
bta=5;

[p,q] = rat( p/q, 1e-12 );  %--- reduce to lowest terms % (usually exact, sometimes not; loses at most 1 second every 10^12 seconds)
if (p==1) && (q==1)
    y = x; 
    return
end
fc = 1/2/max(p,q);
L = 2*N*max(p,q) + 1;
ideal_filter=2*p*fc*sinc(2*fc*(-(L-1)/2:(L-1)/2))/p;
% ideal_filter=firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]);
h = ideal_filter.*kaiser(L,bta)' ;
h = p*h/sum(h);
Lhalf = (L-1)/2;
Lx = size(x, 1);

% Need to delay output so that downsampling by q hits center tap of filter.
nz = floor(q-mod(Lhalf,q));
h = [zeros(1,nz) h];  % ensure that h is a row vector.
Lhalf = Lhalf + nz;

% Number of samples removed from beginning of output sequence 
% to compensate for delay of linear phase filter:
delay = floor(ceil(Lhalf)/q);

% Need to zero-pad so output length is exactly ceil(Lx*p/q).
nz1 = 0;
while ceil( ((Lx-1)*p+length(h)+nz1 )/q ) - delay < ceil(Lx*p/q)
    nz1 = nz1+1;
end
h = [h zeros(1,nz1)];

% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y = upfirdn(x,h,p,q);

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length
% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)
y(1:delay,:) = [];
y(Ly+1:end,:) = [];


