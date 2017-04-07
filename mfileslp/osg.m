function g=osg(aR,theta)
t=linspace(-pi/2,pi/2,100000);
w=aR*cos(t)+(1-aR)+i*aR*sin(t);
% keyboard
g=max(log(abs(w))+log(cos(theta-atan2(imag(w),real(w)))));
end