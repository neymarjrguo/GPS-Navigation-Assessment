A=importdata('D:\AAE\Assignment\Assignment\Data\eph.dat');
B=importdata('D:\AAE\Assignment\Assignment\Data\rcvr.dat');
A=A([2:end 1],:);
sqrta=A(:,10);%squre root of semi-major axis
sqrta2=sqrta.^2;
Wedot=7.2921151467e-5;%value of earth's rotation rate
c=299792458.0;%speed of light
mu=3.986005e+14;%value of earth's universal gravitation constant
F=(-4.442807633e-10);% Relativistic correction term constant
rcvr_tow=A(:,1);%receiver time of weeks
svid=A(:,2);%satellite PRN number;
toc=A(:,3);%reference time of clocak parameters
toe=A(:,4);%reference time of ephemeris parameters
af0=A(:,5);%clock correction coefficient-group delay
af1=A(:,6);%clock correction coefficient
af2=A(:,7);%clock correction coefficient(s/s/s)
ura=A(:,8);%user range accuracy
e=A(:,9);%eccentricity
dn=A(:,11);%mean motion correction
m0=A(:,12);%mean anomaly at reference time
w=A(:,13);%argument of perigee
omg0=A(:,14);%right ascension
i0=A(:,15);%inclination angle at reference time
odot=A(:,16);%rate of right ascension
idot=A(:,17);%rate of inclination angle
cus=A(:,18);%argument of latitude correction, sine
cuc=A(:,19);%argument of latitude correction,cosine
cis=A(:,20);%inclination of correction, sine
cic=A(:,21);%inclination of correction,cosine
crs=A(:,22);%radius correction,sine
crc=A(:,23);%radius correction,cosine
iod=A(:,24);%issue of data number
t=ones(8,1)*440992;%t=440992
tk=t-toe;% time from ephemeris re
mu=mu.*ones(8,1);%earth's universal gravitation constant
n_0=(mu./sqrta2.^3).^(1/2);%Coputed mean motion
n=n_0+dn;%corrected mean motion
mk=m0+n.*tk;%mean anomaly
for i=1:8;
    if mk(i)>2*pi;
        mk(i)=mk(i)-2*pi;
    else
        if mk(i)<0;
            mk(i)=mk(i)+2*pi;
        end
    end
end
E_old=mk;
E_new=mk+e.*sin(E_old);
i=1;
if abs(E_new-E_old)>1e-8
i=i+1;
E_old=E_new;
E_new=mk+e.*sin(E_old);
end
ek=E_new;%Eccentric Anomaly
cosvk=(cos(ek)-e)./(1-e.*cos(ek));
sinvk=(((1-e.^2).^(1/2)).*sin(ek))./(1-e.*cos(ek));
vk=atan2(sinvk,cosvk);%true anomaly
pk=vk+w;%pk:argument of latitude
uk=pk+cus.*sin(2*pk)+cuc.*cos(2*pk);%uk: Correct Argument of Latitude
rk=sqrta2.*(1-e.*cos(ek))+crs.*sin(2*pk)+crc.*cos(2*pk);%rk:Correct Radius
ik=i0+idot.*tk+cis.*sin(2*pk)+cic.*cos(2*pk);%ik:Correct Inclination
xk1=rk.*cos(uk);
yk1=rk.*sin(uk);
zk1=0;% Position in orbital plane
wedot2=Wedot.*ones(8,1);
omgk=omg0+(odot-wedot2).*tk-wedot2.*toe;%corrected longtude of ascending node
xk=xk1.*cos(omgk)-yk1.*cos(ik).*sin(omgk);
yk=xk1.*sin(omgk)+yk1.*cos(ik).*cos(omgk);
zk=yk1.*sin(ik); % Earth-fixed coordinates
tdk=af0+af1.*(t-toe)+af2.*(t-toe).*(t-toe);% The broadcast satellite clock error
pr=B(:,3);% Pseudorange
pr1=pr+c.*tdk;% range which has fixed with the clock error
%Linearization
xu0=-2694685.473;
yu0=-4293642.366;
zu0=3857878.924;% Initial position
tu0=0;%clcok bias
pr0=((xk-xu0).^2+(yk-yu0).^2+(zk-zu0).^2).^0.5;%Approx pseudorange
ax0=(xu0-xk)./((xk-xu0).^2+(yk-yu0).^2+(zk-zu0).^2).^0.5;
ay0=(yu0-yk)./((xk-xu0).^2+(yk-yu0).^2+(zk-zu0).^2).^0.5;
az0=(zu0-zk)./((xk-xu0).^2+(yk-yu0).^2+(zk-zu0).^2).^0.5;% elevation angle
c0=c.*ones(8,1);
H=[ax0,ay0,az0,c0];%Matrix H
deltp0=pr1-pr0;
deltX=H\deltp0;
x0=[xu0,yu0,zu0,0]';% Initial position
eps=1.0e-4;
count=0;
while abs(deltX(1))>eps||abs(deltX(2))>eps||abs(deltX(3))>eps||abs(deltX(4)>eps);
    count=count+1;
x1=x0+deltX;
pr0=((xk-x1(1,:)).^2+(yk-x1(2,:)).^2+(zk-x1(3,:)).^2).^0.5+c.*x1(4,:);%Approx pseudorange
deltp0=pr1-pr0;
ax1=(x1(1,:)-xk)./((xk-x1(1,:)).^2+(yk-x1(2,:)).^2+(zk-x1(3,:)).^2).^0.5;
ay1=(x1(2,:)-yk)./((xk-x1(1,:)).^2+(yk-x1(2,:)).^2+(zk-x1(3,:)).^2).^0.5;
az1=(x1(3,:)-zk)./((xk-x1(1,:)).^2+(yk-x1(2,:)).^2+(zk-x1(3,:)).^2).^0.5;
count=count+1;
H=[ax1,ay1,az1,c0];
deltX=H\deltp0;
x0=x1;
end
