import builtins
from operator import sub
from numpy.random.mtrand import noncentral_chisquare
from scipy.fft import fft,fftfreq,ifft
from scipy.signal import butter, lfilter
import matplotlib.pyplot as plt
import numpy as np
from numpy import array, sinc,cos,pi,abs,power,sin
from matplotlib.pyplot import plot, subplot,show
from matplotlib.pyplot import ylabel, xlabel
fs = 1000      # sampling rate
ts = 1.0/fs     # sampling interval
duration=5
t = np.arange(-duration,duration,ts)
N= 2*fs*duration
def ngovao():
    n = int(input("Số tín hiệu cần ghép kênh là: "))
    Ngovao=[]
    Tanso=[]
    Congsuat=[]
    Dang=[]
    for i in range(0,n):
        print("Ngõ vào thứ ",i+1,'dạng: ',end='')
        dang=input()
        A = float(input("Biên độ: "))
        f = float(input("Tần số: ")) #KHz
        if dang =='sinc':
            x= A*sinc(f*t)
            P = sum((abs((fft(x)))**2)/len(x))
        elif dang =='sinc binh' or dang=='sincbinh' or dang=='sinc^2':
            x= A*(sinc(f*t)**2)
            P = sum((abs((fft(x)))**2)/len(x))
        elif dang =='sin':
            x=A*sin(2*pi*f*t)
            P=(A**2)/2
        elif dang == 'cos':
            x=A*cos(2*pi*f*t)
            P=(A**2)/2
        Ngovao.append(x)
        Tanso.append(f)
        Congsuat.append(P)
        Dang.append(dang)
    fc=[2*Tanso[0]]
    for i in range(1,n):
        fc.append(fc[i-1]+2*Tanso[i]+Tanso[i-1]+5)
    for i in fc:
        a=[]
    if i<50:
        fcx=50
    else:
        a.append(i)
        fcx=min(a)
    if fc[n-1]+Tanso[n-1] < 500: 
        if fc[n-1]+Tanso[n-1]+fcx < 500:
            return(Ngovao,Tanso,Congsuat,n,fc,fcx)
    else:
        print("Tín hiệu của bạn không phù hợp với kênh!")
        print("Hãy nhập lại tín hiệu khác!")
        return(ngovao())
def dieuche(x,fc,P,kieu):
    Ac=4
    songmang=Ac*cos(2*pi*fc*t)
    xc=[None]*len(songmang)
    if kieu=='AM' or kieu=='am':
        for i in range(0,len(xc)):
            xc[i]=songmang[i]*(1+u*x[i])
        Sd=((Ac**2)*(1+(u**2)*P))/2
    elif kieu=='DSB' or kieu=='dsb':
        for i in range(0,len(xc)):
            xc[i]=songmang[i]*(x[i])
        Sd=((Ac**2)*P)/2
    return(xc,Sd)
def giaidieuche(xc,fc):
    Ac=4
    songmang=Ac*cos(2*pi*fc*t)
    y=xc*songmang
    return(y)
def biendoiFourier(x,n):
    X= (fft(x))/len(x)
    f= fftfreq(n,1/fs) #Tính tần số ở tâm của mỗi khối fft
    return(X,f)
def BPF(x,fc,f):
    nyq=fs/2
    fl=(fc-f)/nyq
    fh=(fc+f)/nyq
    b,a=butter(3,[fl, fh],'band')
    y=lfilter(b,a,x)
    return(y)
def LPF(x,f):
    nyq=fs/2
    fcut=f/nyq
    b,a=butter(3,fcut,'low')
    y=lfilter(b,a,x)
    return(y)
def noise():
    #giả sử kênh truyền có nhiễu ngẫu nhiên với N0=1
    nhieu=np.random.normal(0,1,len(t))
    Noi=fft(nhieu)
    Nd=0
    for i in Noi:
        Nd+=(abs(i))**2/len(Noi)
    return(nhieu,Nd)
def SNR(Sd,Nd):
    return(10 * np.log10(Sd/Nd))
#Nhập các tín hiệu cần ghép kênh
x,f,P,n,fc,fcx=ngovao()
plt.figure(1)
plt.suptitle("Dạng sóng tín hiệu ban đầu khi chưa ghép kênh",fontsize=15)
for i in range(0,n):
    subplot(len(x),1,i+1)
    plot(t,x[i])
plt.figure(4)
plt.suptitle("Phổ của tín hiệu ban đầu khi chưa ghép kênh",fontsize=15)
for i in range(0,len(x)):
    Mc,F=biendoiFourier(x[i],N)
    subplot(n,1,i+1)
    plot(F,abs(Mc))
'''xp=[None]*n
for i in range(0,n):
    xp[i]=LPF(x[i],f[i])
for i in range(0,n):
    subplot(len(x),1,i+1)
    plot(t,xp[i])'''
kieu=input("Chọn kiểu điều chế: ")
if kieu =="am" or kieu=='AM':
    u = float(input("Hệ số điều chế u="))
mc=[None]*len(x)
Pc=[None]*len(x)
xb=[]
#Điều chế và ghép kênh
plt.figure(2)
plt.suptitle("Phổ của từng kênh tín hiệu sau khi điều chế",fontsize=15)
for i in range(0,len(x)):
    mc[i],Pc[i]=dieuche(x[i],fc[i],P[i],kieu)
    Mc,F=biendoiFourier(mc[i],N)
    subplot(n,1,i+1)
    plot(F,abs(Mc))
Sb=sum(Pc)
for i in range(0,len(mc[0])):
    a=0
    for j in range(0,len(mc)):
        a+=mc[j][i]
    xb.append(a)
Xb,Fb=biendoiFourier(xb,N)
nhieu,Nd=noise()
#xb+=nhieu
#Điều chế chính
xc,Sd=dieuche(xb,fcx,Sb,kieu)
plt.figure(3)
plt.suptitle("Dạng sóng và phổ của tín hiệu sau khi ghép kênh",fontsize=15)
subplot(211)
plot(t,xc)
xlabel('t')
ylabel('xc(t)')
Xc,Fc=biendoiFourier(xc,N)
subplot(212)
plot(Fc,abs(Xc))
xlabel('f')
ylabel('Xc(f)')
#plt.figure(4)
#Phân kênh
y=giaidieuche(xc,fcx)
'''Y,Fy=biendoiFourier(y,N)
plt.suptitle("Phổ sau khi giải điều chế chính")
plot(Fy,abs(Y))'''
plt.figure(5)
s=[]
for i in range(0,n):
    yl=BPF(y,fc[i],f[i])
    yd=giaidieuche(yl,fc[i])
    xd=LPF(yd,f[i])
    s.append(xd)
plt.suptitle('Tín hiệu sau khi tách kênh của tín hiệu lần lượt là',fontsize=15)
for i in range(0,len(s)):
    subplot(len(s),1,i+1)
    plot(t,s[i])
plt.figure(6)
plt.suptitle("Phổ của tín hiệu sau khi tách kênh",fontsize=15)
for i in range(0,len(s)):
    S,F=biendoiFourier(s[i],N)
    subplot(n,1,i+1)
    plot(F,abs(S))
SNR=SNR(Sd,Nd)
print("SNR=",round(SNR,3),'dB')
print("Tần số điều chế chính fcx=",fcx)
show()