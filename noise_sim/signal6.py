import math
import random
import scipy.fftpack as fftp
import numpy as np
import cmath
import sys
import re
import newton
from scipy.optimize import fmin


def gaussianRnd(sig2=1.0): #function that get gaussian random number
    x1=0.0
    x2=0.0
    while(x1==0.0)and(x2==0.0):
        x1=random.random()
        x2=random.random()
    y1=math.sqrt(-2*sig2*math.log(x1))*math.cos(2*math.pi*x2)
    y2=math.sqrt(-2*sig2*math.log(x1))*math.sin(2*math.pi*x2)
    return [y1,y2]
def gaussian(x,mu=0.0,sig2=1.0,amp=1.0): #Gaussian distribution
    return amp*math.exp(-(x-mu)**2/(2.0*sig2))

#get simulation of noise of atmosphere
def atmosNoise(points=30000,max_amp=30000.0,white_amp=1000.0):
    F=[]
    F.append(0.0)
    for i in range(1,points/2):
        f=max_amp/(i)**0.5*(gaussianRnd()[0]+1j*(gaussianRnd()[0]))+white_amp*(gaussianRnd()[0]+1j*gaussianRnd()[0])
        F.append(f)
    F.append(max_amp/(points/2)**0.5*(gaussianRnd()[0])+white_amp*(gaussianRnd()[0]))
    for i in range(points/2+1,points):
        f=F[points-i].real-1.0j*F[points-i].imag
        F.append(f)
    Y=fftp.ifft(F)
    return [F,Y]

def signaldata(center=15.0,amp=10.0,sig2=20,time=0):
    return gaussian(x=time,mu=center,amp=amp,sig2=sig2)

def signal(points=30000,center=15.0,amp=10.0,sig2=20,time=None):
    Y=[]
    for i in range(0,points):
        Y.append(signaldata(center=center,amp=amp,sig2=sig2,time=time[i]))
    return Y

def switch(points=30000,freq=67,phase=0):
    Y=[]
    for i in range(0,points):
        if ((i+phase)/freq)%2==0:
            Y.append(0)
        else:
            Y.append(1)
    return Y

def constAr(points,const):
    Y=[]
    for i in range(0,points):
        Y.append(const)
    return Y

def zeroAr(points=30000):
    return constAr(points,0)

def mask_test(switch,mask_back=5,mask_go=5):
    Y=zeroAr(len(switch))
    s0=switch[0]
    for i in range(1,len(switch)):
        if s0!= switch[i]:
            for j in range(i-mask_back,i+mask_go-1):
                if j<0 or j>=len(switch):
                    continue
                Y[j]=1
        s0=switch[i]
    return Y

def time(points=30000,maxtime=30.0):
    Y=[]
    for i in range(0,points):
        Y.append(maxtime/points*i)
    return Y

def arReal(cArray):
    newA=[]
    for term in cArray:
        newA.append(term.real)
    return newA

def arImag(cArray):
    newA=[]
    for term in cArray:
        newA.append(term.imag)
    return newA

def combineAr(ar1,ar2):
    ret=[]
    for i in range(0,min(len(ar1),len(ar2))):
        ret.append([ar1[i],ar2[i]])
    return ret

def subtAr(ar1,ar2):
    ret=[]
    for i in range(0,len(ar1)):
        ret.append(ar1[i]-ar2[i])
    return ret

def copyAr(ar):
    ret=[]
    for line in ar:
        ret.append(line)
    return ret

def devideAr(ar):
    ret=[]
    for i in range(0,len(ar[0])):
        ret.append([])
    for item in ar:
        for i in range(0,len(item)):
            ret[i].append(item[i])
    return ret

def absAr(ar):
    ret=[]
    for item in ar:
        ret.append(abs(item))
    return ret

def saveAr(XYarray,filename,yname='y',xname='x'):
    f=open(filename,'w')
    f.write('!%s %s\n'%(xname,yname))
    #print filename
    for line in XYarray:
        f.write('%f\t%f\n'%(line[0],line[1]))
    f.close()

def moveAverage(array,average_halfLenth,weight=None):
    if weight==None:
        weight=constAr(2*average_halfLenth+1,1)
    Y=[]
    for i in range(0,len(array)):
        s=0
        w=0
        for j in range(i-average_halfLenth,i+average_halfLenth+1):
            if j<0 or j>=len(array):
                s+=0
            else:
                s+=array[j]*weight[j-i+average_halfLenth]
            w+=weight[j-i+average_halfLenth]
        Y.append(s/w)
    return Y

def gaussianFit(array):
    s=0
    mu=0
    sig2=0
    dx=(array[len(array)-1][0]-array[0][0])/(len(array)-1)
    for i in range(0,len(array)):
        y=array[i][1]
        x=array[i][0]
        s+=y
        mu+=y*x
        sig2+=y*x*x
    mu=mu/s
    sig2=sig2/s-mu*mu
    if sig2<0:
        return [None,mu,sig2,0]
    amp=s*dx/(2*math.pi*sig2)**0.5
    print [mu,sig2,amp]
    Y=[]
    for i in range(0,len(array)):
        Y.append(gaussian(array[i][0],mu=mu,sig2=sig2,amp=amp))
    return [Y,mu,sig2,amp]

def gaussianFit2(array):
    mu=15
    sig2=10
    amp=10
    A0=[mu,sig2,amp]
    def S(A):
        s=0
        for item in array:
            s+=(item[1]-gaussian(item[0],mu=A[0],sig2=A[1],amp=A[2]))**2
        return s
    A1=newton.solve_minimize(S,A0,0.0001,0.0001)
    print A1
    if A1==None:
        return
    Y=[]
    for i in range(0,len(array)):
        Y.append(gaussian(array[i][0],mu=A1[0],sig2=A1[1],amp=A1[2]))
    return [Y,A1[0],A1[1],A1[2]]

def gaussianFit3(array):
    mu=15
    sig2=10
    amp=10
    A0=[mu,sig2,amp]
    def S(A):
        s=0
        for item in array:
            s+=(item[1]-gaussian(item[0],mu=A[0],sig2=A[1],amp=A[2]))**2
        return s
    A1=fmin(S,A0)
    print A1
    if A1==None:
        return
    Y=[]
    for i in range(0,len(array)):
        Y.append(gaussian(array[i][0],mu=A1[0],sig2=A1[1],amp=A1[2]))
    return [Y,A1[0],A1[1],A1[2]]

def doEmu(OPTION):
    t=time()
    sig=signal(amp=float(OPTION['sig_amp']),sig2=float(OPTION['sig_sig'])**2,time=t)
    swi=switch()
    mask=mask_test(swi)
    noise=atmosNoise()
    output=[]
    sigswi=[]
    output2=[]
    output3=[]
    for i in range(0,30000):
        sigswi.append(sig[i]*swi[i])
        output.append(sig[i]*swi[i]+noise[1][i].real)
        output2.append(sig[i]+noise[1][i].real)
        output3.append(output[i]*(1-mask[i]))
    sigswi_spec=fftp.fft(sigswi)
    swi_spec=fftp.fft(swi)
    DOPLOT={
            "number":range(0,30000),
            "signal":sig,
            "switch":swi,
            "time":t,
            "noise":noise[1],
            "noise_real":arReal(noise[1]),
            "noise_spec":noise[0],
            "output":output,
            "output2":output2,
            "output3":output3,
            "mask":mask,
            "sigswitch_spec_real":arReal(sigswi_spec),
            "sigswitch_spec_imag":arImag(sigswi_spec),
            "sigswitch_spec":sigswi_spec,
            "switch_spec":swi_spec
            }
    remove=[
        'noise','sigswitch_spec','switch_spec','noise_spec'
    ]
    saveAr(combineAr(DOPLOT['number'],absAr(noise[0][:15000])),'noise_spec_abs.list','noise_spec_abs','x')
    for key in DOPLOT.keys():
        if key not in remove:
            saveAr(combineAr(DOPLOT['time'],DOPLOT[key]),'%s.list'%key,key,'time')
    """
    print '!',
    for key in OPTION['plot'].split(','):
        print key,
    print
    for i in range(0,30000):
        for key in OPTION['plot'].split(','):
            print DOPLOT[key][i],
        print
    """
    return DOPLOT

def doAve(DOPLOT,OPTION):
    #mode = Normal Average:
    count_Noise=0
    sum_Noise=0
    count_Sig=0
    sum_Sig=0
    res=[]
    swi_cond=0
    first=True
    starttime=0
    for i in range(0,30000):
        if DOPLOT['switch'][i] > 0:
            if swi_cond==0:
                if first:
                    first=False
                else:
                    if count_Sig>10:
                        res.append([(starttime+DOPLOT['time'][i-1])/2.0,sum_Sig/count_Sig-sum_Noise/count_Noise])
                starttime=DOPLOT['time'][i]
                count_Sig=0
                swi_cond=1
                count_Noise=0
                sum_Noise=0
                sum_Sig=0
            if DOPLOT['mask']==1:
                continue
            sum_Sig+=DOPLOT['output'][i]
            count_Sig+=1
        else:
            if first:
                continue
            if DOPLOT['mask']!=1:
                sum_Noise+=DOPLOT['output'][i]
                count_Noise+=1
                swi_cond=0
            if i==30000-1 and count_Noise>10:
                res.append([(starttime+DOPLOT['time'][i])/2.0,sum_Sig/count_Sig-sum_Noise/count_Noise])
    sampling_time=devideAr(res)[0]
    """
    for item in sampling_time:
        print item
    """
    sampling_data=[]
    for t in sampling_time:
        sampling_data.append(signaldata(amp=float(OPTION['sig_amp']),sig2=float(OPTION['sig_sig'])**2,time=t))
    saveAr(combineAr(sampling_time,sampling_data),'sampling_data.list','sampling_data')
    saveAr(res,'average.list','average')
    saveAr(combineAr(sampling_time,subtAr(devideAr(res)[1],sampling_data)),'average-signal.list','average-signal')
    gf=0
    if OPTION['fit_mode']=='2':
        gf=gaussianFit2(res)
    elif OPTION['fit_mode']=='1':
        gf=gaussianFit(res)
    elif OPTION['fit_mode']=='3':
        gf=gaussianFit3(res)
    if not gf[2]<0:
        saveAr(combineAr(sampling_time,gf[0]),'average_fit.list','average_fit')
    #saveAr(combineAr(sampling_time,absAr(devideAr(res)[1])),'average_abs.list','average_abs')
    return gf
    """
    res2=moveAverage(devideAr(res)[1],5)
    saveAr(combineAr(sampling_time,res2),'average_moveAve.list','average_moveAve')
    res3=moveAverage(res2,5)
    saveAr(combineAr(sampling_time,res3),'average_moveAve2.list','average_moveAve2')
    """

def devideCmd(string):
    return string.split('=')

argv=sys.argv
OPTION={
        'plot':'time,output',
        'sig_amp':'50',
        'sig_sig':'3.5',
        'fit_mode':'2',
        'roop_count':'1',
        'write_mode':'w'
}
match1=''
for key in OPTION.keys():
    match1=match1+key+'|'
match1=re.compile('(%s)=[^=]+'%match1[:-1])

for cmd in argv[1:]:
    if not ( match1.match(cmd) is None):
        tmp=devideCmd(cmd)
        OPTION[tmp[0]]=tmp[1]
f=open('mu-sig2.list',OPTION['write_mode'])
rc=int(OPTION['roop_count'])
for i in range(0,rc):
    gf=doAve(doEmu(OPTION),OPTION)
    print i,
    f.write('%f\t%f\n'%(gf[1],gf[2]))
    if gf[2]<0:
        print 'sig2<0',
    print
f.close()
