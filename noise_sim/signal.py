import math
import random
import scipy.fftpack as fftp
import numpy as np
import cmath
import sys
import re
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

def atmosNoise(points=30000,max_amp=30000.0,white_amp=100.0):
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

def signal(points=30000,center=15000.0,amp=1000.0,sig2=1000000):
    Y=[]
    for i in range(0,points):
        Y.append(gaussian(x=i,mu=center,amp=amp,sig2=sig2))
    return Y

def switch(points=30000,freq=67,phase=0):
    Y=[]
    for i in range(0,points):
        if ((i+phase)/freq)%2==0:
            Y.append(0)
        else:
            Y.append(1)
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

def doEmu(OPTION):
    sig=signal(amp=float(OPTION['sig_amp']),sig2=float(OPTION['sig_sig'])**2)
    swi=switch()
    t=time()
    noise=atmosNoise()
    output=[]
    sigswi=[]
    for i in range(0,30000):
        sigswi.append(sig[i]*swi[i])
        output.append(sig[i]*swi[i]+noise[1][i].real)
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
            "sigswitch_spec_real":arReal(sigswi_spec),
            "sigswitch_spec_imag":arImag(sigswi_spec),
            "sigswitch_spec":sigswi_spec,
            "switch_spec":swi_spec
            }
    print '!',
    for key in OPTION['plot'].split(','):
        print key,
    print
    for i in range(0,30000):
        for key in OPTION['plot'].split(','):
            print DOPLOT[key][i],
        print
    return DOPLOT

def doAve(DOPLOT,OPTION):
    #mode = Normal Average:
    if OPTION['ave_mode']=='normal':
        count_Noiy30yse=0
        sum_Noise=0
        count_Sig=0
        sum_Sig=0
        res=[]
        swi_cond=DOPLOT['switch'][0]
        for i in range(0,30000):
            if DOPLOT['switch'][i] > 0:
                sum_Sig+=DOPLOT['output'][i]
                count_Sig+=1
                if swi_cond==0:
                    swi_cond=1
                    res.append([i-count_Noise/2.0,sum_Noise/count_Noise])
                    count_Noise=0
                    sum_Noise=0
                if i==30000-1:
                    res.append([i-count_Sig/2.0,sum_Sig/count_Sig])
            else:
                sum_Noise+=DOPLOT['output'][i]
                count_Noise+=1
                if swi_cond==1:
                    swi_cond=0
                    res.append([i-count_Sig/2.0,sum_Sig/count_Sig])
                    count_Sig=0
                    sum_Sig=0
                if i==30000-1:
                    res.append([i-count_Noise/2.0,sum_Noise/count_Noise])
        for line in res:
            print '#',line[0],line[1]
        return res
    elif OPTION['ave_mode']=='spec1':
        #create filter:
        filt=[]
        for term in DOPLOT['switch_spec']:
            if abs(term)>500:
                filt.append(1)
            else:
                filt.append(0)
        spec=fftp.fft(DOPLOT['output'])
        new_spec=[]
        for i in range(0,len(spec)):
            new_spec.append(spec[i]*filt[i])
        res=fftp.ifft(new_spec)
        res_r=[]
        for term in res:
            res_r.append(term.real)
        count_Noise=0
        sum_Noise=0
        flat=[]
        swi_cond=DOPLOT['switch'][0]
        for i in range(0,30000):
            if DOPLOT['switch'][i] > 0:
                if swi_cond==0:
                    swi_cond=1
                    flat.append([i-count_Noise/2.0,sum_Noise/count_Noise])
                    count_Noise=0
                    sum_Noise=0
            else:
                sum_Noise+=res_r[i]
                count_Noise+=1
                if swi_cond==1:
                    swi_cond=0
                if i==30000-1:
                    flat.append([i-count_Noise/2.0,sum_Noise/count_Noise])
        nowflat_pos=0
        new_res_r=[]
        for i in range(0,len(res_r)):
            if i <= flat[0][0]:
                new_res_r.append(res_r[i]-flat[0][1])
            elif i >= flat[len(flat)-1][0]:
                new_res_r.append(res_r[i]-flat[len(flat)-1][1])
            else:
                while flat[nowflat_pos+1][0]<=i:
                    nowflat_pos+=1
                nowflat=(flat[nowflat_pos][1]+flat[nowflat_pos+1][1])/2
                new_res_r.append(res_r[i]-nowflat)
        print '#','!','time','output'
        for i in range(0,len(new_res_r)):
            print '#',DOPLOT['time'][i],new_res_r[i]
        return res

def devideCmd(string):
    return string.split('=')

argv=sys.argv
OPTION={
        'plot':'time,output',
        'sig_amp':'50',
        'sig_sig':'4000',
        'ave_mode':'spec1'
}
match1=''
for key in OPTION.keys():
    match1=match1+key+'|'
match1=re.compile('(%s)=[^=]+'%match1[:-1])

for cmd in argv[1:]:
    if not ( match1.match(cmd) is None):
        tmp=devideCmd(cmd)
        OPTION[tmp[0]]=tmp[1]
doAve(doEmu(OPTION),OPTION)
