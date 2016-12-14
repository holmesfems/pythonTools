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
    for i in range(0,len(ar1)):
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

def saveAr(XYarray,filename,yname='y',xname='x'):
    f=open(filename,'w')
    f.write('!%s %s\n'%(xname,yname))
    print filename
    for line in XYarray:
        f.write('%f\t%f\n'%(line[0],line[1]))
    f.close()

def doEmu(OPTION):
    t=time()
    sig=signal(amp=float(OPTION['sig_amp']),sig2=float(OPTION['sig_sig'])**2,time=t)
    swi=switch()
    noise=atmosNoise()
    output=[]
    sigswi=[]
    output2=[]
    for i in range(0,30000):
        sigswi.append(sig[i]*swi[i])
        output.append(sig[i]*swi[i]+noise[1][i].real)
        output2.append(sig[i]+noise[1][i].real)
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
            "sigswitch_spec_real":arReal(sigswi_spec),
            "sigswitch_spec_imag":arImag(sigswi_spec),
            "sigswitch_spec":sigswi_spec,
            "switch_spec":swi_spec
            }
    remove=[
        'noise','sigswitch_spec','switch_spec','noise_spec'
    ]
    for key in DOPLOT.keys():
        if key not in remove:
            saveAr(combineAr(DOPLOT['time'],DOPLOT[key]),'%s.list'%key,key,'time')
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
                    res.append([(starttime+DOPLOT['time'][i-1])/2.0,sum_Sig/count_Sig-sum_Noise/count_Noise])
                starttime=DOPLOT['time'][i]
                count_Sig=0
                swi_cond=1
                count_Noise=0
                sum_Noise=0
                sum_Sig=0
            sum_Sig+=DOPLOT['output'][i]
            count_Sig+=1
        else:
            if first:
                continue
            sum_Noise+=DOPLOT['output'][i]
            count_Noise+=1
            swi_cond=0
            if i==30000-1:
                res.append([(starttime+DOPLOT['time'][i])/2.0,sum_Sig/count_Sig-sum_Noise/count_Noise])
    sampling_time=devideAr(res)[0]
    for item in sampling_time:
        print item
    sampling_data=[]
    for t in sampling_time:
        sampling_data.append(signaldata(amp=float(OPTION['sig_amp']),sig2=float(OPTION['sig_sig'])**2,time=t))
    saveAr(combineAr(sampling_time,sampling_data),'sampling_data.list','sampling_data')
    saveAr(res,'average.list','average')
    saveAr(combineAr(sampling_time,subtAr(devideAr(res)[1],sampling_data)),'average-signal.list','average-signal')
    #output2's average
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
                    res.append([(starttime+DOPLOT['time'][i-1])/2.0,(sum_Sig+sum_Noise)/(count_Sig+count_Noise)])
                starttime=DOPLOT['time'][i]
                count_Sig=0
                swi_cond=1
                count_Noise=0
                sum_Noise=0
                sum_Sig=0
            sum_Sig+=DOPLOT['output2'][i]
            count_Sig+=1
        else:
            if first:
                continue
            sum_Noise+=DOPLOT['output2'][i]
            count_Noise+=1
            swi_cond=0
            if i==30000-1:
                res.append([(starttime+DOPLOT['time'][i])/2.0,(sum_Sig+sum_Noise)/(count_Sig+count_Noise)])
    saveAr(res,'average2.list','average2')
    saveAr(combineAr(sampling_time,subtAr(devideAr(res)[1],sampling_data)),'average2-signal.list','average2-signal')
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
    saveAr(combineAr(DOPLOT['time'],new_res_r),'specfilter.list','specfilter')
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
                    res.append([(starttime+DOPLOT['time'][i-1])/2.0,sum_Sig/count_Sig-sum_Noise/count_Noise])
                starttime=DOPLOT['time'][i]
                count_Sig=0
                swi_cond=1
                count_Noise=0
                sum_Noise=0
                sum_Sig=0
            sum_Sig+=new_res_r[i]
            count_Sig+=1
        else:
            if first:
                continue
            sum_Noise+=new_res_r[i]
            count_Noise+=1
            swi_cond=0
            if i==30000-1:
                res.append([(starttime+DOPLOT['time'][i])/2.0,sum_Sig/count_Sig-sum_Noise/count_Noise])
    saveAr(res,'specfilter_sample.list','specfilter_sample')
    saveAr(combineAr(sampling_time,subtAr(devideAr(res)[1],sampling_data)),'specfilter-signal.list','specfilter-signal')

def devideCmd(string):
    return string.split('=')

argv=sys.argv
OPTION={
        'plot':'time,output',
        'sig_amp':'50',
        'sig_sig':'3.5',
        'ave_mode':'normal'
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
