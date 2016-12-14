#!/usr/bin/env python
import numpy as np
import re
import cmath
import matplotlib.pyplot as plt
import sys
def load(filename,comment='#',delim='\t ',commentInv='false'):
    res=[]
    dds=re.compile('[%s]+'%(delim))
    cmline=re.compile('^[%s]+'%(comment))
    f=open(filename)
    cInv=False
    if commentInv.lower()=='true':
        cInv=True
    for line in f:
        line=line.strip()
        if line=='':
            continue
        cml=cmline.match(line)
        if not cml is None:
            if cInv:
                line=line[len(cml.group()):].strip()
            else:
                continue
        else:
            if cInv:
                continue
        tmp=dds.split(line)
        res.append(tmp)
    f.close()
    return res

def loadnp(filename,comment='#',delim='\t ',labelmark='!',commentInv='false'):
    tmp=load(filename,comment,delim,commentInv)
    Xa=[]
    Yas=[]
    #if the firstline is label:
    colname=[]
    if tmp[0][0]==labelmark:
        col=len(tmp[0][2:])
        colname=tmp[0][1:]
        tmp=tmp[1:]
    elif tmp[0][0].find(labelmark)==0:
        col=len(tmp[0][1:])
        tmp[0][0]=tmp[0][0][1:]
        colname=tmp[0]
        tmp=tmp[1:]
    else:
        col=len(tmp[0][1:])
        colname.append('x')
        for i in range(0,col):
            colname.append('y%d'%(i+1))
    #print col
    i=0
    for i in range(0,col):
        Yas.append([[],colname[i+1]])
    for line in tmp:
        Xa.append(line[0])
        for i in range(0,col):
            Yas[i][0].append(line[1:][i])

    res=[[np.array(Xa),colname[0]]]
    for Ya in Yas:
        res.append([np.array(Ya[0]),Ya[1]])
    return res

def str_to_strarr(string): #split by ','
    string = string.strip()
    return string.split(',')

def str_to_strnumarr(string):
    astr=re.compile('^[0-9]+(\.[0-9]+)?(\,[0-9]+(\.[0-9]+)?)*$') #matches Demicals or Integers
    astr2=re.compile('^([0-9]+(\.[0-9]+)?(\,[0-9]+(\.[0-9]+)?)*)$')
#    split=re.compile('(?<=[0-9]+(\.[0-9]+)?),(?=[0-9]+(\.[0-9]+)?)') 
    if not astr.match(string) is None:
        #spl=split.split(string)
        spl=string.split(',')
        return spl
    elif not astr2.match(string) is None:
        #spl=split.split(string)
        spl=string[1:-1].split(',')
        return spl
    else:
        return None

def str_to_floatarr(string):
    arr=str_to_strnumarr(string)
    res=[]
    for line in arr:
        res.append(float(line))
    return res

def str_to_intarr(string):
    arr=str_to_strnumarr(string)
    res=[]
    for line in arr:
        res.append(int(line))
    return res

#print loadnp('wv05.txt')
def doplot(plots,option):
    ars = loadnp(plots,commentInv=option['commentInv'])
    plt.figure(figsize=(4,3),dpi=80)
    plt.subplot(111)
    if option['xlabel'] is None:
        option['xlabel']=ars[0][1]
    plt.xlabel(option['xlabel'],fontsize=15)
    plt.ylabel(option['ylabel'],fontsize=15)
    if not option['xlim'] is None:
        val=str_to_floatarr(option['xlim'])
        plt.xlim(val[0],val[1])
    if not option['ylim'] is None:
        val=str_to_floatarr(option['ylim'])
        plt.ylim(val[0],val[1])
    plt.grid(which='major',linestyle='-')
    plt.title(option['title'])
    yplot=None
    if not option['yplot'] is None:
        yplot=str_to_strarr(option['yplot'])
    nowplot=0
    if not option['xticks'] is None:
        val = str_to_floatarr(option['xticks'])
        plt.xticks(np.arange(val[0],val[1],val[2]))
    if not option['yticks'] is None:
        val = str_to_floatarr(option['yticks'])
        plt.yticks(np.arange(val[0],val[1],val[2]))
    if not option['xscale'] is None:
        plt.xscale(option['xscale'])
    if not option['yscale'] is None:
        plt.yscale(option['yscale'])
    """if not option['xticks_minor'] is None:
        val = str_to_floatarr(option['xticks_minor'])
        plt.xticks(np.arange(val[0],val[1],val[2]),minor=True)
    if not option['yticks_minor'] is None:
        val = str_to_floatarr(option['yticks_minor'])
        plt.yticks(np.arange(val[0],val[1],val[2]),minor=True)"""
    #plt.grid(which='both')
    #plt.grid(which='minor',alpha=0.2)
    #plt.grid(which='major',alpha=0.5)
    for y in ars[1:]:
        nowplot+=1
        if not yplot is None:
            if str(nowplot) not in yplot and y[1] not in yplot :
                continue
        plt.plot(ars[0][0],y[0],linewidth=1.0,linestyle='-',label=y[1])
    plt.legend()
    if not option['save'] is None:
        plt.savefig(option['save'])
    else:
        plt.show()

#x=0.5 , 1.0 , 3.0 , 5.0
def main(argv):
    option={}
    #option.update({'subplot':111})
    option.update({'xlabel':None})
    option.update({'ylabel':'y'})
    option.update({'xlim':None})
    option.update({'ylim':None})
    option.update({'title':None})
    option.update({'yplot':None})
    option.update({'xticks':None})
    option.update({'yticks':None})
    option.update({'xscale':None})
    option.update({'yscale':None})
    option.update({'save':None})
    option.update({'commentInv':'false'})
    #option.update({'xticks_minor':None}) #not available
    #option.update({'yticks_minor':None})
    if len(argv)>1:
        plots=''
        is_readfile=False
        for read in argv[1:]:
            check1=re.compile('[^=]+$')
            check2=re.compile('[^=]+=[^=]+')
            if not check1.match(read) is None: #filename
                if not is_readfile:
                    plots=read
                    is_readfile=True
                    if option['title']==None:
                        option['title']=read
                else:
                    print 'Filename multityped!'
                    return
            elif not check2.match(read) is None:
                spl=read.split('=')
                if spl[0] in option.keys():
                    option[spl[0]]=spl[1]
                else:
                    print 'can\'t recognize OPTION "%s"'%(spl[0])
                    return
        doplot(plots,option)
    else:
        print 'usage:plotfile.py [FILENAME] [OPTION=VALUE]'
        print 'OPTIONs now available:'
        for key in option.keys():
            print key

if __name__=='__main__':
    main(sys.argv)
