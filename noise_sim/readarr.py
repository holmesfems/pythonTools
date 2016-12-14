#!/usr/bin/env python
import numpy as np
import re
import sys
def loadarr(filename,comment='#',delim='\t ',commentInv='false',form=float):
    res=[]
    dds=re.compile('[%s]+'%(delim))
    cmline=re.compile('^[%s]+'%(comment))
    f=open(filename)
    cInv=False
    first=True
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
        if first:
            for i in range(0,len(tmp)):
                res.append([])
            first=False
        for i in range(0,len(tmp)):
            res[i].append(form(tmp[i]))
    f.close()
    return res
