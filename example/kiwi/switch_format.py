import datetime
import numpy as np
fname="KIW06-Mid_6m.ios"
outname="kiw06.dat"
f = open(fname,'r')
g = open(outname,'w')
lines = f.readlines()
rnum = np.random.randn(len(lines))*4.

for i,line in enumerate(lines):
    num,dat,tm,u,v = line.strip().split()
    date_string = "%s %s" % (dat,tm)
    dt =datetime.datetime.strptime(date_string, "%Y/%m/%d %H:%M:%S.00")
    outstring=dt.strftime("%d%m%Y%H%M")
    outstring=dt.strftime("%Y-%m-%d %H:%M")
    u = float(u)
    v = float(v)
    up = u + 0.1
    vp = v + .2
    #if u > -990 and v > -990:
    #g.write("%s %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n" % (outstring,u,v,up,vp,u,v))
    g.write("%s %9.4f %9.4f\n" % (outstring,u/10.,v/10.))
f.close()
g.close()

