'''
Created on Jan 21, 2015

@author: s3cha
'''
import os
import sys

f = open(sys.argv[1],'r')
for line in f:
    exe_line = line
#     exe_line = f.readline()
    data = exe_line.split(' ')
    check = False
    print data
    for item in data:
        print item,check
        if check:
            if not os.path.exists(os.path.split(item)[0]):
                os.system('mkdir '+os.path.split(item)[0])
            break
        if item == '-o':
            check = True
    os.system(exe_line)
f.close()
# s = open(sys.argv[2]+'/'+sys.argv[1]+'_dummy.txt','w')
# s.close()
print 'END'
    