#!/usr/bin/python

import re
import os
import sys

batch = [""]

cl = re.compile('Parameters from commandline')
tl = re.compile('T = ')
dl = re.compile('==> Ground State Binding Energy : ')

for id in batch:
    output = []
    fm = re.compile("zrun\d+%s.log" % id)
    dirList = os.listdir('.')
    for fname in dirList:
        fd = False
        m1 = fm.match(fname)

	print("m1 = %s" % m1)

        if m1:
            file = open(fname,"r")
	    fp = False
            for line in file:
                m2 = cl.match(line)

		print("m2 = %s" % m2)

	        if m2:
		    fp = True
	        m3 = tl.match(line)

		print("m3 = %s" % m3)

	        if (m3 and fp):
	          (junk, t) = line.split('= ')
	          t = eval(t.strip())
	          print("t = %s" % t)

		m4 = dl.match(line)

		print("m4 %s" % m4)

                if m4:
                    (junk, data) = line.split(": ")
                    data = data.strip()
		    (rep,imp) = data.split(',')
		    rep = -eval(rep[1:])
		    imp = eval(imp[:-1])
                    fd = True
            if fd:
                output += [[t,rep,imp]]
            file.close()

    soutput = sorted(output)
    of = open("be%s.txt" % id,"w")	
    for item in soutput:
       of.write("%f\t%f\t%f\n" % (item[0],item[1],item[2]))
    of.close()
    
    #TEMPORARY CODE
    of2 = open("be2%s.txt" % id,"w")
    for item in output:
       of2.write("%f\t%f\t%f\n" % (item[0],item[1],item[2]))
    of2.close()

