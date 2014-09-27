# runModel.py

import sys
import subprocess
import os
import shutil

paramFile = open(sys.argv[1], 'r')
for params in paramFile:
	if not params.startswith("#"):
		data = params.split(" ")
		k = int(data[0].strip())
		m = int(data[1].strip())
		n = int(data[2].strip())
		p1 = float(data[3].strip())
		p2 = float(data[4].strip())
		cmdline = str(k) + " " + str(m) + " " + str(n) + " " + str(p1) + " " + str(p2)
		p = subprocess.Popen('java Model ' + cmdline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		seen = []
		data = p.stdout.readlines()
		for time in data: # Should only be a single line
			time = time.strip()
			print >> sys.stderr, str(k) + "," + str(m) + "," + str(n) + "," + str(p1) + "," + str(p2) + "," + time
			print(str(k) + "," + str(m) + "," + str(n) + "," + str(p1) + "," + str(p2) + "," + time)
