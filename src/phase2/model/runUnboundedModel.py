# runUnboundModel.py

import datetime
import sys
import subprocess
import os
import shutil

k = sys.argv[1]
m = sys.argv[2]
p1 = sys.argv[3]
p2 = sys.argv[4]
n = 2

while True:
	cmdline = str(k) + " " + str(m) + " " + str(n) + " " + str(p1) + " " + str(p2)
	start = datetime.datetime.now()
	p = subprocess.Popen('java Model ' + cmdline, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	seen = []
	data = p.stdout.readlines()
	for time in data: # Should only be a single line
		time = time.strip()
		print >> sys.stderr, str(k) + "," + str(m) + "," + str(n) + "," + str(p1) + "," + str(p2) + "," + time
		print(str(k) + "," + str(m) + "," + str(n) + "," + str(p1) + "," + str(p2) + "," + time)
		end = datetime.datetime.now()
		print >> sys.stderr, str(end - start)
		print str(end - start)

	# Next iteration
	n = n + 1