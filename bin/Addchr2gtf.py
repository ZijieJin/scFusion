import sys

infile = open(sys.argv[1])

for line in infile.readlines():
	if line.startswith('#'):
		print(line, end='')
		continue
	if not line.startswith('chr'):
		print('chr' + line, end='')
	else:
		print(line, end='')

infile.close()
