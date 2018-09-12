import sys
import os
import re
sys.path.append('/home/ennisb/me/explore-cgi/src/rules/')
import const

infile = sys.argv[1]
output = sys.argv[2]
GEMINI_DIR = const.GEMINI_DIR

output = GEMINI_DIR + output

with open(infile, 'r') as inf, open(output, 'w') as of: 
    for line in inf:
        new_line = re.split('\t |:|-', line)
        new_line = str('\t'.join(new_line))
        of.write(new_line)


