#!/usr/bin/ipython3
import re

with open('./parsed_lib.yml','w') as output:
	lib =  open(r'./library.yml','r')
	tags = []
	temp_div = None
	dividers = []
	for line in lib:
		outline = ''
		if 'DIVIDER' in line:
			temp_div = re.search('(?<=:).*$',line).group(0)[1:]
			continue
		if re.match('\ *-',line) != None:
			outline = re.match('^\ *',line).group(0)
			g = re.search('(?<=-).*?(?=:)', line)
			if g != None:
				if temp_div != None:
					dividers.append([g.group(0),temp_div])
					temp_div = None
				if len(tags) >= 2:
					if g.group(0) == tags[-2]:
						if tags[-1] == dividers[-1][0]:
							dividers.pop()
						tags.pop()
						outline += ''*len(tags)
						outline += re.search('(?<=:).*$',line).group(0)[1:] + ':'
					else:
						if(tags[-1] != g.group(0)):
							tags.append(g.group(0))
						outline += ''*len(tags)
						outline += re.search('(?<=:).*$',line).group(0)[1:] + ':'
				else:
					tags.append(g.group(0))
					outline += ''*len(tags)
					outline += re.search('(?<=:).*$',line).group(0)[1:] + ':'
			outline += '\n'
			if len(dividers) != 0:
				if g.group(0) == dividers[-1][0]:
					outline += re.match('^\ *',line).group(0) + '  ' + 'division: ' + dividers[-1][1]+'\n'

		else:
			outline += ''*len(tags)
			outline += line
		print(line)
		#print(tags)
		#outline += '\n'
		output.write(outline)
