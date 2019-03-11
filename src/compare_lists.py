#!/usr/bin/python
import sys

infile1=open(sys.argv[1],'r')
infile2=open(sys.argv[2],'r')
outfile1=open("common.txt",'w')
outfile2=open("diff.txt",'w')

list1=[]
list2=[]
for line in infile1:
	line=line.replace("\n","")
	list1.append(line)
for line in infile2:
	line=line.replace("\n","")
	list2.append(line)

print len(list1)
print len(list2)

for i in list1:
	enter=False
	for j in list2:
		if i==j:
			enter=True
	if enter==True:
		outfile1.write(i)
		outfile1.write("\n")
	else:
		outfile2.write(i)
		outfile2.write("\n")
	
	
