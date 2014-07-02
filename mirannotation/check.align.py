#!/usr/bin/python

#load simulated reads
sim=open('sim.20.hsa.fa','r')
data={}
dataseq={}
lendata={}
name=""
for line in sim:
	line=line.strip()
	if (line.find(">")>=0): 
		name=line.replace(">","")
		slot=name.split("_")[0].split("-")
		data[name]="-".join(slot[0:3])
		#print name
	else:
		lendata[name]=len(line)
		dataseq[line]=name
sim.close()

#load miraligner results
check={}
mir=open('miraligner/sim.20.hsa.mirna')
for line in mir:
	cols=line.split("\t")
	slot=cols[2].split("-")
	add=line.find("add")
	mut=line.find("mut")
	if (line.find("hsa")>=0):
		if (line.find("miRNA")>=0) and (not check.has_key(cols[1])):
			check[cols[1]]=1
			print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner" %(cols[1],"-".join(slot[0:3]),data[cols[1]],add,mut,lendata[cols[1]])
		elif (line.find("precursor")>=0) and (not check.has_key(cols[1])):
			check[cols[1]]=1
			slot=cols[2].split("-")
			print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tmiraligner" %(cols[1],"-".join(slot[0:3]),data[cols[1]],add,mut,lendata[cols[1]])
mir.close()

#if not aligned, print them here
for k in data.keys():
	if not check.has_key(k):
		add=k.find("add")
		mut=k.find("mut")
		print "%s\t%s\tno\tNA\t%s\t%s\t%s\tmiraligner" % (k,data[k],add,mut,lendata[k])

#load bwotie2 results
check={}
mir=open('bowtie2/sim.20.hsa.sam')
for line in mir:
	cols=line.split("\t")
	if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
		check[cols[0]]=1
		slot=cols[2].split("-")
		add=cols[0].find("add")
		mut=cols[0].find("mut")
		print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tbowtie2" %(cols[0],"-".join(slot[0:3]),data[cols[0]],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
	if not check.has_key(k):
		add=k.find("add")
		mut=k.find("mut")
		print "%s\t%s\tno\tNA\t%s\t%s\t%s\tbowtie2" % (k,data[k],add,mut,lendata[k])

#load novoaligner results
check={}
mir=open('novo/sim.20.novo.sam')
for line in mir:
	cols=line.split("\t")
	if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
		slot=cols[2].split("-")
		add=cols[0].find("add")
		mut=cols[0].find("mut")
		check[cols[0]]=1
		print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tnovocraft" %(cols[0],"-".join(slot[0:3]),data[cols[0]],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
	if not check.has_key(k):
		add=k.find("add")
		mut=k.find("mut")
		print "%s\t%s\tno\tNA\t%s\t%s\t%s\tnovocraft" % (k,data[k],add,mut,lendata[k])


#load srnabench results
check={}
mir=open('srnabench/reads_orig.fa')
srna={}
name=""
for line in mir:
	line=line.strip()
	if (line.find(">")>=0): 
		name=line.replace(">","")
	else:
		srna[name]=line
mir.close()

mir=open('srnabench/hairpin.parsed')
for line in mir:
	cols=line.split("\t")
	seq=srna[cols[0]]
	if (cols[2].find("hsa")>=0) and (not check.has_key(dataseq[seq])):
		slot=cols[2].split("-")
		add=dataseq[seq].find("add")
		mut=dataseq[seq].find("mut")
		check[dataseq[seq]]=1
		print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tsrnabench" %(dataseq[seq],"-".join(slot[0:3]),data[dataseq[seq]],add,mut,lendata[dataseq[seq]])
mir.close()

#if not aligned, print them here
for k in data.keys():
	if not check.has_key(k):
		add=k.find("add")
		mut=k.find("mut")
		print "%s\t%s\tno\tNA\t%s\t%s\t%s\tsrnabench" % (k,data[k],add,mut,lendata[k])


#load GEM results
check={}
mir=open('gem/sim.20.hsa.sam')
for line in mir:
	if line.find("@")<0:
		cols=line.split("\t")
		if (cols[2].find("hsa")>=0) and (not check.has_key(cols[0])):
			slot=cols[2].split("-")
			add=cols[0].find("add")
			mut=cols[0].find("mut")
			check[cols[0]]=1
			print "%s\t%s\tyes\t%s\t%s\t%s\t%s\tGEM" %(cols[0],"-".join(slot[0:3]),data[cols[0]],add,mut,lendata[cols[0]])
mir.close()

#if not aligned, print them here
for k in data.keys():
	if not check.has_key(k):
		add=k.find("add")
		mut=k.find("mut")
		print "%s\t%s\tno\tNA\t%s\t%s\t%s\tGEM" % (k,data[k],add,mut,lendata[k])
