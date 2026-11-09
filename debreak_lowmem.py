"""
Low-memory implementations of the memory-heavy DeBreak steps.

Nothing in the original modules is modified. The driver (`debreak`) uses the
functions below only when --lowmem is given; without the flag the original code
path runs unchanged. Every function here writes the same output files as its
original counterpart, it only avoids holding genome-sized or read-sequence-sized
data in RAM.

  detect_sortbam          replaces debreak_detect.detect_sortbam
                          - split reads no longer keep their full query sequence
                            in the `segmentreads` dict (see LazySeq)
                          - clip info is streamed to disk instead of collected

  split_raw_calls         replaces the rawinscall/rawdelcall/... dicts of the driver
                          - one streaming pass over sv_raw_calls/*.debreak.temp
                            writing per-chromosome, per-type files
  cluster_file,
  cluster_ins_file,
  merge_translocation_file
                          pool workers that read their own split file instead of
                          receiving the whole raw-call list through pickling

  cluster, cluster_ins    replace debreak_merge_cluster.cluster / cluster_ins
                          - the `[0]*contiglength` per-base coverage list (8 bytes
                            per base, ~2 GB for chr1, allocated by every worker)
                            is replaced by an event sweep that only materialises
                            stretches with coverage >= 3
                          - region assignment uses bisect instead of an n*m scan

  identify_duplication    replaces debreak_resdup_selfalignment.identify_duplication
                          - reference slices come from an indexed FASTA
                            (pysam.FastaFile, .fai created on first use) instead of
                            loading the whole genome three times over
"""

import bisect
import os
import shutil
import time

import pysam

import debreak_detect
import debreak_merge_cluster
import debreak_merge_contig
import debreak_resdup_selfalignment


THRESHOLD=3   # `threshold=3` of debreak_merge_cluster.cluster / cluster_ins


# ---------------------------------------------------------------------------
# raw SV signal detection
# ---------------------------------------------------------------------------

class LazySeq(object):
	"""Stand-in for the primary alignment's query sequence stored at index 7 of a
	`segmentreads` entry. debreak_detect.segmentdeletion only slices it to extract
	the inserted sequence of an I-segment call, so the sequence is fetched from
	the BAM on first use instead of being kept for every split read."""
	__slots__=('samfile','chrom','start','name','seq')

	def __init__(self,samfile,chrom,start,name):
		self.samfile=samfile
		self.chrom=chrom
		self.start=start
		self.name=name
		self.seq=None

	def _load(self):
		seq=''
		for align in self.samfile.fetch(self.chrom,self.start,self.start+1):
			if align.reference_start==self.start and align.query_name==self.name and not align.is_supplementary and not align.is_secondary:
				seq=align.query_sequence
				break
		self.seq=seq if seq is not None else ''

	def __getitem__(self,item):
		if self.seq is None:
			self._load()
		return self.seq[item]

	def __len__(self):
		if self.seq is None:
			self._load()
		return len(self.seq)


def detect_sortbam(filename,writepath,min_size,max_size,chrom,chromosomes,record_clip,ifrescuedup,iftumor):
	# Copy of debreak_detect.detect_sortbam; the differences are marked "# lowmem".
	print ('Start to detect SV from '+filename)
	samfile=pysam.AlignmentFile(filename,"rb")
	seqfile=pysam.AlignmentFile(filename,"rb")   # lowmem: separate handle for LazySeq look-ups
	tempfile=open(writepath+'sv_raw_calls/'+filename.split('/')[-1]+'-'+chrom+'.debreak.temp','w')
	clipfile=None   # lowmem: clip info is written as it is found instead of collected in a list
	if record_clip:
		clipfile=open(writepath+'debreak_ins_workspace/readinfo_start_end_'+chrom,'w')
	allreads=samfile.fetch(chrom)
	segmentreads={}
	totalmaplength=0
	for align in allreads:
		if align.is_unmapped or align.is_secondary:
			continue
		readname=align.query_name
		flag=align.flag
		position=align.reference_start+1
		mappingquality=align.mapping_quality
		cigar=align.cigar
		readclipinfo=chrom
		readstart=align.reference_start
		readstop=align.reference_end
		if record_clip and ( cigar[0][0] in [4,5] and cigar[0][1]>=200 ) or (cigar[-1][0] in [4,5] and  cigar[-1][1]>=200):
			if (cigar[0][0]==4 or cigar[0][0]==5) and  cigar[0][1]>=200:
				readclipinfo+='\t0\t'+str(readstart)
			else:
				readclipinfo+='\t'+str(readstart)+'\t0'
			if (cigar[-1][0]==4 or cigar[-1][0]==5) and  cigar[-1][1]>=200:
				readclipinfo+='\t0\t'+str(readstop)
			else:
				readclipinfo+='\t'+str(readstop)+'\t0'
			if clipfile is not None:   # lowmem
				clipfile.write(readclipinfo+'\n')

		if align.is_supplementary:
			cigarinfo=[0,0,0]
			cigarleft=align.cigar[0]
			cigarright=align.cigar[-1]
			if cigarleft[0]==4 or cigarleft[0]==5:
				cigarinfo[0]=cigarleft[1]
			if cigarright[0]==4 or cigarright[0]==5:
				cigarinfo[2]=cigarright[1]
			cigarinfo[1]=align.query_alignment_length
			refend=align.reference_end+1
			readinfo=[readname,flag,chrom,int(position),refend,cigarinfo,mappingquality,'']
			if readname in segmentreads:
				segmentreads[readname]+=[readinfo]
			else:
				segmentreads[readname]=[readinfo]
		else:
			cigar=align.cigarstring
			readseq=align.query_sequence
			cigarinfo=debreak_detect.cigardeletion(flag,chrom,position,cigar,min_size,max_size)
			cigarsv=cigarinfo[0]
			totalmaplength+=align.query_length
			rawsvid=0
			for d in cigarsv:
				rawsvid+=1
				if 'I-cigar' in d:
					insertseq=readseq[d[4]:d[4]+d[2]]
					if  ifrescuedup:
						rawseq=readseq[:d[4]]+readseq[d[4]+d[2]:]
					else:
						rawseq=''
					tempfile.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'_cigar'+str(rawsvid)+'\t'+str(flag)+'\t'+str(mappingquality)+'\t'+str(d[4])+'\t'+insertseq+'\t'+rawseq+'\n')
				else:
					tempfile.write(d[0]+'\t'+str(d[1])+'\t'+str(d[2])+'\t'+d[3]+'\t'+readname+'_cigar'+str(rawsvid)+'\t'+str(flag)+'\t'+str(mappingquality)+'\n')

			if align.has_tag("SA"):
				refend=position+cigarinfo[1]
				readinfo=[readname,flag,chrom,position,refend,cigarinfo[2],mappingquality,LazySeq(seqfile,chrom,align.reference_start,readname)]   # lowmem
				if readname in segmentreads:
					segmentreads[readname]+=[readinfo]
				else:
					segmentreads[readname]=[readinfo]
				suppaligns=align.get_tag("SA").split(';')[:-1]
				for c in suppaligns:
					if c.split(',')[0]==chrom or c.split(',')[0] not in chromosomes:
						continue
					if c.split(',')[2]=='+':
						suppflag=2048
					else:
						suppflag=2064
					cigarinfo=debreak_detect.simplifycigar(c.split(',')[3])
					suppinfo=[readname,suppflag,c.split(',')[0],int(c.split(',')[1]),int(c.split(',')[1])+cigarinfo[1],cigarinfo,c.split(',')[4]]
					segmentreads[readname]+=[suppinfo]
	for readgroup in segmentreads:
		if len(segmentreads[readgroup])<2 or len(segmentreads[readgroup])>20:
			continue
		if iftumor:
			segmentsv=debreak_detect.segmentdeletion_tumor(segmentreads[readgroup],min_size,max_size)
		else:
			segmentsv=debreak_detect.segmentdeletion(segmentreads[readgroup],min_size,max_size)
		for d in segmentsv:
			tempfile.write(d+'\n')
	tempfile.close()
	if clipfile is not None:
		clipfile.close()
	seqfile.close()
	samfile.close()
	print (chrom+' is done with maplength: '+str(totalmaplength))
	if totalmaplength!=0:
		f=open(writepath+'map_depth/maplength_'+chrom,'w')
		f.write(str(totalmaplength)+'\n')
		f.close()
	return True


# ---------------------------------------------------------------------------
# raw calls -> per-chromosome / per-type files (replaces the driver's dicts)
# ---------------------------------------------------------------------------

def split_raw_calls(writepath,vcflist,chromosomes):
	"""One streaming pass over the raw-call files. For every chromosome it writes
	exactly the lines the driver used to keep in rawinscall / rawdelcall /
	rawdupcall / rawinvcall / rawtracall into
	    <writepath>sv_raw_calls/lowmem_split/<chrom>.<ins|del|dup|inv|tra>
	and returns that directory."""
	splitdir=writepath+'sv_raw_calls/lowmem_split/'
	if os.path.exists(splitdir):
		shutil.rmtree(splitdir)
	os.makedirs(splitdir)
	handles={}
	def out(chrom,svtype):
		key=(chrom,svtype)
		if key not in handles:
			handles[key]=open(splitdir+chrom+'.'+svtype,'w')
		return handles[key]
	for filename in vcflist:
		path=writepath+'sv_raw_calls/'+filename
		if not os.path.exists(path):
			continue
		f=open(path,'r')
		for event in f:
			event=event.rstrip('\n')
			if event=='':
				continue
			event_chrom=event.split('\t')[0]
			# classification mirrors the driver: substring tests on the whole line;
			# an insertion line is reduced to its first 7 columns (no sequences)
			if 'I-' in event:
				event=event.split('\t')
				out(event_chrom,'ins').write(event[0]+'\t'+event[1]+'\t'+event[2]+'\t'+event[3]+'\t'+event[4]+'\t'+event[5]+'\t'+event[6]+'\n')
				continue
			if 'D-' in event:
				out(event_chrom,'del').write(event+'\n')
			if 'DUP-' in event:
				out(event_chrom,'dup').write(event+'\n')
			if 'INV-' in event:
				out(event_chrom,'inv').write(event+'\n')
			if 'TRA-' in event:
				event=event.split('\t')
				out(event_chrom,'tra').write(event[0]+'\t'+event[1]+'\t'+event[2]+'\t'+event[3]+'\t'+event[7]+'\t'+event[5]+'\n')
		f.close()
	for h in handles.values():
		h.close()
	return splitdir


def remove_split(splitdir):
	shutil.rmtree(splitdir,ignore_errors=True)


def _read_split(splitdir,chrom,svtype):
	path=splitdir+chrom+'.'+svtype
	if not os.path.exists(path):
		return []
	f=open(path,'r')
	lines=[c.rstrip('\n') for c in f]
	f.close()
	return lines


def cluster_file(splitdir,outpath,chrom,contiglength,mins,maxdepth,svtype,minmapq):
	return cluster(outpath,_read_split(splitdir,chrom,svtype),chrom,contiglength,mins,maxdepth,svtype,minmapq)


def cluster_ins_file(splitdir,outpath,chrom,contiglength,mins,maxdepth,minmapq):
	return cluster_ins(outpath,_read_split(splitdir,chrom,'ins'),chrom,contiglength,mins,maxdepth,minmapq)


def merge_translocation_file(splitdir,minsupp,outpath,chrom,upper_bound,minqual):
	return debreak_merge_contig.merge_translocation(minsupp,outpath,_read_split(splitdir,chrom,'tra'),chrom,upper_bound,minqual)


# ---------------------------------------------------------------------------
# clustering without a per-base list of the whole contig
# ---------------------------------------------------------------------------

def _slice_bounds(lo,hi,length):
	"""Normalise [lo:hi] exactly like Python does for a list of `length`
	elements (the original relies on that when an interval starts before
	position 1 or runs past the contig end). Returns None for an empty slice."""
	if lo<0:
		lo+=length
		if lo<0:
			lo=0
	elif lo>length:
		lo=length
	if hi<0:
		hi+=length
		if hi<0:
			hi=0
	elif hi>length:
		hi=length
	if hi<=lo:
		return None
	return (lo,hi)


def _coverage_runs(intervals,contiglength):
	"""Sparse equivalent of the dense `genomeposition` list. Yields
	(offset, depths) for every maximal stretch of positions covered by at
	least THRESHOLD intervals; `depths` is the per-base coverage of that
	stretch. Outside those stretches the original block-detection loop can
	neither start nor continue a block, so nothing else has to exist."""
	events=[]
	for lo,hi in intervals:
		b=_slice_bounds(lo,hi,contiglength)
		if b is None:
			continue
		events.append((b[0],1))
		events.append((b[1],-1))
	events.sort()
	depth=0
	prev=0
	run_off=None
	run=[]
	i=0
	n=len(events)
	while i<n:
		pos=events[i][0]
		if pos>prev:
			if depth>=THRESHOLD:
				if run_off is None:
					run_off=prev
				run.extend([depth]*(pos-prev))
			elif run_off is not None:
				yield run_off,run
				run_off=None
				run=[]
		while i<n and events[i][0]==pos:
			depth+=events[i][1]
			i+=1
		prev=pos
	if run_off is not None:
		yield run_off,run


def _close_block(svregion,start,end,maxdep,localdep,maxdepth):
	# identical to the block-closing branch of debreak_merge_cluster.cluster
	if maxdep<=maxdepth:
		peakpos=localdep.index(maxdep)
		peakleftsize=0
		for i in range(peakpos):
			if localdep[peakpos-i-1]>=maxdep/10.0:
				peakleftsize+=1
			else:
				break
		svregion+=[(start+peakpos-peakleftsize,end,maxdep)]


def _find_regions(intervals,contiglength,maxdepth):
	"""Same regions as the dense sweep in debreak_merge_cluster.cluster /
	cluster_ins, computed only over stretches with coverage >= THRESHOLD."""
	svregion=[]
	for offset,dep in _coverage_runs(intervals,contiglength):
		inblock=False
		start=0; maxdep=0; localdep=[]
		for i in range(len(dep)):
			d=dep[i]
			if inblock:
				if d>=max(maxdep/10.0,THRESHOLD):
					localdep+=[d]
					if d>maxdep:
						maxdep=d
				else:
					inblock=False
					_close_block(svregion,offset+start,offset+i,maxdep,localdep,maxdepth)
			else:
				if d>THRESHOLD:
					inblock=True
					localdep=[d]
					start=i
					maxdep=d
		if inblock:
			# the position right after the run has depth < THRESHOLD, which closes
			# the block there - unless the run reaches the end of the contig, where
			# the original loop simply ends without ever closing the block
			end=offset+len(dep)
			if end<contiglength:
				_close_block(svregion,offset+start,end,maxdep,localdep,maxdepth)
	return [c for c in svregion if c[2]<maxdepth]


def _assign_to_regions(allsv,svregion,span):
	"""Same result as the original n*m overlap scan (every sv goes to every
	region it overlaps, in allsv order). Regions are disjoint and sorted, so
	bisect finds the candidates; the full scan is kept as a safety net."""
	groups=[[] for d in svregion]
	ends=[d[1] for d in svregion]
	sorted_ok=all(svregion[k][0]<svregion[k+1][0] and ends[k]<=ends[k+1] for k in range(len(svregion)-1))
	for c in allsv:
		start,end=span(c)
		if sorted_ok:
			j=bisect.bisect_right(ends,start)
			while j<len(svregion) and svregion[j][0]<end:
				d=svregion[j]
				if min(end,d[1])-max(d[0],start)>0:
					groups[j].append(c)
				j+=1
		else:
			for j in range(len(svregion)):
				d=svregion[j]
				if min(end,d[1])-max(d[0],start)>0:
					groups[j].append(c)
	return groups


def _cluster_large(allsv,mins,svtypeinfo):
	# identical to the "Large DEL" / "Large INS" section of the original
	largesv=[c for c in allsv if  int(c.split('\t')[2])>2000]
	window=1600
	largesv.sort(key=debreak_merge_cluster.counttimesort)
	largedel=[]
	start=0
	candi=[]
	for event in largesv:
		if int(event.split('\t')[1])<=start+window:
			candi+=[event]
			continue
		if len(candi)>=mins:
			largedel+=debreak_merge_cluster.mergeinfo_insertion(candi,mins,svtypeinfo)
		candi=[event]
		start=int(event.split('\t')[1])
	if len(candi)>=mins:
		largedel+=debreak_merge_cluster.mergeinfo_insertion(candi,mins,svtypeinfo)
	return largedel


def _drop_large_covered(largedel,sv):
	# identical to the de-duplication of large SVs against region SVs in the original
	newsv=[]
	for c in largedel:
		testif=0
		for d in sv:
			if min(int(c.split('\t')[1])+int(c.split('\t')[2]), int(d.split('\t')[1])+int(d.split('\t')[2])) - max(int(c.split('\t')[1]),int(d.split('\t')[1]))>0 and 0.8*int(d.split('\t')[2])<=int(c.split('\t')[2])<=int(d.split('\t')[2])/0.8:
				testif=1; break
		if testif==0:
			newsv+=[c]
	return newsv


def _write_cluster(path,mode,newsv,minmapq):
	if newsv!=[]:
		f=open(path,mode)
		for c in newsv:
			if minmapq==None:
				if float(c.split('\t')[5])<1 or float(c.split('\t')[5])>=10:
					f.write(c+'\n')
			else:
				if float(c.split('\t')[5])>=minmapq:
					f.write(c+'\n')
		f.close()


def cluster(outpath,allsv,chrom,contiglength,mins,maxdepth,svtype,minmapq):
	if allsv==[]:
		return 0
	t1=time.time()
	svtypeinfo={'del':'Deletion','dup':'Duplication','inv':'Inversion'}[svtype]
	largedel=_cluster_large(allsv,mins,svtypeinfo)

	allsv=[c for c in allsv if int(c.split('\t')[2])<=3000]
	intervals=[]
	for c in allsv:
		start=int(c.split('\t')[1])
		end=int(c.split('\t')[1])+int(c.split('\t')[2])
		intervals+=[(start-1,end-1)]
	svregion=_find_regions(intervals,contiglength,maxdepth)

	def span(c):
		start=int(c.split('\t')[1])
		return (start,start+int(c.split('\t')[2]))
	sv=[]
	for svinfo in _assign_to_regions(allsv,svregion,span):
		sv+=debreak_merge_cluster.mergeinfo_insertion(svinfo,mins,svtypeinfo)

	newsv=_drop_large_covered(largedel,sv)
	newsv+=sv
	newsv.sort(key=debreak_merge_cluster.counttimesort)

	t2=time.time()
	print (chrom,'-',svtype,' finished, time cost: ',t2-t1)
	_write_cluster(outpath+svtype+'-merged-cluster-'+chrom,'w',newsv,minmapq)
	return 0


def cluster_ins(outpath,allsv,chrom,contiglength,mins,maxdepth,minmapq):
	if allsv==[]:
		return 0
	t1=time.time()
	svtypeinfo='Insertion'
	largedel=_cluster_large(allsv,mins,svtypeinfo)

	allsv=[c for c in allsv if  int(c.split('\t')[2])<=3000]
	intervals=[]
	for c in allsv:
		start=int(c.split('\t')[1])-100
		end=int(c.split('\t')[1])+100
		intervals+=[(start-1,end-1)]
	svregion=_find_regions(intervals,contiglength,maxdepth)

	def span(c):
		start=int(c.split('\t')[1])-50
		return (start,start+100)
	sv=[]
	groups=_assign_to_regions(allsv,svregion,span)
	for k in range(len(svregion)):
		c=svregion[k]
		mergedins=debreak_merge_cluster.mergeinfo_insertion(groups[k],mins,svtypeinfo)
		for m in mergedins:
			sv+=[m+'\t'+chrom+'\t'+str(c[0])+'\t'+str(c[1])+'\t'+str(c[2])]

	newsv=_drop_large_covered(largedel,sv)
	newsv+=sv
	newsv.sort(key=debreak_merge_cluster.counttimesort)

	t2=time.time()
	print (chrom,'-ins finished, time cost: ',t2-t1)
	_write_cluster(outpath+'ins-merged-cluster-'+chrom,'a',newsv,minmapq)
	return 0


# ---------------------------------------------------------------------------
# DUP rescue with an indexed reference
# ---------------------------------------------------------------------------

def _ref_slice(fasta,chrom,a,b):
	"""fasta[chrom][a:b] with the same semantics as slicing the in-memory
	chromosome string in the original (negative or overflowing bounds)."""
	bounds=_slice_bounds(a,b,fasta.get_reference_length(chrom))
	if bounds is None:
		return ''
	return fasta.fetch(chrom,bounds[0],bounds[1])


def identify_duplication(vcflist,writepath,refpath):
	# Copy of debreak_resdup_selfalignment.identify_duplication; the only change is
	# that reference slices are fetched from an indexed FASTA (pysam creates the
	# .fai next to the reference if it is missing) instead of loading the whole
	# genome into memory.
	ins=open(writepath+'insertion-merged','r').read().split('\n')[:-1]

	reffasta=pysam.FastaFile(refpath)   # lowmem

	os.system("mkdir "+writepath+"debreak_resdup_insertseq/")
	os.system("mkdir "+writepath+"debreak_resdup_refseq/")
	os.system("mkdir "+writepath+"debreak_resdup_map_space/")

	insread={}
	for sv in ins:
		insinfo=sv.split('\t')[0]+'__'+sv.split('\t')[1]+'__'+sv.split('\t')[2]+'__'+sv.split('\t')[3]
		readnames=sv.split('\t')[6].split(';')
		for d in readnames:
			insread[d]=insinfo

		f=open(writepath+"debreak_resdup_refseq/"+insinfo+'.fa','w')
		f.write('>'+insinfo+'\n')
		f.write(_ref_slice(reffasta,sv.split('\t')[0],int(sv.split('\t')[1])-2000-int(sv.split('\t')[2]),int(sv.split('\t')[1])+2000+int(sv.split('\t')[2]))+'\n')   # lowmem
		f.close()
	reffasta.close()

	for vcffile in vcflist:
		f=open(writepath+'sv_raw_calls/'+vcffile,'r')
		a=f.readline()
		while a!='':
			if 'I-cigar' not in a or a.split('\t')[4] not in insread:
				a=f.readline(); continue
			readname=a.split('\t')[0]+'_'+a.split('\t')[1]+'_'+a.split('\t')[2]+'_'+a.split('\t')[4]
			insinfo=insread[a.split('\t')[4]]
			g=open(writepath+'debreak_resdup_insertseq/'+insinfo+'.fa','a')
			insseq=a.split('\t')[8]
			g.write('>'+readname+'\n'+insseq+'\n')
			g.close()
			a=f.readline()
		f.close()

	rescueddup=[]
	trueins=[]

	ff=open(writepath+'debreak_resdup_map_space/maprateinfo','w')
	for sv in ins:
		insinfo=sv.split('\t')[0]+'__'+sv.split('\t')[1]+'__'+sv.split('\t')[2]+'__'+sv.split('\t')[3]
		try:
			os.system('minimap2 -k 9 '+writepath+'debreak_resdup_refseq/'+insinfo+'.fa  '+writepath+'debreak_resdup_insertseq/'+insinfo+'.fa -N 100 > '+writepath+'debreak_resdup_map_space/'+insinfo+'.paf')
			(maprate,bpcorrection)=debreak_resdup_selfalignment.check_dup(insinfo,writepath)
		except:
			(maprate,bpcorrection)=(0.0,0)
		ff.write(insinfo+'\t'+str(maprate)+'\n')
		if maprate>=0.5:
			sv=sv.split('\t')
			sv[1]=str(int(sv[1])-bpcorrection)
			sv='\t'.join(sv)
			rescueddup+=[sv]
		else:
			trueins+=[sv]
	ff.close()

	f=open(writepath+'insertion-merged','w')
	for sv in trueins:
		f.write(sv+'\n')
	f.close()
	dup=open(writepath+'duplication-merged','r').read().split('\n')[:-1]
	rescueddup=[cc.replace('Insertion','Duplication') for cc in rescueddup]
	dup+=rescueddup
	dup.sort(key=debreak_resdup_selfalignment.sortdup)
	f=open(writepath+'duplication-merged','w')
	for sv in dup:
		f.write(sv+'\n')
	f.close()

	return 0
