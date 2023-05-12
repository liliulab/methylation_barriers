test.lmer = function(chr, chr.folder,chr.methyl.na.logit.all){
	print(paste(chr.folder, '.test.start', sep=''))
	load(paste(chr.folder, 'chr.tss.island.repeat.2k', sep=''))
   

	un=unique(chr.tss.island.repeat.2k[,c('gene.id','ENST','tss','repeat.start','repeat.end','island.start','island.end')])
	un$tss200.left = as.numeric(un$tss)-200
	un$tss200.right = as.numeric(un$tss)+200

	lmer.topo.chr=c()
	for(i in 1:nrow(un)){
	#cat('processing...', i, '...', i, '\n'); flush.console(); 
		tss200L = un[i,8]
		tss200R = un[i,9]
		repeat.start = un[i,4]
		repeat.end = un[i,5]
		island.start = un [i,6]
		island.end = un [i,7]
		pp = un[i, 4:9]
		temp = chr.methyl.na.logit.all[which(chr.methyl.na.logit.all$pos >= min(pp) & chr.methyl.na.logit.all$pos <= max(pp)), ]
		this.tss200.me=temp[which(temp$pos>=tss200L&temp$pos<=tss200R),]
		this.repeat.me=temp[which(temp$pos>=repeat.start &temp$pos<=repeat.end),]
		this.island.me = temp[which(temp$pos>=island.start &temp$pos<=island.end),]
		if(nrow(this.tss200.me)>0 & nrow(this.repeat.me)>0 & nrow(this.island.me)>0){
			this.tss200.me$location ='tss200'
			this.repeat.me$location = 'repeats'
			this.island.me$location = 'island'
			this.methyl = rbind(this.tss200.me, this.repeat.me,this.island.me)
			this.methyl.long =reshape(this.methyl, 
				direction = "long",
				varying = list( c(colnames(this.methyl[4:32]))),
				v.names = 'methyl.value',
				times = c(colnames(this.methyl[4:32])),
				timevar ='clinic'
			)
			this.methyl.long$samples = this.methyl.long$clinic
			this.methyl.long$clinic = 'cancer'	
			this.methyl.long[grep('N.methyl', this.methyl.long$samples), 'clinic'] = 'normal'
			this.methyl.long$pos = factor(this.methyl.long$pos)
			samples = this.methyl.long[,'samples']
			this.methyl.long$patient.id =gsub("[^[:digit:]]", "", this.methyl.long$samples)
			this.methyl.long$samples=factor(this.methyl.long$samples)
			this.methyl.long$patient.id=factor(this.methyl.long$patient.id)
			a = this.methyl.long[which(!is.na(this.methyl.long$methyl.value)), ]	
			tb = table(a$location, a$clinic)
				
			if (nrow(tb) > 2 & ncol(tb) > 1&sum(colSums(tb ==0))==0 & length(unique(a$methyl.value)) > 5 & length(which(tb < 4)) == 0){
				this.lmer = lmer(methyl.value~location+clinic+location:clinic+(1|patient.id),data = a)
		
				s1= data.frame(anova(this.lmer))
			
				p.location = s1[1,6]
				p.clinic = s1[2,6]
				p.inter = s1[3,6]
				
				s2= data.frame(lsmeans(this.lmer, pairwise ~ clinic|location,lmer.df = "asymp")[2])
				s3= data.frame(lsmeans(this.lmer, pairwise ~ location|clinic,lmer.df = "asymp")[2])
				   
				s = cbind(p.location, p.clinic, p.inter, s2[1,c(3,7)],s2[3,c(3,7)],s3[1,c(3,7)],s3[3,c(3,7)],s3[4,c(3,7)],s3[6,c(3,7)])
				this.lmer.sum=cbind(chr, un[i,],s)
				lmer.topo.chr= rbind(lmer.topo.chr,this.lmer.sum)
					}
			}
	}

	#dim(lmer.topo.chr)
	
	colnames(lmer.topo.chr) = c('chr',colnames(un),'p.location','p.clinic','p.inter','is.ca_is.nor.esti','is.ca_is.nor.p','ts.ca_ts.nor.esti','ts.ca_ts.nor.p','is.ca_rep.ca.esti',
	'is.ca_rep.ca.p','rep.ca_ts.ca.esti','rep.ca_ts.ca.p','is.nor_rep.nor.esti','is.nor_rep.nor.p', 'rep.nor_ts.nor.esti','rep.nor_ts.nor.p'  )

	lmer.topo.chr$adj.p.location=p.adjust(lmer.topo.chr$p.location,method='BH',n=length(lmer.topo.chr$p.location))
	lmer.topo.chr$adj.p.clinic=p.adjust(lmer.topo.chr$p.clinic,method='BH',n=length(lmer.topo.chr$p.clinic))
	lmer.topo.chr$adj.p.inter=p.adjust(lmer.topo.chr$p.inter,method='BH',n=length(lmer.topo.chr$p.inter))

	#dim(lmer.topo.chr)

	outfile = paste(chr.folder, 'lmer.out.logit.RData', sep='')
	save(lmer.topo.chr, file=outfile)
	print(paste(chr, '.test.finish', sep=''))
	return(lmer.topo.chr)
	
}


## input: 
## output: 
filter.lmer = function(lmer.topo.chr, chr.folder) {
	print(paste(chr.folder, '.filter.start', sep=''))
	u.sig.lmer.chr = lmer.topo.chr[which(lmer.topo.chr$adj.p.location <=0.1 & lmer.topo.chr$adj.p.clinic<=0.1 & lmer.topo.chr$adj.p.inter<=0.1), ]
	a=u.sig.lmer.chr[which(u.sig.lmer.chr$ts.ca_ts.nor.esti>0 & u.sig.lmer.chr$ts.ca_ts.nor.p < 0.05),]
	b = a[which(a$is.ca_is.nor.esti>0& a$is.ca_is.nor.p < 0.05),]
	c =b[which(b$rep.nor_ts.nor.esti >0&b$rep.nor_ts.nor.p <0.05),]
	d= c[which(c$is.nor_rep.nor.esti<0 &c$is.nor_rep.nor.p<0.05),]
	e = d[which((d$rep.ca_ts.ca.esti>0 & d$rep.ca_ts.ca.p < 0.05) | d$rep.ca_ts.ca.p > 0.05),]
	f= e[which((e$is.ca_rep.ca.esti<0 & e$is.ca_rep.ca.p<0.05) |e$is.ca_rep.ca.p >0.05),]
	lmer.topo.filtered = f
	
	outfile = paste(chr.folder, 'lmer.filtered.logit.RData', sep='')
	save(lmer.topo.filtered, file=outfile)
	print(paste(chr.folder, '.filter.finish', sep=''))
	return(lmer.topo.filtered)
}

fold.lmer = function(lmer.topo.filtered, fc.cutoff=4, chr.folder,chr.methyl.na) {
	print(paste(chr.folder, '.fold.start', sep=''))
	u.sig.2=lmer.topo.filtered
	methyl.dif.4=c()
	fc.cutoff = 4
	for(u in 1:nrow(u.sig.2)) {
		#cat('processing', u, '...\n')
		tss200L=u.sig.2[u,9]
		tss200R= u.sig.2[u,10]
		repeat.start=u.sig.2[u,5]
		repeat.end=u.sig.2[u,6]
		island.start=u.sig.2[u,7]
		island.end=u.sig.2[u,8]
		
		normal.col = grep('N.methyl', colnames(chr.methyl.na))
		cancer.col = grep('.methyl', colnames(chr.methyl.na))
		cancer.col = cancer.col[-which(cancer.col %in% normal.col)]
		
		pp = u.sig.2[u, 5:10]
		temp = chr.methyl.na[which(chr.methyl.na$pos >= min(pp) & chr.methyl.na$pos <= max(pp)), ]
		
		normal.tss200.avg = apply(temp[which(temp$pos>=tss200L&temp$pos<=tss200R), normal.col], 1, mean, na.rm=T)
		normal.island.avg = apply(temp[which(temp$pos>=island.start&temp$pos<=island.end), normal.col], 1, mean, na.rm=T)
		normal.repeats.avg = apply(temp[which(temp$pos>=repeat.start&temp$pos<=repeat.end), normal.col], 1, mean, na.rm=T)
		cancer.island.avg = apply(temp[which(temp$pos>=island.start&temp$pos<=island.end), cancer.col], 1, mean, na.rm=T)
		cancer.repeats.avg = apply(temp[which(temp$pos>=repeat.start&temp$pos<=repeat.end), cancer.col], 1, mean, na.rm=T)
		if(!is.na(mean(cancer.island.avg, na.rm=T))
			&&!is.na(mean(normal.island.avg, na.rm=T))
			&& mean(cancer.island.avg, na.rm=T) >= (fc.cutoff * mean(normal.island.avg, na.rm=T)) 
			&&!is.na(mean(normal.repeats.avg, na.rm=T))
			&&!is.na(mean(normal.tss200.avg, na.rm=T))
			&& mean(normal.repeats.avg, na.rm=T) >= (fc.cutoff * mean(normal.tss200.avg, na.rm=T))) 
			{
			this.methyl=u.sig.2[u,]
			#cat(u,'this.methyl', unlist(this.methyl[,c(3:7)]), '\n'); flush.console()
			methyl.dif.4 = rbind(methyl.dif.4, this.methyl)
			}
	}
	methyl.dif.4.lmer.topo.chr=methyl.dif.4
	
	outfile = paste(chr.folder, 'lmer.filtered.fc.logit.RData', sep='')
	save(methyl.dif.4.lmer.topo.chr, file=outfile)
	print(paste(chr.folder, '.fold.finish', sep=''))
	return(methyl.dif.4.lmer.topo.chr)
}
#####Failed in batch with 3 steps above, so I used 'batch topo neg'code##############
neg.lmer = function(lmer.topo.chr, chr.folder,chr.methyl.na) {
	print(paste(chr, '.neg.lmer', sep=''))
	u.sig.2 = lmer.topo.chr[which((lmer.topo.chr$adj.p.location >0.1 & lmer.topo.chr$adj.p.clinic >0.1 & lmer.topo.chr$adj.p.inter >0.1)|(lmer.topo.chr$p.location < 0.001 & lmer.topo.chr$p.clinic>0.1  & lmer.topo.chr$p.inter>0.1 & lmer.topo.chr$rep.nor_ts.nor.esti < 0 & lmer.topo.chr$is.nor_rep.nor.esti > 0 & lmer.topo.chr$rep.ca_ts.ca.esti < 0 &lmer.topo.chr$is.ca_rep.ca.esti > 0)), ]
	neg.frag.topo=c()
	cut.low = 0.8
	cut.high = 1.2
	for(u in 1:nrow(u.sig.2)) {
		#cat('processing', u, '...\n')
		tss200L=u.sig.2[u,9]
		tss200R= u.sig.2[u,10]
		repeat.start=u.sig.2[u,5]
		repeat.end=u.sig.2[u,6]
		island.start=u.sig.2[u,7]
		island.end=u.sig.2[u,8]
		
		normal.col = grep('N.methyl', colnames(chr.methyl.na))
		cancer.col = grep('.methyl', colnames(chr.methyl.na))
		cancer.col = cancer.col[-which(cancer.col %in% normal.col)]
		
		pp = u.sig.2[u, 5:10]
		temp = chr.methyl.na[which(chr.methyl.na$pos >= min(pp) & chr.methyl.na$pos <= max(pp)), ]
		
		normal.tss200.avg = apply(temp[which(temp$pos>=tss200L&temp$pos<=tss200R), normal.col], 1, mean, na.rm=T)
		normal.island.avg = apply(temp[which(temp$pos>=island.start&temp$pos<=island.end), normal.col], 1, mean, na.rm=T)
		normal.repeats.avg = apply(temp[which(temp$pos>=repeat.start&temp$pos<=repeat.end), normal.col], 1, mean, na.rm=T)
		cancer.island.avg = apply(temp[which(temp$pos>=island.start&temp$pos<=island.end), cancer.col], 1, mean, na.rm=T)
		cancer.repeats.avg = apply(temp[which(temp$pos>=repeat.start&temp$pos<=repeat.end), cancer.col], 1, mean, na.rm=T)
		if(!is.na(mean(cancer.island.avg, na.rm=T))
		  &&!is.na(mean(normal.island.avg, na.rm=T))
		  && mean(cancer.island.avg, na.rm=T) >= (cut.low * mean(normal.island.avg, na.rm=T))
		  && mean(cancer.island.avg, na.rm=T) <= (cut.high * mean(normal.island.avg, na.rm=T))
		  &&!is.na(mean(normal.repeats.avg, na.rm=T))
		  &&!is.na(mean(normal.tss200.avg, na.rm=T))
		  && mean(normal.repeats.avg, na.rm=T) >= (cut.low * mean(normal.tss200.avg, na.rm=T))
		  && mean(normal.repeats.avg, na.rm=T) <= (cut.high* mean(normal.tss200.avg, na.rm=T)))	  {
			this.methyl=u.sig.2[u,]
			neg.frag.topo = rbind(neg.frag.topo, this.methyl)
	}
}
	outfile = paste(chr.folder, 'neg.lmer.topo.logit.RData', sep='')
	save(neg.frag.topo, file=outfile)
	print(paste(chr, '.neg.lmer.topo.finish', sep=''))
	return(neg.frag.topo)
}



all.steps = function(chr) {
	chr.folder = paste(chr, '/', sep='')
	varname = paste(chr, '.methyl', sep='')
	load(paste(chr.folder, varname, sep=''))
	chr.methyl = get(varname)
	#cat(chr, 'chr.tss.island.repeat', nrow(chr.tss.island.repeat), '\n'); flush.console()
	#head(chr.tss.island.repeat,2)
	chr.methyl.na = chr.methyl[rowSums(!is.na(chr.methyl)) >= 6,]
	
	chr.methyl.na.logit = data.frame(logit(chr.methyl.na[,4:32],adjust = 0.001))
	chr.methyl.na.logit.all = cbind( chr.methyl.na[,1:3],  chr.methyl.na.logit)
	
	
	#test = T
	#if(test) {
		#cat(chr, '\n');
		#temp = data.frame(chr=chr, message='test')
		#write.table(temp, paste(chr, 'test.txt', sep='.'))
	#} else {
		lmer.topo.chr = test.lmer(chr=chr, chr.folder=chr.folder,chr.methyl.na.logit.all=chr.methyl.na.logit.all )
		lmer.topo.filtered = filter.lmer(lmer.topo.chr, chr.folder=chr.folder)
		lmer.topo.fold = fold.lmer(lmer.topo.filtered, fc.cutoff=4, chr.folder=chr.folder, chr.methyl.na=chr.methyl.na)
		neg.frag.topo = neg.lmer(lmer.topo.chr, chr.folder=chr.folder,chr.methyl.na=chr.methyl.na)
	}



library(parallel)
library(Matrix)
library(lme4) 
library(lmerTest)
library(emmeans)
library(lsmeans)
library(carData)
library(car)
setwd('/scratch/jshu4/Protector/sum')
chrs = as.list(c(paste(rep('chr',22), 1:22, sep=''), 'chrX','chrY'))
m = mclapply(X=chrs, FUN=all.steps, mc.cores=24)


# sbatch batch.trio.sh











