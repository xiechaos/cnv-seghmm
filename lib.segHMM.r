check.library <- function(hmm, cbs, refine)
{
	if(hmm | refine == 1) library(RHmm)
	if(hmm) library(plyr)
	if(cbs) library(DNAcopy)
}

hmm.it <- function(data, log2.threshold, min.window) 
{
	require(RHmm)
	require(plyr)
	require(MASS)

	data$hmm <- 0
	data$hmm.size <- NA
	data$hmm.log2 <- NA

	cap <- 21 
	nn.state <- floor(2 ^ (max(abs(data$log2)) + 1)) * 2 - 3; 
	n.state <- max(5, min(cap, nn.state)); 

	normal.state <- (n.state + 1) / 2
	copy <- c(rep(2, normal.state), seq(3, by = 1, length = n.state - normal.state))
	if(n.state == cap)
	{
		copy[n.state] = n.state + 0.7 * (nn.state - n.state) 
	}

#   combine non-significant states
#	combined <- abs(log2(copy/rev(copy))) < abs(log2.threshold) 
#	combined[normal.state] = FALSE
#	n.state <- n.state - length(which(combined))
#	normal.state <- (n.state + 1) / 2
#	copy <- copy[!combined]

	# transision probability
	trans <- matrix(rep(0, n.state ^ 2), nrow = n.state)
	for(i in 1:n.state) trans[i,i] <- 10000000  # self to self
	trans[, normal.state] <- 0.1   # from other state to normal
	trans[normal.state, ] <- 1    # from normal to other state
	trans[normal.state, normal.state] <- 1000000000 * sum(trans[normal.state,])   # from normal to normal
	trans <- aaply(trans, 1, function(x){x/sum(x)})

	# initial probability
	initial <- rep(0, n.state)
	initial[normal.state] <- 1
	initial <- initial / sum(initial)

	# emission
	m <- c()
	v <- c()
	n <- mean(median(data$test), median(data$ref))
	for(i in 1:n.state)
	{
		tcopy <- copy[i]
		rcopy <- rev(copy)[i]
		num <- 100000
		test <- rpois(num, tcopy * n)
		# data must be noisy, so
		test <- test + rnorm(num, 0, tcopy * n/cap/3 * min(tcopy, cap) )
		zero <- which(test<=0)
		if(length(zero))
		{
			test[zero] <- runif(length(zero), 0.001, 0.4)
		}
		ref <- rpois(num, rcopy * n)
		# data must be noisy, so
		ref <- ref + rnorm(num, 0, rcopy * n/cap/3 * min(rcopy, cap) )
		zero <- which(ref<=0)
		if(length(zero))
		{
			ref[zero] <- runif(length(zero), 0.001, 0.4)
		}
		tf <- log2(test/ref)
		est <- fitdistr(tf, 'normal')$estimate
		m <- c(m, est['mean'])
		v <- c(v, est['sd'] ^ 2)
	}

	state2log2 <- log2(copy/rev(copy))
	state2test <- copy
	state2ref <- rev(copy)
	hmm <- HMMSet(initial, trans, 'NORMAL', mean = m, var = v)
#	print(hmm)
	states <- c()

	if(length(data$log2) < 50000)
	{
		# make sure the first observation is in nromal.state
		states <- viterbi(hmm, c(0, data$log2))$states[-1]
	}else
	{
		frag <- 40000
		n.cycle <- 3
#		n.cycle <- 1
		out <- data.frame(matrix(NA, ncol = n.cycle, nrow = length(data$log2)))
		for(cycle in 1:n.cycle)
		{
			start <- 1 + cycle * frag/(n.cycle + 1)
			remaining <- length(data$log2) - start + 1
			starts <- c(1, start + frag * 0:(remaining %/% frag))
			maxs <- c(start - 1, rep(length(data$log2), length(starts)))
			for(i in 1:length(starts))
			{
				x <- starts[i]
				y <- min(starts[i] + frag - 1, maxs[i]) 
				if(y - x < 50) 
				{
					out[x:y, cycle] <- normal.state
					next
				}
				piece <- data$log2[x:y] 
				# make sure the first observation is in nromal.state
				piece.states <- viterbi(hmm, c(0, piece))$states[-1]
				out[x:y, cycle] <- piece.states
			}
		}
		states <- apply(out, 1, median)
	}

	data$hmm.log2 <- state2log2[states]
	for(s in 1:n.state)
	{
		if(s == normal.state) next
#		label <- patch_label(data$hmm.log2 == state2log2[s], 1, max(data$hmm))
		label <- patch_label(data$hmm.log2 == state2log2[s], min.window, max(data$hmm))
		data$hmm <- data$hmm + label
	}
	if(max(data$hmm)>0)
	{
		for(id in seq(1,max(data$hmm)))
		{
			sub <- subset(data, hmm==id)
			start <- ceiling(mean(c(min(sub$start), min(sub$position))))
			end <- floor(mean(c(max(sub$end), max(sub$position))))
			size <- end - start +1
			data[rownames(sub), 'hmm.size'] <- size
		}
	}

	data
}

cnv.print <- function(data, type, file="")
{
	require(plyr)
	field <- c(type, 'chromosome', paste(type, 'size', sep='.'), paste(type, 'log2', sep='.'))
	temp <- data[data[,type]>0, c(field, 'start', 'end')]
	temp <- ddply(temp, field, summarize, start = min(start), end = max(end))
	temp[,type] <- paste(type, temp[,type], sep='_')
	temp$chromosome <- paste('chr', temp$chromosome, sep='')
	temp <- temp[, c(1, 2, 5, 6, 3, 4)]
	write.table(temp, sep="\t", quote=F, file=file, row.names=F)
}

cw.it <- function(data, log2.threshold, minimum.window, norm.factor)
{
	data$cw <- 0
	data$cw.size <- NA
	data$cw.log2 <- 0 
	p.label <- patch_label(data$log2 >= log2.threshold, minimum.window, max(data$cw)) 
	data[,'cw'] <- data[,'cw'] + p.label
	n.label <- patch_label(data$log2 <= -1 * log2.threshold, minimum.window, max(data$cw)) 
	data[,'cw'] <- data[,'cw'] + n.label
	for(i in unique(data$cw))
	{
		if(i == 0) next
		temp <- subset(data, cw == i)
		rows <- rownames(temp)
		data[rows, 'cw.log2'] = log2(sum(temp$test) / sum(temp$ref)) - norm.factor
		data[rows, 'cw.size'] = max(temp$end) - min(temp$start) + 1
	}
	data
}
patch_label <- function(bool, n, l)
{
	rle.bool <- rle(bool)
	good <- rle.bool$values & rle.bool$lengths >= n
	good.id <- which(good)
	rle.bool$values <- rep(0, length(rle.bool$values))
	rle.bool$values[good.id] <- l + seq_along(good.id)
	return(inverse.rle(rle.bool));
}

z2t <- function(z, lambdax, lambday)
{
	(lambday*z-lambdax)/sqrt(lambday*z^2+lambdax)
}

re.cap <- function(data)
{
	# re-value Inf and -Inf values:
	my.max <- max(data[which(data <  Inf)])
	my.min <- min(data[which(data > -Inf)])
	data[which(data ==  Inf)] <- my.max
	data[which(data == -Inf)] <- my.min
	return(data)
}

cnv.cal <- function(file, log2.threshold, cw=FALSE, minimum.window=4, hmm=TRUE, cbs=FALSE)
{
	if(hmm) require('RHmm')
	if(hmm) require('plyr')
	if(cbs) require('DNAcopy')
	data <- read.delim(file, colClasses = c('factor', 'integer', 'integer', 'integer', 'integer'))
	# remove window that is too small or too large
	# like the box plot
	win <- data$end - data$start
	q <- quantile(win, c(0.25, 0.75))
	iqr <- diff(q)
	data <- data[win >= q[1]-iqr & win <= q[2]+iqr, ]
	

	data$p.value <- NA
	data <- transform(data, position = round((end+start)/2), log2 = log2(test/ref))
	norm.factor <- median(data$log2)
	data$log2 <- re.cap(data$log2 - norm.factor)

	t.median <- median(data$test)
	r.median <- median(data$ref)
	rows <- rownames(subset(data, log2 >= 0))
	if(length(rows))
	{
		temp <- data[rows,]
		data[rows, 'p.value'] <- 2 * pnorm(z2t(temp$test/temp$ref,t.median,r.median), lower.tail=FALSE)
	}
	rows <- rownames(subset(data, log2 < 0))
	if(length(rows))
	{
		temp <- data[rows,]
		data[rows, 'p.value'] <- 2 * pnorm(z2t(temp$test/temp$ref,t.median,r.median), lower.tail=TRUE)
	}

	if(cw)
	{
		cat('Simple consequtive windows ...\n')
		data <- cw.it(data, log2.threshold, minimum.window, norm.factor)
		cat(sprintf("found %d CNV segments\n", max(data$cw)))
	}

	if(hmm)
	{
		cat('HMM-ing ...\n')
		data <- hmm.it(data, log2.threshold, minimum.window);
		cat(sprintf("found %d CNV segments\n", max(data$hmm)))
	}
	if(cbs)
	{
		cat('CBS-ing ...\n')
		data <- cbs.it(data, log2.threshold, minimum.window);
		cat(sprintf("found %d CNV segments\n", max(data$cbs)))
	}
	data
}


cbs.it <- function(data, log2.threshold)
{
	require('DNAcopy');
	cna.obj <- CNA(cbind(data$log2), data$chrom, data$position, data.type='logratio')
	smoothed <- smooth.CNA(cna.obj)
	segments <- segment(smoothed, verbose=0)
	data$cbs.log2 <- rep(segments$output$seg.mean, segments$output$num.mark)
	data$cbs.size <- NA
	data$cbs <- 0
	id <- 0
	for(i in which(abs(segments$out$seg.mean) >= log2.threshold))
	{
#		if(segments$output$num.mark[i] < minimum.window) next
		id <- id + 1
		chr <- segments$output$chrom[i]
		start <- segments$output$loc.start[i]
		end  <- segments$output$loc.end[i]
		sm  <- segments$output$seg.mean[i]
		index <- which(data$chrom == chr & data$position >= start & data$position <= end)
		data[index, 'cbs'] <- id
		data[index, 'cbs.size'] <- end - start + 1
	}
	data
}

hmm.bound <- function(cand.file, test.file, ref.file, total.file, out.file='', unreliable.file = '', header=FALSE) 
{
	require(RHmm)

	cand <- read.delim(cand.file)
	test <- read.delim(test.file)
	ref  <- read.delim( ref.file)
	total <- read.delim(total.file)

	out <- data.frame(id = c(), chromosome = c(), start = c(), end = c(), size= c(), log2 = c())
	unreliable <- data.frame(id = c(), chromosome = c(), start = c(), end = c(), size= c(), log2 = c())
	chr <- unique(cand$chromosome)
# one line only
	n.factor <- total$test / total$ref;
	for(myid in unique(cand$id))
	{
		bad <- F
		bound <- c()
		log2.here <- c()
		start.old <- subset(cand, id == myid & side == 'left')$mid
		end.old   <- subset(cand, id == myid & side == 'right')$mid
		for(myside in c('left', 'right'))
		{
			info <- subset(cand, id == myid & side == myside)
			log2.here <- info$log2.a
			t.read <- subset(test, location >= info$from & location <= info$to)
			r.read <- subset(ref,  location >= info$from & location <= info$to)
			if(nrow(t.read) == 0 | nrow(r.read) == 0)
			{
				bad <- T
				warning('No mapped reads in this region, check the *.unreliable file')
			}else
			{
				bound <- c(bound, hmm.resolve(info, n.factor, t.read, r.read))
			}
		}
		if(length(bound) > 0)
		{
			temp <- data.frame()
			start.now <- bound[1] 
			end.now <- bound[2]
			if(is.na(start.now))
			{
				bad = T
				start.now <- start.old
			}
			if(is.na(end.now))
			{
				bad = T
				end.now <- end.old
			}

			temp  <- data.frame(
				id = myid, 
				chromosome = paste('chr', chr, sep=''), 
				start = start.now, 
				end = end.now, 
				size = end.now - start.now + 1, 
				log2 = log2.here
			)

			if(bad)
			{
				unreliable <- rbind(unreliable, temp)
				warning('can not determine bound reliably, check the *.unreliable file')
			}else
			{
				out <- rbind(out, temp)
				cat(paste(temp[1,]), fill=T)
			}
		}
	}

	head <- c('id', 'chromosome', 'start', 'end', 'size', 'log2')
	if(header) cat(head, sep = '\t', fill = T, file = out.file)
	if(header) cat(head, sep = '\t', fill = T, file = unreliable.file)
	write.table(out, file=out.file, row.names=F, sep='\t', quote=F, append = T, col.names = F)
	write.table(unreliable, file=unreliable.file, row.names=F, sep='\t', quote=F, append = T, col.names = F)
}
			
hmm.resolve <- function(info, n.factor, t.read, r.read)
{
#	print(info)
	# initial probability
	initial <- c(1, 0, 0)

	# emission probability
	za <- 2^info$log2.a * n.factor
#	print(za)
	dpa <- c(za/(1+za), 1/(1+za), 0)
	zb <- 2^info$log2.b * n.factor
	dpb <- c(zb/(1+zb), 1/(1+zb), 0)
	dpe <- c(0, 0, 1)

#	print(dpa)
#	print(dpb)
#	print(dpe)

	# transition probability
	na <- length(c(which(t.read$location <= info$mid[1]), which(r.read$location < info$mid[1])))
	nb <- nrow(t.read) + nrow(r.read) - na
	trans <- matrix(c(
		(na-1)/na, 1/na,       0,
		0,         (nb-1)/nb,  1/nb,
		0,         0,          0), 
		byrow=T, nrow=3)

	labs <- c('a', 'b', 'e') # for the three states
#	print(summary(t.read))
#	print(summary(r.read))
	t.read$obs <- 'a'
	r.read$obs <- 'b'
	data <- rbind(t.read, r.read)
	# randomize the order of reads at the same spot
	data <- data[sample(1:nrow(data)), ]
	data <- data[order(data$location), ]
	max.obs <- 40000
	if(nrow(data) > max.obs)
	{
		offset <- as.integer((nrow(data) - max.obs)/2)
		data <- data[offset:(offset + max.obs -1), ]
	}

	obs <- c(data$obs, 'e')
#	print(summary(as.factor(obs)));

	hmm <- HMMSet(initial, trans, 'DISCRETE', proba=list(dpa, dpb, dpe), labels=labs)
	states <- viterbi(hmm, obs)$states
	
	bound <- c()
	al <- length(which(states == 1))
	bl <- length(which(states == 2))
	if(al/(al+bl) > 0.05 & bl/(al+bl) > 0.05)
	{
		a <- max(which(states == 1))
		b <- min(which(states == 2))
		bound <- as.integer(mean(data[c(a,b), 'location']))
	}else
	{
		bound <- NA
	}

	return(bound)
}

plot.bound <- function(cand.dir, cnv,  i)
{
	this <- subset(cnv, id == i)
	print(this)
	if(nrow(this))
	{
		require(ggplot2)

		cand.f <- list.files(cand.dir, sprintf('^%d.cand', i), full = T)
		test.f <- list.files(cand.dir, sprintf('^%d.test', i), full = T)
		ref.f <- list.files(cand.dir, sprintf('^%d.ref', i), full = T)
		bound <- read.delim(cand.f)
		test <- read.delim(test.f)
		ref <- read.delim(ref.f)
		test$sample <- 'test'
		ref$sample <- 'ref'
		read <- rbind(test, ref)

		left <- subset(bound, side == 'left')$to
		right <- subset(bound, side == 'right')$from

		read <- transform(read,
			left = ifelse(location <= left, T, F),
			right = ifelse(location >= right, T, F)
		)
		read <- melt(read, measure = c('left', 'right'))
		read <- read[read$value, ]

		df <- data.frame(variable = c('left', 'right'), bound = c(this$start, this$end))

		bplot <- ggplot(read) + 
			geom_density(aes(location, fill = sample), alpha = I(0.5), adjust = 0.2) + 
			facet_wrap(~variable, scales = 'free', ncol = 1) + 
			geom_linerange(aes(bound), df, ymin = 0, ymax = Inf,  colour = 'red')

		bplot + theme_bw()

	}else
	{
		warning(sprintf("CNV id %s can not be found", i))
	}
}
plot.cnv <- function(raw, cnv, i)
{
	this <- subset(cnv, id == i)
	print(this)
	if(nrow(this))
	{
		require(ggplot2)

		size <- this$size

		left <- this$start - size
		temp <- cnv[cnv$end < this$start, 'end']
		if(length(temp))
		{
			L <- max(cnv[cnv$end < this$start, 'end'])
			L <- L + (left + size - L) / 3
			left <- max(left, L)
		}

		right <- this$end + size
		temp <- cnv[cnv$start > this$end, 'start']
		if(length(temp))
		{
			R <- min(cnv[cnv$start > this$end, 'start'])
			R <- R - (R - (right - size)) / 3
			right <- min(right, R)
		}

		data <- subset(raw, end >= left & start <= right)
		cplot <- ggplot() + 
			geom_point(aes(position, log2), data) +
			geom_vline(x = c(this$start, this$end), colour = 'red') +
			geom_segment(aes(x = start, y = log2, xend = end, yend = log2), this, colour = 'blue') 
		cplot + theme_bw()
	}else
	{
		warning(sprintf("CNV id %s can not be found", i))
	}
}
