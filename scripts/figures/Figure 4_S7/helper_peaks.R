source('helper.R')

peakflank <- 0
windowsize <- 1000

# data assumes chr1, pos1, chr2, pos2, svtype, sample_id
svwindow <- function(data, chr, start, end, window=windowsize, extra=peakflank, endsonly=T) {
	stopifnot((end-start)%%window==0)
	stopifnot(sum(is.na(data$chr1))==0)
	stopifnot(sum(colnames(data)!=c('chr1','pos1','chr2','pos2','svtype','sample_id'))==0)
	steps <- (end-start)/window
	output <- data.frame()
	for(i in 1:steps) {
		if(i %% 500 == 0) {
			print(paste(i, '/', steps))
		}
		left <- start+((i-1)*window)
		right <- left+window-1
		output[i,'chrom'] <- chr
		output[i,'start'] <- left
		output[i,'end'] <- right
		onechr <- is.na(data$chr2) | is.na(data$pos2)
		samechr <- data$chr1==data$chr2 & !onechr
		diffchr <- data$chr1!=data$chr2 & !onechr
		if(endsonly) {
			diffchr <- diffchr | samechr
			samechr <- rep(F,length(samechr))
		}
		left <- left-extra
		right <- right+extra
		rowsone <- onechr & inrange(data$pos1,left,right) & data$chr1==chr
		rowssame <- samechr & overlap(data$pos1,data$pos2,left,right) & data$chr1==chr
		rowsdiff <- diffchr &
		((inrange(data$pos1,left,right) & data$chr1==chr) |
		(inrange(data$pos2,left,right) & data$chr2==chr))
		
		datarows <- rowsone | rowssame | rowsdiff
		stopifnot(sum(is.na(datarows))==0)
		output[i,'value'] <- length(unique(data[datarows,'sample_id']))
	}
	return(output)
}


# data assumes chrom, start, end, sample_id
slidingwindow <- function(data, chr, start, end, bysample=T, window=windowsize, extra=peakflank) {
	stopifnot((end-start)%%window==0)
	if(bysample) {
		stopifnot(sum(colnames(data)!=c('chrom','start','end','sample_id'))==0)
	} else {
		stopifnot(sum(colnames(data)!=c('chrom','start','end','value'))==0)
	}
	
	steps <- (end-start)/window
	output <- data.frame()
	for(i in 1:steps) {
		if(i %% 500 == 0) {
			print(paste(i, '/', steps))
		}
		left <- start+((i-1)*window)
		right <- left+window-1
		output[i,'chrom'] <- chr
		output[i,'start'] <- left
		output[i,'end'] <- right

		left <- left-extra
		right <- right+extra
		stopifnot((right-left+1) == (window+extra+extra))
		rows <- overlap(data$start,data$end,left,right) & data$chrom==chr
		stopifnot(sum(is.na(rows))==0)
		if(sum(rows)>0) {
			dataslice <- data[rows,]
			if(bysample) {
				output[i,'value'] <- length(unique(dataslice$sample_id))
			} else {
				#dataslice$start <- pmax(dataslice$start,start)
				#dataslice$end <- pmin(dataslice$start,end)
				#meansize <- mean(dataslice$end-dataslice$start+1)
				#dataslice$weight <- (dataslice$end-dataslice$start+1)/meansize
				#output[i,'value'] <- sum(dataslice$value*dataslice$weight,na.rm=T)
				output[i,'value'] <- sum(dataslice$value,na.rm=T)
			}
		} else {
			output[i,'value'] <- 0
		}
	}
	return(output)
}

findpeaks <- function(sw_bed, threshold, minwidth=10000, within=10000, minsize=3) {
	stopifnot(sum(colnames(sw_bed)!=c('chrom','start','end','value'))==0)
	results <- data.frame()
	inpeak <- F
	currentpeak <- 0
	
	if(max(sw_bed$value)<threshold) {
		return(NULL)
	}
	
	for(i in 1:dim(sw_bed)[1]) {
		stopifnot(!is.na(currentpeak))
		val <- sw_bed[i,'value']
		stopifnot(!is.na(val))
		if(val>=threshold) {
			if(!inpeak) { #new peak?
				inpeak <- T
				newpeakstart <- sw_bed[i,'start']
				if(currentpeak==0) { #first peak
					currentpeak <- currentpeak+1
					results[currentpeak,'start'] <- newpeakstart
					results[currentpeak,'max'] <- val
				} else if(newpeakstart-results[currentpeak,'end']>within) { #new peak!
					currentpeak <- currentpeak+1
					results[currentpeak,'start'] <- newpeakstart
					results[currentpeak,'max'] <- val
				} else {
					#keep current peak, and continue
					results[currentpeak,'end'] <- NA
				}
			} else { 
				#continue peak
				results[currentpeak,'max'] <- max(results[currentpeak,'max'],val)
			}
			if(i==dim(sw_bed)[1]) { #stop if last one
				results[currentpeak,'end'] <- sw_bed[i,'end']
				inpeak <- F
			}
		} else {
			if(inpeak) { #done with peak
				results[currentpeak,'end'] <- sw_bed[i-1,'end']
				inpeak <- F
			} else {
				#do nothing and keep going
			}
		}
		#print(results)
	}
	stopifnot(sum(is.na(results$start))==0)
	stopifnot(sum(is.na(results$end))==0)
	stopifnot(sum(is.na(results$max))==0)
	results$width <- results$end-results$start
	results <- results[results$width>=minwidth & results$max>=minsize,]
	return(results)
}















