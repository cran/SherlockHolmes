#library(stringr)

#' grabFunctionParameters
#' 
#' @description retrieve capture all of the parameter names and values passed in
#' 
#' @details copied and pasted from
#' https://stackoverflow.com/questions/66329835/using-r-how-to-get-all-parameters-passed-into-a-function-with-their-values

#' @returns a list whose components are the symbolic names of the function parameters,
#'  and their values.

#' @export
grabFunctionParameters <- function() {
  pf <- parent.frame()    
  args_names <- ls(envir = pf, all.names = TRUE, sorted = FALSE)
  if("..." %in% args_names) {
    dots <- eval(quote(list(...)), envir = pf)
  }  else {
    dots = list()
  }
  args_names <- sapply(setdiff(args_names, "..."), as.name)
  if(length(args_names)) {
    not_dots <- lapply(args_names, eval, envir = pf) 
  } else {
    not_dots <- list()
  }
  out <- c(not_dots, dots)
  out[names(out) != ""]                                  # remove unnamed values in ... (if any)
}   


#' Sherlock

#' @import stringr
#' @import qpdf
#' @import tableHTML
#' @import dpseg
#' @import plotrix
#' @import zoo
#' @import stargazer
#' @import utils
#' @import graphics
#' @import grDevices
#' @import stats
#' @import textBoxPlacement
#' @import plot.matrix
#' @import devtools

#' @description This function is the driver that organizes
#' the computation of concordances in Sherlock Holmes stories

#' @param titles is a character string containing the full path name
#'  for a text file containing the titles of the stories
#'   in the same order that they appear in the texts file.
#'   If titles=="NONE", treat the entire book as one story.
#' @param texts is a character string containing the full path name for a text file containing the full texts of all of the stories
#' @param patterns is a vector containing the search patterns
#' @param toupper is a Boolean TRUE if the titles should be converted to upper case
#' @param odir is a character string containing the full path name of the output directory
#' @param concord Boolean if TRUE invoke concordance()
#' @param minl is an integer param passed to dpseg::dpseg
#' @param P is a numeric param passed to dpseg::dpseg
#' @param verbose Boolean if TRUE print informative or diagnostic messages to console
#' 
#' @examples
#' titles<-system.file("extdata/contents3.txt",package="SherlockHolmes")
#' texts<-system.file("extdata/processed_download3.txt",package="SherlockHolmes")
#' SH<-Sherlock(titles=titles,texts=texts,patterns=patterns[1],
#'  toupper=TRUE,odir=tempdir(),concord=FALSE,minl=100,P=0.00001,
#'  verbose=FALSE)
#'  
#' 
#' @return returns no value but has side effect of driving the concordance computations
#'
#' @export
Sherlock<-
	function(titles="NONE",texts,patterns,toupper,odir,concord=FALSE,
	         minl=100,P=0.00001,verbose=FALSE) {

	  if(!file.exists(odir))
			dir.create(odir)
	  message(sprintf("OUTPUT DIRECTORY: %s",odir))
		metadata<-sprintf("%s/metadata.txt",odir)
		
		# archive the calling parameters for this run
		# https://stackoverflow.com/questions/11885207/get-all-parameters-as-list
		sink(metadata)
		mc<-match.call()
		print("MATCH CALL:\n")
		print(mc)
		gfp<-grabFunctionParameters()
		print("\nGRAB FUNCTION PARAMETERS:\n")
		print(gfp)
		print(sprintf("OUTPUT DIRECTORY: %s",odir))
		sink()

		texts.vec<-readLines(texts)

		# where does each story start?
		if(titles=="NONE") { # treat the entire book as one story
		  titles.vec<-sprintf("Book Title = [\"%s\"]",basename(texts))
		  starts<-1
		}
		else {
	  	titles.vec<-readTitles(titles)
	  	starts<-startLine(titles.vec,texts.vec,toupper)
		}
		
		# frequency of patterns in each story
		freqs<-frequency(texts.vec,starts,patterns)

		freqDir<-sprintf("%s/freqs/",odir)
		dir.create(freqDir)
		histDir<-sprintf("%s/histograms",freqDir)
		dir.create(histDir)
		chronDir<-sprintf("%s/chronology",freqDir)
		dir.create(chronDir)
		dir.create(sprintf("%s/plots",chronDir))
		dir.create(sprintf("%s/archive",chronDir))
		lengthDir<-sprintf("%s/length",freqDir)
		dir.create(lengthDir)
		dir.create(sprintf("%s/plots",lengthDir))
		dir.create(sprintf("%s/archive",lengthDir))
		
		# histogram of frequencies
		freqHist(patterns,starts,titles.vec,freqs,histDir)
		
		# Plotted in order of date
		# this depends on the titles being given in order of date
		chronology(titles.vec,patterns,starts,freqs,chronDir,overlay=FALSE)
		chronology(titles.vec,patterns,starts,freqs,chronDir,overlay=TRUE)
		coChronology(titles.vec,patterns,starts,freqs,chronDir)

		# Plotted in order of story length
		lengths(titles.vec,patterns,starts,freqs,lengthDir)
		
	  distributions(freqs,titles.vec,minl,P,odir)

		rolling(freqs,titles.vec,windowPct=0.10,odir,verbose)
		
		if(concord)
		  concordance(freqs,titles.vec,texts.vec,starts,window=5,odir)
	}

#' contingency
#' 
#' @description compute chisq value for a 2 x 2 contingency table
#' 
#' @param inside numeric vector of raw counts
#' @param outside numeric vector of raw counts
#' 
#' @examples
#' con<-contingency(inside=c(4,5),outside=c(20,7))
#' 
#' @returns numeric vector of chisq.test() p.values

#' @export
contingency<-
  function(inside,outside) {
    m<-matrix(nrow=2,ncol=2)
    colnames(m)<-c("in","out")
    rownames(m)<-c("target","non_target")
    
    result<-list()
    
    sin<-sum(inside)
    sout<-sum(outside)
    
    for(i in 1:length(inside)) {
      m["target","in"]<-inside[i]
      m["non_target","in"]<-sin-inside[i]
      m["target","out"]<-outside[i]
      m["non_target","out"]<-sout-outside[i]

      w1<-which(is.infinite(m))
      w2<-which(m<0)
      if((length(w1>0)) | (length(w2>0))) {
        print("negative or non-finite input to chisq.test")
        browser(m)
      }
      
      result[i]<-chisq.test(m)$p.value
    }
    
    return(result)	
  }


#' rolling

#' @description compute rolling average of ratio of number of occurrences of query string divided by total number of words

#' @param freqs return value of frequency()
#' @param titles.vec character vector containing the titles of the stories
#' @param windowPct a numeric control size of plot window
#' @param odir character string containing the full path name for the output directory
#' @param verbose Boolean if TRUE print informative or diagnostic messages to console

#' @examples
#' rol<-rolling(freqs,titles.vec,windowPct=0.10,odir=tempdir(),verbose=FALSE)

#' @returns returns noo value, but has side effect of generating graphs
#'
#' @export
rolling<-
  function(freqs,titles.vec,windowPct=0.10,odir,verbose) {
    
    ltv<-length(titles.vec)
    roll.dir<-sprintf("%s/roll.dir",odir)
    dir.create((roll.dir))
    
    # determine number of patterns
    f<-freqs[[titles.vec[1]]]
    npat<-length(names(f$patterns))
    # adjust the size of the pdf accordingly
    scaleFactor<-.5 # extra height in inches per pattern
    height<-7+npat*scaleFactor
    width<-height
    pdf(sprintf("%s/%s.pdf",roll.dir,"compilation"),width=width,height=height)

    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar=c(5,5,6,2)+0.1) # slightly increase top and left margins
    par(cex.lab=1.5) # slightly increase size of axis labels
    
    for(i in 1:ltv) {
      f<-freqs[[titles.vec[i]]]
      
      # compute cumulative sums
      csw<-cumsum(f$wPerLine)

      ymax<-vector("numeric")
      rmsList<-list()
    
      # determine y-axis range for overlay plot
      for(pattern in names(f$patterns)) {
        f1<-f$patterns[[pattern]]
        ppl<-f1$pPerLine
        window<-as.integer(length(ppl)*windowPct)
        # compute rolling average
        rm<-zoo::rollmean(ppl,window)
        if(max(rm)!=min(rm)) { # exclude if just a constant function
          rmsList[[pattern]]<-rm
          ymax[pattern]<-max(rmsList[[pattern]])
        }
      } # for(pattern in names(f$patterns))
    
      # in order to unclutter the overlay plot, we need to stagger the graphs
      # the offset for each graph will be the sum of the max values for all
      # of the preceding graphs.
      # so the stack of staggered graphs will have max y (ie, yymax) equal to
      # the sum of the max's
      ylim<-c(0,sum(ymax))
    
      maxy<-max(ymax)
      sorted<-sort.list(ymax,decreasing=TRUE)
      
      main<-titles.vec[i]
      xlab<-"Lines of Text"
      ylab<-"Rolling Average"

      labs<-list()
      labs$main<-main
      labs$xlab<-xlab
      labs$ylab<-ylab
      positionTextBoxDriverDriverDriver(yList=rmsList,textList=names(rmsList),
              nApprox=10,labs=labs,sortB=TRUE,verbose=verbose)
    } # for(i in 1:ltv)
    
    dev.off()
  }

#' retrieveLmStats

#' @description This function retrieves intercept, slope,
#'  r.squared, and adj.r.squared from lm()

#' @param x is second argument to lm()
#' @param y is first argument to lm()
#' @return returns a list containing the return value of lm,
#' intercept, slope, r.squared, and adj.r.squared
#'
#' @examples
#' retr<-retrieveLmStats(1:10,runif(10,0,1))
#'
#' @export
retrieveLmStats<-
  function(x,y) {
    l<-list()

    lm<-lm(y ~ x)
    
    l$lm<-lm
    # https://stackoverflow.com/questions/66771929/intercept-and-slope-functions-in-r
    cf<-coef(l)
    l$intercept<-cf[1]
    l$slope<-cf[2]
    # https://www.statology.org/extract-r-squared-from-lm-in-r/
    l$r.squared<-summary(lm)$r.squared
    l$adj.r.squared<-summary(lm)$adj.r.squared
    
    return(l)
  }

#' plot_dpseg2

#' @description Alternative plot procedure for dpseg,
#' special function provided personally by dpseg curator.
#' I made a few custom tweeks
#' Including option to overlay multiple plots
#'
#' @import dpseg
#'
#' @param x dpseg object to plot
#' @param delog Boolean use log scale if TRUE
#' @param col color
#' @param main character title of graph
#' @param xlab character label for x axis
#' @param ylab character label for y axis
#' @param res numeric resolution
#' @param vlines Boolean if FALSE suppress vertical lines in graph
#' @param overlay Boolean if TRUE this plot is an overlay of previous plot
#' @param textX numeric x position for text box
#' @param textY numeric y position for text box
#' @param textLabel character string to label the points in the graph
#' @param textLabel character string to label the points in the graph
#' @param ylim numeric vector ylim for plot

#' @return returns no value but has side effect of producing a graph
#' 
#' @examples
#' pdp<-plot_dpseg2(segs,overlay=FALSE,xlab="xaxis",
#'   ylab="yaxis",vlines=FALSE,textX=2000,textY=20,
#'   textLabel="label",ylim=c(0,60))
#' 
#' @export
plot_dpseg2  <- function(x, delog=FALSE, col, main, xlab, ylab, res=10,
                         vlines,overlay,textX,textY,textLabel,ylim) {
  segs <- x$segments
  xy <- x$xy
  
  if ( !"numeric" %in% class(xy$y) ) {
    stop("The results do not contain the original data and can not",
         " be plotted. Use `addLm(dpseg, x=x, y=y)` to add the
original",
         " data and linear regression coefficients.")
  }
  
  if ( delog ) xy$y <- exp(xy$y)
  
  ## colors: generate if not present as argument or column
  if ( missing(col) ) col <- seq_len(nrow(segs))
  if ( !"col"%in%colnames(segs) )
    segs <- cbind(segs,col=col)
  
  cols <- rep(8, length(xy$x))
  for ( i in seq_len(nrow(segs)) ) # note: start overwrites previous end
    cols[segs[i,"start"]:segs[i,"end"]] <- segs[i,"col"]
  
  ## main text: parameters
  if ( missing(main) )
    main <- paste(paste(names(x$parameters),
                        x$parameters,sep=": "),collapse="; ")
  
  if(overlay)
    points(xy, col=cols, main=main)
  else 
    plot(xy, col=cols, main=main, xlab=xlab, ylab=ylab, ylim=ylim)
  
  # https://www.tutorialspoint.com/how-to-rotate-text-in-base-r-plot
  text(textX,textY,textLabel,adj=c(0,-1),srt=+60,cex=2)
  
  ## print predict result as lines
  #  equispaced between min(x) and max(x) at res*length
  xout <- seq(min(xy$x), max(xy$x), length.out=length(xy$x)*res)
  test <- try(pr <- predict(x, xout=xout), silent=TRUE)
  
  if ( !"try-error" %in% class(test) ) {
    if ( delog ) pr$y <- exp(pr$y)
    lines(pr) # TODO: by segment and colored?
  } else warning(test[1])
  
  if ( vlines ) {
    abline(v=xy$x[segs[seq_len(nrow(segs)),"end"]])
    abline(v=xy$x[segs[seq_len(nrow(segs)),"start"]],col=2,lty=2)
  }
  
  invisible(segs)
}

#' segments

#' @description reformat seqs$segments as a legend to insert into segment plot

#' @param segs return value of dpseg::dpseg()
#'
#' @returns reformatted matrix suitable for printing
#'
#' @examples
#' seg<-segments(segs)
#'
#' @export
segments<-
  function(segs) {
    s<-segs$segments
    nrow<-length(s$x1)
    cn<-c("x1","intercept","slope","r2")
    ncol<-length(cn)
    m<-matrix(nrow=nrow,ncol=ncol)
    colnames(m)<-cn
    
    for(r in 1:nrow) {
      m[r,"x1"]<-s$x1[r]
      m[r,"intercept"]<-round(s$intercept[r],2)
      m[r,"slope"]<-round(s$slope[r],4)
      m[r,"r2"]<-round(s$r2[r],2)
    }
    return(m)
  }

#' distributions

#' @description compute distribution of ratio of number of occurrences of query string divided by total number of words

#' @param freqs return value of frequency()
#' @param titles.vec character vector containing the titles of the stories
#' @param minl is an integer param passed to dpseg::dpseg
#' @param P is a numeric param passed to dpseg::dpseg
#' @param odir character string containing the full path name for the output directory
#'
#' @examples
#' dis<-distributions(freqs,titles.vec[1],minl=100,P=0.00001,tempdir())
#'
#' @returns returns no value but has side effect of generating graphs
#'
#' @export	
distributions<-
  function(freqs,titles.vec,minl,P,odir) {
    # file for archiving the lm and segs regression results
    cum.dir<-sprintf("%s/cum.dir",odir)
    dir.create(cum.dir)
    pdf.dir<-sprintf("%s/pdf.dir",cum.dir)
    dir.create(pdf.dir)
    archive.dir<-sprintf("%s/archive.dir",cum.dir)
    dir.create(archive.dir)
   
    ltv<-length(titles.vec)
    pdf(sprintf("%s/%s.pdf",pdf.dir,"compilation"),width=10,height=10)
    for(i in 1:ltv) {
      f<-freqs[[titles.vec[i]]]
      
      # compute cumulative sums
      csw<-cumsum(f$wPerLine)
      
      overlay<-FALSE
      pcount<-0
      nlabels<-length(names(f$patterns))
      xrange<-length(csw)
      
      # pre-compute ylim
      ymax<-0
      for(pattern in names(f$patterns)) {
        f1<-f$patterns[[pattern]]
        csp<-cumsum(f1$pPerLine)
        ymax<-max(ymax,max(csp))
      }
      ylim<-c(0,ymax)
        
      
      # since the curves are more bunched together near 0
      # we will space the labels between xrange*3 and xrange*.9
      del<-xrange*.6/(nlabels-1)
      for(pattern in names(f$patterns)) {
        f1<-f$patterns[[pattern]]
        ofile<-sprintf("%s/cum.txt",archive.dir)
        write(sprintf("CUMULATIVE SUM RESULTS FOR THE TEXT: \"%s\"",titles.vec[i]),ofile,append=TRUE)
        write(sprintf("CUMULATIVE SUM RESULTS FOR THE SEARCH PATTERN: \"%s\"",pattern),ofile,append=TRUE)
        
        # compute cumulative sums
        csp<-cumsum(f1$pPerLine)
        
        # compute text label position
        xVal<-floor(xrange*.3+del*pcount)
        # max val for del*pcount = xrange*.6*(nlabels-1)/(nlabels-1)
        pcount<-pcount+1
        textX<-csw[xVal]
        textY<-csp[xVal]

        # first perform naive linear regression
        lms<-retrieveLmStats(csw, csp)
        lms.tab<-vector("numeric",4)
        lms.tab<-c(lms$lm$coefficients[1],lms$lm$coefficients[2],lms$r.squared,lms$adj.r.squared)
        names(lms.tab)<-c("INTERCEPT","SLOPE","R.SQUARED","ADJUSTED R.SQUARED")
        
        write("\nSIMPLE LINEAR REGRESSION:",ofile,append=TRUE)
        st<-capture.output(stargazer::stargazer(lms.tab,type="text",summary=FALSE))
        cat(paste(st,"\n"),file=ofile,append=TRUE)
        
        # second perform piecewise linear regression
        segs<-dpseg::dpseg(x=csw,y=csp,minl=minl,P=P,verb=0)
        # reformat seqs$segments as a legend to insert into segment plot
        segsm<-segments(segs)
        
        write("\nPIECEWISE LINEAR REGRESSION:",ofile,append=TRUE)
        st<-capture.output(stargazer::stargazer(segsm,type="text",summary=FALSE))
        cat(paste(st,"\n"),file=ofile,append=TRUE)
        
        # plotting
        
        main<-sprintf("%s\nFraction = %f", titles.vec[i],f1$fraction)
        sub<-sprintf("type %s minl %d maxl %d P %f jumps %d",
              segs$parameters$type,
              segs$parameters$minl,segs$parameters$maxl,
              segs$parameters$P,segs$parameters$jumps)
        
        cat(paste(main,"\n"),file=ofile,append=TRUE)
        cat(paste(sub,"\n"),file=ofile,append=TRUE)

        write("\n==============================================\n",ofile,append=TRUE)
        
        xlab<-"Cumulative Total Words"
        ylab<-"Cumulative Pattern"
        
        # special function provided personally by dpseg curator
        
        if(length(f$patterns)==1 | !overlay)
          vlines<-FALSE
        
        if(length(f$patterns)>1) {
          main<-titles.vec[i]
          sub<-""
        }
        
        oldpar <- par(no.readonly = TRUE)
        on.exit(par(oldpar))
        par(mar=c(5,5,6,2)+0.1) # slightly increase top and left margins
        par(cex.lab=1.5) # slightly increase size of axis labels
        plot_dpseg2(segs, vlines=vlines, main=sprintf("%s\n%s",main,sub),
                    xlab=xlab,ylab=ylab,overlay=overlay,textX=textX,textY=textY,
                    textLabel=pattern,ylim=ylim)
        # restore margins to original
        on.exit(par(oldpar))
        
        if(length(f$patterns)==1) {
          # if there is just one search pattern
          # then there is enough room to
          # add legend to segment plot
          # if too long, it is truncated
          # but the complete version is in the archive file
        
          ns<-min(c(8,nrow(segsm))) # too many lines will clutter up the graph 
          ypos<-csp[length(csp)]*.6
          if(nrow(segsm)>1) # generates error when nrow=1
            plotrix::addtable2plot(0,ypos,segsm[1:ns,],hlines=TRUE,vlines=TRUE,xpad=.3)
          
          # prepare legend for naive linear regression results
          # https://stat.ethz.ch/pipermail/r-help/2005-September/079089.html
          labels<-sprintf("SIMPLE LINEAR REGRESSION\n\nINTERCEPT
                          %f\nSLOPE                               %f\nR.SQUARED:                     %f\nADJUSTED R.SQUARED: %f",lms$lm$coefficients[1],lms$lm$coefficients[2],lms$r.squared,lms$adj.r.squared)
          xpos<-max(csw)/2
          ypos<-0
          text(xpos,ypos,adj=c(0,0),labels=labels)
          
        } # if(length(f$patterns)==1)

        # overlay regression line for simple linear regression
        abline(lms$lm)
        overlay<-TRUE
      } # for(pattern in names(f[[pattern]]))
    } # for(i in 1:ltv)
    dev.off()
  }

#' mergeTables
#' 
#' @description merge (inner join) the results in 2 tables generated from 2 vectors
#' 
#' @param tv first table
#' @param tw second table
#' @param cnv character name for column coming from v
#' @param cnw character name for column coming from w
#' 
#' @examples
#' mt<-mergeTables(inside,outside,"in","out")[1:10,]
#' 
#' @returns numeric matrix generated from merging tables from v and w
#' @export
mergeTables<-
  function(tv,tw,cnv,cnw) {
    rn<-setdiff(names(tv),c(""," "))
    m<-matrix(0,nrow=length(rn),ncol=2)
    colnames(m)<-c(cnv,cnw)
    rownames(m)<-rn
    for(n in names(tv)) {
      if(n=="")
        next
      m[n,1]<-tv[n]
      w<-which(names(tw)==n)
      if(length(w)==1)
        m[n,2]<-tw[n]
    }	

    return(m)
  }

#' lengths

#' @description frequencies plotted in order of story length

#' @param titles.vec character vector containing the titles of the stories
#' @param patterns vector of character string query patterns
#' @param starts integer vector of starting positions
#' @param freqs return value of frequency()
#' @param lengthDir character string full path name for output directory
#'
#' @examples
#' freqDir<-tempdir()
#' lengthDir<-sprintf("%s/length",freqDir)
#' dir.create(lengthDir)
#' print(lengthDir)
#' dir.create(sprintf("%s/plots",lengthDir))
#' dir.create(sprintf("%s/archive",lengthDir))
#' le<-lengths(titles.vec,patterns,starts,freqs,lengthDir) 
#'
#' @returns returns no value, but has side effect generating graph
#'
#' @export
lengths<-
  function(titles.vec,patterns,starts,freqs,lengthDir) {
    l<-length(titles.vec)
    z<-1:l
    xlab<-"Story Length"
    ylab<-"Frequency"
    header<-"Frequency Plotted in Story Length Order"
    
    for(pattern in patterns) {
      main<-sprintf("%s\n\"%s\"",header,pattern)
      ofile<-sprintf("%s/archive/%s.txt",lengthDir,pattern)
      v<-vector("numeric",length(starts))
      x<-vector("numeric",length(starts))
      for(i in 1:length(starts)) {
        t<-titles.vec[i]
        v[i]<-freqs[[t]]$patterns[[pattern]]$fraction
        x[i]<-freqs[[t]]$wordSum
        write(sprintf("%s\t%d\t%f",t,x[i],v[i]),ofile,append=TRUE)
      }
      pdf(sprintf("%s/plots/%s_lengths.pdf",lengthDir,pattern))
      plot(x,v,main=main,xlab=xlab,ylab=ylab)
      
      dev.off()
    } # for(pattern in patterns)
  }

#' chronology

#' @description frequencies plotted in order of date
#'  (if the titles are given in order of date)

#' @param titles.vec character vector containing the titles of the stories
#' @param patterns vector of character string query patterns
#' @param starts integer vector of starting positions
#' @param freqs return value of frequency()
#' @param chronDir character string full path name for output directory
#' @param overlay Boolean if TRUE overlay the chronolgy for multiple search patterns
#' 

#' @examples
#' freqDir<-tempdir()
#' chronDir<-sprintf("%s/chronology",freqDir)
#' dir.create(chronDir)
#' dir.create(sprintf("%s/plots",chronDir))
#' dir.create(sprintf("%s/archive",chronDir))
#' print(chronDir)
#' chr<-chronology(titles.vec,c("Holmes","Watson"),starts,freqs,chronDir)
#'
#' @returns returns no value, but has side effect generating graph
#'
#' @export
chronology<-
  function(titles.vec,patterns,starts,freqs,chronDir,overlay=FALSE) {
    l<-length(titles.vec)
    z<-1:l
    xlab<-"Story Number in Chronological Order"
    ylab<-"Frequency"
    header<-"Frequency Plotted in Chronological Order"
    
    first<-TRUE
    if(overlay) {
      COLORS<-c("red","blue","black","orange","green","gray")
      coleurs<-vector("character",length(patterns))
      names(coleurs)<-patterns
      i<-1
      for(pattern in patterns) {
        coleurs[pattern]<-COLORS[i]
        i<-i+1
      } # for(pattern in patterns)
      
      pdf(sprintf("%s/plots/%s_chronology.pdf",chronDir,"composite"))
      mx<-0
      for(pattern in patterns) {
        v<-vector("numeric",length(starts))
        for(i in 1:length(starts)) {
          t<-titles.vec[i]
          v[i]<-freqs[[t]]$patterns[[pattern]]$fraction
        } # for(i in 1:length(starts))
        mx<-max(max(v),mx)
      } # for(pattern in patterns)
    } # if(overlay)
    
    for(pattern in patterns) {
      if(!overlay)
        main<-sprintf("%s\n\"%s\"",header,pattern)
      else
        main<-sprintf("%s",header)
      ofile<-sprintf("%s/archive/%s.txt",chronDir,pattern)
      v<-vector("numeric",length(starts))
      for(i in 1:length(starts)) {
        t<-titles.vec[i]
        v[i]<-freqs[[t]]$patterns[[pattern]]$fraction
        # change to append=FALSE to avoid double-listing
        # when we do both overlay and !overlay
        if(!overlay)
          write(sprintf("%s\t%f",t,v[i]),ofile,append=TRUE)
      }
      
      if(!overlay) {
        pdf(sprintf("%s/plots/%s_chronology.pdf",chronDir,pattern))
        plot(v,main=main,xlab=xlab,ylab=ylab,pch=19)
      }
      
      if(overlay & first) {
        plot(v,main=main,xlab=xlab,ylab=ylab,ylim=c(0,mx),col=coleurs[pattern],pch=19)
        first<-FALSE
      }
      if(overlay & !first)
        points(v,col=coleurs[pattern],pch=19)
      
      if((length(v)>1) & !overlay) {
        lms<-retrieveLmStats(z,v)
        abline(lms$lm)
        
        # prepare legend for naive linear regression results
        # https://stat.ethz.ch/pipermail/r-help/2005-September/079089.html
        labels<-sprintf("SIMPLE LINEAR REGRESSION\n\nINTERCEPT                   %f\nSLOPE                                 %f\nR.SQUARED:                     %f\nADJUSTED R.SQUARED: %f",lms$lm$coefficients[1],lms$lm$coefficients[2],lms$r.squared,lms$adj.r.squared)
        ymax<-max(v)
        text(1,ymax*.8,adj=0,labels=labels)
      } # if(length(v)>1)
      if(!overlay)
        dev.off()
    } # for(pattern in patterns)
    if(overlay) {
      legend("topleft",legend=patterns,fill=coleurs[patterns])
      dev.off()
    }
  }

#' coChronology

#' @description graphical indicator of search patterns within stories

#' @param titles.vec character vector containing the titles of the stories
#' @param patterns vector of character string query patterns
#' @param starts integer vector of starting positions
#' @param freqs return value of frequency()
#' @param chronDir character string full path name for output directory
#' 
#' @examples
#' freqDir<-tempdir()
#' chronDir<-sprintf("%s/chronology",freqDir)
#' dir.create(chronDir)
#' dir.create(sprintf("%s/plots",chronDir))
#' dir.create(sprintf("%s/archive",chronDir))
#' print(chronDir)
#' coch<-coChronology(titles.vec,c("Holmes","Watson"),starts,freqs,chronDir)

#' @returns returns an integer matrix whose rows are search patterns
#'   and columns are stories, value of 1 indicates the presence of the
#'   corresponding search pattern in the corresponding story
#'
#' @export
coChronology<-
  function(titles.vec,patterns,starts,freqs,chronDir) {
    l<-length(titles.vec)
    xlab<-"Story Number in Chronological Order"
    ylab<-"Co-Occurrences of Frequency"
    header<-"Co-Occurrences of Frequency Plotted in Chronological Order"
    
    m<-matrix(0,nrow=length(patterns)+1,ncol=length(starts))
    colnames(m)<-titles.vec
    rownames(m)<-c(patterns,"total")
    
    for(pattern in patterns) {
      main<-sprintf("%s\n\"%s\"",header,pattern)
      ofile<-sprintf("%s/coArchive/%s.txt",chronDir,pattern)
      for(i in 1:length(starts)) {
        t<-titles.vec[i]
        if(freqs[[t]]$patterns[[pattern]]$fraction>0) {
          m[pattern,i]<-1
          m["total",i]<-m["total",i]+1
        } # if(freqs[[t]]$patterns[[pattern]]$fraction>1)
      } # for(i in 1:length(starts))
    } # for(pattern in patterns) 

    # the x axis labels and the title each take up about 5 rows
    # of height for a total of 10 rows
    # so the grand total of rows = length(titles.vec) + 10
    # for Sherlock, there are 60 stories + 10 rows = 70 rows
    # this requires a height of 20 inches, or 20/70 inches per row
    height<-(length(titles.vec) + 10)*20/70
    
    # for Sherlock, when using 5 search patterns, the y axis labels
    # require the same width as the 5 columns for the search patterns
    # we also require 1.5 column widths for the right margin
    # so the grand total of columns = length(patterns) + 6.5
    # for Sherlock, there are 5 search patterns + 6.5 columns = 
    # 11.5 columns
    # this requires a width of 8 inches, or 8/11.5 inches per column
    width<-(length(patterns) + 7.5)*8/11.5
    
    pdf(sprintf("%s/plots/colorChronology.pdf",chronDir),
        width=width,height=height)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar=c(6.1, 20.1, 4.1, 4.1)) # increase bottom and left margins
    # plot() in next line uses the version in package plot.matrix to
    # generate a color map of red and white squares
    plot(t(m[1:(nrow(m)-1),]),col=c("white","red"),axis.row=list(side=2,las=2),
         axis.col=list(side=1,las=2),key=NULL,
         main="At least One Instance of\nSearch Pattern in Text",
         ylab="", xlab="")
    dev.off()
  
    rs<-rowSums(m)         
    return(cbind(m,rs))
  }

#' freqHist

#' @description histogram of frequencies

#' @param patterns vector of character string query patterns
#' @param starts integer vector of starting positions
#' @param titles.vec character vector containing the titles of the stories
#' @param freqs return value of frequency()
#' @param histDir character string full path name for output directory
#'
#' @examples
#' fh<-freqHist(patterns,starts,titles.vec,freqs,histDir=tempdir())
#'
#' @returns returns no value, but has side effect generating histogram
#' @export
freqHist<-
  function(patterns,starts,titles.vec,freqs,histDir) {
    # histogram of frequencies
    for(pattern in patterns) {
      v<-vector("numeric",length(starts))
      for(i in 1:length(starts)) {
        t<-titles.vec[i]
        v[i]<-freqs[[t]]$patterns[[pattern]]$fraction
      } # for(i in 1:length(starts))
      
      pdf(sprintf("%s/%s_histogram.pdf",histDir,pattern))
      # histogram of frequencies
      header<-"Histogram of Frequencies (Patterns/Total Words)"
      main<-sprintf("%s\n\"%s\"",header,pattern)
      xlab<-"Frequency"
      ylab<-"Histogram of Frequencies"
      hist(v,main=main,xlab=xlab,ylab=ylab)
      dev.off()
    } # for(i in 1:length(starts))
  }


#' concordance

#' @description retrieve words that are close to occurrences of pattern

#' @param freqs return value of frequency()
#' @param titles.vec character vector containing the titles of the stories
#' @param texts.vec character vector of entire text
#' @param starts integer vector of starting positions
#' @param window integer number of lines to take before and after the pattern match
#' @param odir character string containing the full path name for the output directory
#'
#' @examples
#' \donttest{
#' con<-concordance(freqs,titles.vec[3],texts.vec,starts,window=2,odir=tempdir())
#' }
#'
#' @returns returns no value but has side effect of generating graphs
#'
#' @export
concordance<-
  function(freqs,titles.vec,texts.vec,starts,window,odir) {
    cdir<-sprintf("%s/concord.dir",odir)
    dir.create(cdir)
    
    ltv<-length(titles.vec)
    # process each title in the text
    arch.dir<-sprintf("%s/arch.dir",cdir)
    dir.create(arch.dir)

    for(i in 1:ltv) {
      # file for archiving the concordance() results
      arch.txt<-sprintf("%s/%s.txt",arch.dir,titles.vec[i])

      start<-starts[i]+1 # add 1 to avoid taking the current title line
      if(i<length(starts))
        end<-starts[i+1]-1 # subtract 1 to avoid taking the next title line
      else
        end<-length(texts.vec)
      
      t<-texts.vec[start:end]
      
      f<-freqs[[titles.vec[i]]]
      
      # process each search pattern
      for(pattern in names(f$patterns)) {
        write(sprintf("CONCORDANCE RESULTS FOR THE TEXT: \"%s\"",titles.vec[i]),arch.txt,append=TRUE)
        write(sprintf("CONCORDANCE RESULTS FOR THE SEARCH PATTERN: \"%s\"",pattern),arch.txt,append=TRUE)
        f1<-f$patterns[[pattern]]
        write(sprintf("window=%d",window),arch.txt,append=TRUE)
        
        pp<-f1$pPerLine
        lpp<-length(pp)
        u<-vector("integer")
        for(j in 1:lpp) {
          if(pp[j]>0) { # is the search pattern found in the jth line of text?
            startW<-max((j-window),1)
            endW<-min((j+window),lpp)
            u<-union(u,(startW:endW))
            write(c("\n[line number in individual story] matching line =",j,"   window =",startW,endW),
                  ncolumns=50,arch.txt,append=TRUE)
            write(c("\n[line number in concatenated text file] matching line =",j+start-1,"   window =",start+startW-1,start+endW-1),
                  ncolumns=50,arch.txt,append=TRUE)
            
            write(c("\n",t[startW:(j-1)]),ncolumns=1,arch.txt,append=TRUE)
            write(c("==>",t[j],"<=="),ncolumns=3,arch.txt,append=TRUE)
            write(c(t[(j+1):endW]),ncolumns=1,arch.txt,append=TRUE)
          } # if(pp[j]>0)
        } # for(j in 1:lpp)
        
        not_u<-setdiff(1:lpp,u)
        # https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html
        #inside<-sort(table(unlist(strsplit(t[u], "\\W+"))),decreasing = TRUE)
        #outside<-sort(table(unlist(strsplit(t[not_u], "\\W+"))),decreasing = TRUE)
        
        inside<-strSplitTab(t[u])
        outside<-strSplitTab(t[not_u])

        mt<-mergeTables(inside,outside,"in","out")
        mt1<-cbind(mt,mt[,"in"]/sum(mt[,"in"]))
        mt2<-cbind(mt1,mt[,"out"]/sum(mt[,"out"]))
        ratio<-round(mt2[,3]/mt2[,4],3)
        mtr<-cbind(mt2,ratio)
        colnames(mtr)<-c("in.raw","out.raw","in.norm","out.norm","ratio.norm")
        
        write("\n",arch.txt,append=TRUE)
        write(sprintf("SEARCH PATTERN: %s\n",pattern),arch.txt,append=TRUE)
        write("COMPARISON OF WORD COUNTS WITHIN VS. OUTSIDE OF THE CONTEXT WINDOW:\n",arch.txt,append=TRUE)
        write("RATIO.NORM IS IN.NORM/OUT.NORM\n",arch.txt,append=TRUE)
        write("CHISQ PVAL IS LESS THAN 0.O5 FOR CANDIDATE CORRELATED WORDS\n",arch.txt,append=TRUE)
        write("AND FOR WHICH RATIO.NORM IS GREATER THAN 1",arch.txt,append=TRUE)
        write("THE MOST LIKELY LEADS (IF ANY) ARE THE TERMS IN THE TOP FOCUS TABLE\n",arch.txt,append=TRUE)
        write("!!BEWARE OF FALSE POSITIVES!!",arch.txt,append=TRUE)
        write("AFTER YOU IDENTIFY A GOOD LEAD, YOU CAN USE THE LEAD TO PERFORM",arch.txt,append=TRUE)
        write("A TEXT SEARCH OF THE CONTEXT PASSAGES LISTED ABOVE",arch.txt,append=TRUE)
        
        w<-which(mtr>0)
        #if((mtr[,"in.raw"]+mtr[,"out.raw"])>0) {
        if(length(w>0)) {
          pval<-contingency(mtr[,"in.raw"],mtr[,"out.raw"])
        
          mtr<-cbind(mtr,pval)
  
          mtr.srt<-mtr[order(unlist(mtr[,"pval"]),decreasing=FALSE),]
        
          w1<-which(mtr.srt[,"ratio.norm"]>1)
          w2<-which(mtr.srt[,"pval"]<=.05)
          w<-intersect(w1,w2)
          
          # https://stackoverflow.com/questions/54113706/how-to-fix-error-in-if-nchartext-matrixr-c-max-lengthreal-c-miss
          # The error is caused because underscores (_) need to be escaped in Latex and that confuses stargazer
          # (even if you're printing text or html).
          rn<-gsub("_","[underscore]",rownames(mtr.srt),fixed=TRUE)
          rownames(mtr.srt)<-rn
          
          if(length(w)>=2) {
            mtr.focus<-mtr.srt[w,]

            st<-capture.output(stargazer::stargazer(mtr.focus,type="text",summary=FALSE))
            cat(paste(st,"\n"),file=arch.txt,append=TRUE)
            write("\n",arch.txt,append=TRUE)
          }
        }
        write("\n==============================================\n",arch.txt,append=TRUE)
        
      } # for(pattern in names(f[[pattern]]))
    } # for(i in 1:ltv)
  }

#' strSplitTab
#' 
#' @description use strsplit to parse words from text t, delete the empty string
#'  from the result, and compile into a sorted table of word frequencies
#' 
#' @param t vector of character strings representing lines of the orginal text
#'  
#' @examples sst<-strSplitTab(texts.vec)
#'  
#' @returns a sorted table of raw word counts
#'  
#' @export
strSplitTab<-
  function(t){
    u<-unlist(strsplit(t, "\\W+"))
    v<-sort(table(u[u!=""]),decreasing = TRUE)
    return(v)
  }

#' frequency

#' @description compute ratio of number of occurrences of query string divided by total number of words

#' @param patterns vector of character string query patterns
#' @param starts integer vector of starting positions
#' @param texts.vec character vector of entire text

#' @examples
#' fr<-frequency(texts.vec,starts,patterns)

#' @returns a list whose components are sub-lists
#'  \itemize{ # indexed by the titles of the stories
#' 	  \item start integer starting line in text
#' 	  \item end integer ending line in text
#' 	  \item wPerLine integer words perline
#' 	  \item wordSum integer sum of wPerLine
#' 	  \item patterns a sub-list 
#'	    \itemize{
#'	      \item integer pPerLine integer patterns per line
#' 		    \item patSum integer total of pPerLine
#' 		    \item fraction numeric ratio of patSum/wordSum
#' 		  }  
#'	}
#'
#' @export
frequency<-
	function(texts.vec,starts,patterns) {
		freqs<-list()
			
		for(i in 1:length(starts)) {
		  title<-names(starts[i])
		  freqs[[title]]<-list()
		 
			start<-starts[i]+1 # add 1 to avoid taking the current title line
			freqs[[title]]$start<-start
			freqs[[title]]$end<-end
			if(i<length(starts))
				end<-starts[i+1]-1 # subtract 1 to avoid taking the next title line
			else
				end<-length(texts.vec)
			
			t<-texts.vec[start:end]
			
			# https://regexone.com/lesson/whitespaces
			# "\\S+" run of one or more non-whitespace characters
			wPerLine<-stringr::str_count(t,"\\S+")
			wordSum<-sum(wPerLine)
			freqs[[title]]$wPerLine<-wPerLine
			freqs[[title]]$wordSum<-wordSum
			freqs[[title]]$patterns<-list()
			
			for(pattern in patterns) {
			  freqs[[title]]$patterns[[pattern]]<-list()
			  pPerLine<-stringr::str_count(t,pattern)
			  patSum<-sum(pPerLine)
			
			  freqs[[title]]$patterns[[pattern]]$pPerLine<-pPerLine
			  freqs[[title]]$patterns[[pattern]]$patSum<-patSum
			  freqs[[title]]$patterns[[pattern]]$fraction<-
			    round(patSum/wordSum,5)
			} # for(pattern in patterns)
		} # for(i in 1:length(starts))

		return(freqs)
	}

#' readTitles

#' @description	read and edit titles to remove blank lines and white space

#' @param titles is a character string containing the full path name for a text file containing the titles of the stories in the same order that thney appear in the texts file

#' @examples
#' titles<-system.file("extdata/contents3.txt",package="SherlockHolmes")
#' rt<-readTitles(titles)

#' @returns a character vector of titles
#' @export
readTitles<-
	function(titles) {
		data<-readLines(titles)
		# https://stackoverflow.com/questions/63851857/removing-elements-with-empty-character-in-string-vector
		data[nzchar(data )]
		
		return(stringr::str_trim(data))
		}

#' startLine

#' @description	where does each story start?

#' @param titles.vec is a character string containing the full path name for a text file containing the titles of the stories in the same order that they appear in the texts file
#' @param texts.vec is a character string containing the full path name for a text file containing the full texts of all of the stories
#' @param toupper is a Boolean TRUE if the titles should be converted to upper case

#' @details each title in titles.vec must appear on a single line
#'    in titles.vec and texts.vec -
#'    a title cannot be split across multiple lines.
#'    each title must only appear one time within titles.vec and texts.vec

#' @examples
#' sl<-startLine(titles.vec,texts.vec,toupper=TRUE)

#' @returns an integer vector of the starting lines of each story
#' @export
startLine<-
  function (titles.vec,texts.vec,toupper) {
    v<-vector("numeric",length(titles.vec))
    names(v)<-stringr::str_trim(titles.vec)
    for(t in names(v)) {
      if(toupper)
        g<-grep(toupper(t),texts.vec,fixed=TRUE,value=FALSE)
      else
        g<-grep(t,texts.vec,fixed=TRUE,value=FALSE)
      
      l<-length(g)
      
      if(l!=1) {
        print("startLine() DETECTED AN ERROR IN texts.vec")
        print(c("THE TITLE SEARCH PATTERN",t,"WAS FOUND",l,"TIMES"))
        print("TYPE Q TO QUIT BROWSER")
        browser()
      } # if(l!=1)
      v[t]<-g
    } # for(t in names(v))
    return(v)
  }

