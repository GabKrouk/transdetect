TRANSDETECT.1.0.gb<-function(Target,data1,data2,TF,cor.method='pearson',disable.pval.cutoff=FALSE,threshold.1=1.5,threshold.2=0.75)
{

### makes an output matrix named out that sum up the correlations obtained in the learnt and 
out<-matrix(ncol=4)
Network<-matrix(ncol=2)

#### Find out the TF that potentially bind in the target promoter
print(paste('looking for TF that could bind in the',Target,'promoter'))
cat(paste('looking for TF that could bind in the',Target,'promoter'),file=paste(Target,'Transdetect.rslt.txt'),fill = T)

as.character(Ntwrk.1BS[grep(Target,Ntwrk.1BS[,3]),1])->c

#### Normalize the data
Norm.by.colX.gb(as.matrix(data1))->data1
Norm.by.colX.gb(as.matrix(data2))->data2

proc.time()[1]->To

#### 
C<-list(cor.predict=0,cor.fitted=0)
B<-list(cor.predict=0,cor.fitted=0)

#### potential regulators
row.names(data1)->a
row.names(data2)->b

intersect(intersect(a,b),TF)->d
intersect(d,c)->TF
print(paste('number of potential regulators',length(TF)))
cat(paste('number of potential regulators',length(TF)),file=paste(Target,'Transdetect.rslt.txt'),fill = T,append=TRUE)
#### end potential regulators

### Makes a PDF holding the results
pdf(file=paste(Target,'Transdetect.rslt.pdf'), width=10, height=10)

#### enter in the loop of calculation
for(i in 1:length(TF))
  {
  	print(paste('Processing TF ', i, ' of ', length(TF)))
  	for(j in i:length(TF))
  	{
  		TF[i]->x
  		TF[j]->y
  		
  		if(x==y)
  		{}	
  		else{
  			### this is an external (provided below) script that makes the actual linear fit
  			Connection.tester.gb(Target,x,y,data1,data2,cor.method=cor.method)->A
  			
  			rbind(out,c(A$cor.predict,A$cor.fitted,A$TF1,A$TF2))->out
  			
  			
  			
  			#### rules of the model selection that keep the best model.
  			if((as.numeric(A$cor.predict)+as.numeric(A$cor.fitted)>as.numeric(C$cor.predict)+as.numeric(C$cor.fitted)))
  			{
  				####control the significance of the model based on pvalues of x y (contribution of both TF) and x:y the rule being : (x & y)|(x:y) , see in "connection" script.
  				if(A$significance|disable.pval.cutoff)
  				{
  				
  					
  				
  				
  					### if TRUE replace the model with the better one.
  					A->C
  					
  					#### print the discovery.
  					cat(paste('Better!','sum:',as.numeric(C$cor.predict+C$cor.fitted),'-predicted:',as.numeric(C$cor.predict),'-fitted:',as.numeric(C$cor.fitted),'TF1:',C$TF1,'TF2:',C$TF2,"Significance:",C$significance, 'Time:',as.numeric(proc.time()[1]-To)),file=paste(Target,'Transdetect.rslt.txt'),fill = T,append=TRUE)
  					
  					
  
  				}
  			}
  				
  				#### print attributes of the "good models" even if not "best"
  			
  			if(((as.numeric(A$cor.predict)+as.numeric(A$cor.fitted))>threshold.1)&A$significance&as.numeric(A$cor.fitted)>threshold.2)
  			{
  				A->B
  				
  				cat(paste('Good!','sum:',as.numeric(B$cor.predict+B$cor.fitted),'-predicted:',as.numeric(B$cor.predict),'-fitted:',as.numeric(B$cor.fitted),'TF1:',B$TF1,'TF2:',B$TF2,"Significance:",B$significance, 'Time:',as.numeric(proc.time()[1]-To)),file=paste(Target,'Transdetect.rslt.txt'),fill = T,append=TRUE)
  				
  				#makes the boolean Network
  				rbind(Network,c(B$TF1,paste('AND',i,j)))->Network
  				rbind(Network,c(B$TF2,paste('AND',i,j)))->Network
  				rbind(Network,c(paste('AND',i,j),Target))->Network
  				####
  				
  				
  				#### control plots that are printed in the PDF along the calculation
  				par(mfrow=c(2,2))
  				
  				B$model->model
  				
  				plot(data1[Target,],xlab="learnt data points",ylab=paste(Target,'mRNA level'),main=paste(B$TF1,'*',B$TF2,'-sig:',B$significance),col.main='red',type='b')
  				#lines(data1[B$TF1,],col=colors()[15])
  				#lines(data1[B$TF2,],col=colors()[16])
  				lines(predict(model),col='red')
  				
  				cor(predict(model),data1[Target,],method=cor.method)->c1
  				round(c1,digit=3)->c1
  				plot(predict(model),data1[Target,],xlab='Fit',ylab='observed',main=paste('learnt data set, Cor=',c1))
  				
  				as.numeric(data2[B$TF1,])->x1
  				as.numeric(data2[B$TF2,])->y1
  
  				plot(as.numeric(data2[Target,]),xlab="externat data points",ylab=paste(Target,'mRNA level'),main=paste(B$TF1,'*',B$TF2),col.main='red',type='b')
  				#lines(data2[B$TF1,],col=colors()[15])
  				#lines(data2[B$TF2,],col=colors()[16])
  				lines(predict(model,list(x=x1,y=y1)),col='red')
  				
  				cor(predict(model,list(x=x1,y=y1)),as.numeric(data2[Target,]),method=cor.method)->c2
  				round(c2,digit=3)->c2
  				plot(predict(model,list(x=x1,y=y1)),as.numeric(data2[Target,]),xlab='predicted',ylab='observed',main=paste('external data set, Cor=',c2))
  				##### end of control ploting
  				}
  			}
  		}
  	}

	C->B
		
	#### control plots Final model
	par(mfrow=c(2,2))
	
	B$model->model
	
	plot(data1[Target,],xlab="learnt data points",ylab=paste(Target,'mRNA level'),main=paste(B$TF1,'*',B$TF2,'-sig:',B$significance),col.main='blue',type='b')
	#lines(data1[B$TF1,],col=colors()[15])
	#lines(data1[B$TF2,],col=colors()[16])
	lines(predict(model),col='red')
	
	cor(predict(model),data1[Target,],method=cor.method)->c1
	round(c1,digit=3)->c1
	plot(predict(model),data1[Target,],xlab='Fit',ylab='observed',main=paste('learnt data set, Cor=',c1))
	
	as.numeric(data2[B$TF1,])->x1
	as.numeric(data2[B$TF2,])->y1

	plot(as.numeric(data2[Target,]),xlab="externat data points",ylab=paste(Target,'mRNA level'),main=paste(B$TF1,'*',B$TF2),col.main='blue',type='b')
	#lines(data2[B$TF1,],col=colors()[15])
	#lines(data2[B$TF2,],col=colors()[16])
	lines(predict(model,list(x=x1,y=y1)),col='red')
	
	cor(predict(model,list(x=x1,y=y1)),as.numeric(data2[Target,]),method=cor.method)->c2
	round(c2,digit=3)->c2
	plot(predict(model,list(x=x1,y=y1)),as.numeric(data2[Target,]),xlab='predicted',ylab='observed',main=paste('external data set, Cor=',c2))
	##### end of control ploting final model
		
		
		
		
	proc.time()[1]->Tf
	cat(paste(as.numeric(Tf-To),'seconds of running time'),file=paste(Target,'Transdetect.rslt.txt'),fill = T,append=TRUE)
	print(summary(B$model))
	
	
	
	
	par(mfrow=c(1,1))
	
	
	#####
	dev.off()
	
	if(length(out)>2){
	out[2:length(out[,1]),]->out
	}
	else{}
	
	
	if(length(Network)>2){
		Network[2:length(Network[,1]),]->Network
	}
	else{}
	
		Final<-list(model=C,out=out,Network=Network)
		return(Final)
   }


######################
#######################
############


Norm.by.colX.gb<-function (data, colones) 
{
    col <- length(data[1, ])
    row <- length(data[, 1])
    out <- matrix(nrow = row, ncol = col)
    for (i in 1:length(data[, 1])) {
        out[i, ] <- as.numeric(data[i, ])/mean(as.numeric(data[i, 
            colones]))
    }
    row.names(out) <- row.names(data)
    colnames(out) <- colnames(data)
    return(out)
}

######################
#######################
############

Connection.tester.gb<-function (Target, X, Y, data1, data2, cor.method = "pearson") 
{
    x <- data1[X, ]
    y <- data1[Y, ]
    x1 <- as.numeric(data2[X, ])
    y1 <- as.numeric(data2[Y, ])
    model <- lm(data1[Target, ] ~ x * y)
    corelation <- cor(predict(model), data1[Target, ], method = cor.method)
    prediction <- predict(model, list(x = x1, y = y1))
    corelation2 <- cor(prediction, as.numeric(data2[Target, ]), 
        method = cor.method)
    b1 <- as.data.frame(summary(model)[4])[4]["x", ] < 0.01
    b2 <- as.data.frame(summary(model)[4])[4]["y", ] < 0.01
    b3 <- as.data.frame(summary(model)[4])[4]["x:y", ] < 0.01
    binary <- (b1 & b2) | b3
    The.list <- list(target = Target, TF1 = X, TF2 = Y, model = model, 
        cor.fitted = corelation, cor.predict = corelation2, significance = binary)
    return(The.list)
}

######################
#######################
############

Transdetect.counting.gene.effect.gb<-function(out.data)

## It takes results out of "Transdetect.gb"
## and get the effects of Genes and how many times they are involed in a "good" model.
{
	
	
	
	out.data[[3]]->b
		
	
	union(as.character(b[,1]),as.character(b[,1]))->net.nodes
	net.nodes[grep('At',net.nodes, ignore.case=T)]->genes
	
	##find the target gene
	union(as.character(b[,2]),as.character(b[,2]))->yo
	yo[grep('At',yo, ignore.case=T)]->target
	
	out<-matrix(nrow=length(genes), ncol=2)
	
for(i in 1:length(genes))
{
	
	length(grep(genes[i],as.character(b[,1])))->out[i,2]
	out[i,1]<-genes[i]
	
	
	
}		

	print(paste('your target is',target))
	as.data.frame(out)->out
	out[order(out[,2],decreasing=T),]->out
	colnames(out)<-c('AGI','# of models')
	
	return(out)
	
	
	
	
	
	
}
