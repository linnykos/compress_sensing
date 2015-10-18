tmp = list.files()
files = tmp[grep("loqo_res.*out",tmp)]

for(i in 1:length(files)){
	d = read.table(files[i],sep="\n")
	d2 = as.character(unlist(d))

	x = grep(":*beta0",d2)[1]
	d2 = d2[(x+1):length(d2)]

	d3 = unlist(strsplit(d2,"\\s+"))
	
	if(grepl("sparse",files[i])){
		d3 =  d3[-seq(1,length(d3),4)]
		d3 =  d3[-seq(1,length(d3),3)]
	} else {
		d3 = d3[-seq(1,length(d3),3)]
	}
	d4 = matrix(d3,length(d3)/2,2,byrow=TRUE)
	write.table(d4,file=paste(unlist(strsplit(files[i],".out")),".csv",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep=",")
}
