rdsdir <- paste(getwd(), "rds", sep = "/")
filelist <- list.files(rdsdir)

for (filenum in 1:length(filelist)) {
	file <- paste(rdsdir, filelist[filenum], sep = "/")
	vecin <- readRDS(file)
	print(vecin)
}

save.image("fsim.RData")