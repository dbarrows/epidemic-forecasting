## Used to write data out to files for cuIF2
writedata <- function(data, neinum, neibmat, outdir, thread) {

	datafile <- paste(outdir, "/", "data-", thread, ".txt", sep = "")
	write(format(t(data), scientific = FALSE), file = datafile, ncolumns = dim(data)[2])

	neinumfile <- paste(outdir, "/", "neinum-", thread, ".txt", sep = "")
	write(neinum, file = neinumfile, ncolumns = 1)

	neibmatfile <- paste(outdir, "/", "neibmat-", thread, ".txt", sep = "")
	write(t(neibmat), file = neibmatfile, ncolumns = dim(neibmat)[2])

}