library("matconv")
flist <- list.files('GeTallele/toolbox', '*.m$')

newflist <- paste0('R/', gsub('.m', '.r', flist, fixed=TRUE))

sapply(flist, function(x){
    mat2r(inMat=paste0('GeTallele/toolbox/', x), pathOutR =  paste0('R/', gsub('.m', '.r', x, fixed=TRUE)))
})
