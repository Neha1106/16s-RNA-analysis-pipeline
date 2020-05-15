file=as.data.frame(list.files())
dir = as.data.frame(getwd())
file$path = dir$`getwd()`
file$direction = "forward"
colnames(file) = c("sample-id", "absolute-filepath", "direction")
write.csv (file,"/manifist.dir/manifest.csv")