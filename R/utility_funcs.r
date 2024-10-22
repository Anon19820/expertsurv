
gen_pathway<- function(string){
  folder_bits <- strsplit(string, "\\/")[[1]]
  
  file_extensions <- c("R", "png", "xlsx", "csv", "jpg", "svg", "txt", "pptx","pdf")
  file_extensions_grep <- paste0(paste0("\\.",file_extensions), collapse ="|")
  
  
  exclude.last <- grepl(file_extensions_grep,folder_bits[length(folder_bits)], ignore.case = T)
  
  
  if(exclude.last){
    folder_num <- length(folder_bits) -1
    
  }else{
    folder_num <- length(folder_bits) 
  }
  
  for(i in 1:folder_num){
    folder_path_temp <- paste0(folder_bits[1:i], collapse = "/")
    
    if(!file.exists(folder_path_temp)){
      dir.create(folder_path_temp)
    }
  }
  
  return(string) #Need to return string for actual pasting
  
}
pow <- function(x,y){
  x^y
}