`%!in%` = Negate(`%in%`)
import_pattern <-function (filepath){
  if (nchar(readChar(filepath, file.info(filepath)$size)) == 
      0) {
    warning(paste(filepath, "imported 0 characters."))
  }
  readChar(filepath, file.info(filepath)$size)
}

extract_pattern <- function (filepath, keyword, preserve = FALSE){
  stopifnot(is.logical(preserve), length(keyword) == 1)
  x <- import_pattern(filepath)
  keyword <- as.character(keyword)
  if (!grepl(keyword, x)) 
    stop("Couldn't find keyword in provided file.")
  locations <- gregexpr(keyword, x)[[1]]
  if (length(locations) < 2) 
    stop("Keyword only appears once in provided file.")
  if (length(locations) > 2) 
    warning("Keyword found in more than two places. \n                                     extract_pattern will only pull text between\n                                     the first two occurances.")
  if (preserve) {
    substr(x, locations[[1]], locations[[2]] + nchar(keyword))
  }
  else {
    substr(x, locations[[1]] + nchar(keyword), locations[[2]] - 
             1)
  }
}
use_parameters <- function (template, names, init.params = TRUE, is.file = FALSE){
  stopifnot(is.logical(init.params), is.logical(is.file), length(template) == 1)
  if (is.file) {
    content <- import_pattern(template)
    x <- extract_pattern(template, "---", preserve = TRUE)
  }
  else {
    content <- template
    locations <- gregexpr("---", content)[[1]]
    if (length(locations) < 2) {
      stop("Not detecting YAML markers in the provided template. Did you mean\n           to use is.file = TRUE?")
    }
    x <- substr(content, locations[[1]], locations[[2]] + 
                  nchar("---"))
  }
  header_length <- nchar(x)
  dots <- names
  new_params <- paste0(vapply(lapply(dots, rlang::as_name), 
                              function(x) paste0("  ", x, ": NA\n"), character(1)), 
                       collapse = "")
  if (init.params) {
    init_params <- paste0(vapply(lapply(dots, rlang::as_name), 
                                 function(x) paste0(x, " <- params$", x, "\n"), character(1)), 
                          collapse = "")
  }
  if (grepl("params:", x)) {
    param_start <- regexpr("params:", x)[[1]][[1]]
    param_end <- regexpr("\n[[:alpha:]]", substr(x, param_start, 
                                                 header_length))[[1]][[1]] + param_start - 1
    if (init.params) {
      paste0(substr(x, 1, param_end), new_params, substr(x, 
                                                         param_end + 1, header_length), "```{r, echo = FALSE}\n", init_params, 
             "```\n", substr(content, header_length, nchar(content)))
    }
    else {
      paste0(substr(x, 1, param_end), new_params, substr(x, 
                                                         param_end + 1, header_length), substr(content, 
                                                                                               header_length, nchar(content)))
    }
  }
  else {
    if (init.params) {
      paste0(substr(x, 1, header_length - 4), "params:\n", 
             new_params, "---\n", "```{r, echo = FALSE}\n\n", init_params, 
             "```\n", substr(content, header_length, nchar(content)))
    }
    else {
      paste0(substr(x, 1, header_length - 4), "params:\n", 
             new_params, "---\n", substr(content, header_length, 
                                         nchar(content)))
    }
  }
}


m_default_gen <- function(){ #Might add arguments to this function
  
  m_default <- matrix(nrow = 3, ncol = 2)
  colnames(m_default) <- c("Cum Prob", "Expert_1")
  #rownames(m_default) <- rep("Expert_1",3)
  m_default[,1] <- c(0.025, 0.5, 0.975)
  return(m_default)
}

m_default_gen2 <- function(){ #Might add arguments to this function
  
  m_default <- matrix(nrow = 1, ncol = 1)
  colnames(m_default) <- c("Expert_1")
  rownames(m_default) <- "MLV"
  return(m_default)
}


return_pooled_info <- function(input_mat, St_indic = 1,dist = "best", mode =NULL){
  #dist_considered <- c("normal","t","gamma", "lognormal", "beta") 
  
  if(St_indic == 1){
    lower_bound = 0
    upper_bound = 1
  }else{
    lower_bound = -Inf
    upper_bound = Inf
  }
  
  
  fit.eval <- fitdist_mod(input_mat[,2:ncol(input_mat), drop = F],
                          probs = input_mat[,1], upper = upper_bound, lower = lower_bound, 
                          expertnames = paste0("Expert_",1:(ncol(input_mat)-1)),
                          mode = mode, trunc = St_indic)
  # browser()
  
  plts_pool <- makePoolPlot(fit= fit.eval,
                            xl =lower_bound,
                            xu =upper_bound,
                            d = dist,
                            w = 1,
                            lwd =1,
                            xlab = "x",
                            ylab =expression(f[X](x)),
                            legend_full = TRUE,
                            ql = NULL,
                            qu = NULL,
                            nx = 200,
                            addquantile = FALSE,
                            fs = 12,
                            expertnames = paste0("Expert_",1:(ncol(input_mat)-1)),
                            St_indic =St_indic)
  
  dfs_pool <-  plts_pool[["data"]]
  if(dist == "best"){
    selc_fit <- fit.eval$best.fitting[,"best.fit"]
  }else{
    selc_fit <- rep(dist, length(fit.eval$best.fitting[,"best.fit"]))
  }
  selc_fit_loc <- sapply(selc_fit, function(x){which(x  == names(fit.eval$ssq))})
  
  pool.df_output <- matrix(nrow = length(selc_fit),ncol = 3)
  colnames(pool.df_output) <- c("param1", "param2", "param3")
  
  for(j in 1:length(selc_fit_loc)){
    pool.df_output[j,1:length(fit.eval[[selc_fit_loc[j]]][j,])] <-  as.numeric(as.vector(fit.eval[[selc_fit_loc[j]]][j,]))
  }
  dfs_expert <- data.frame(dist = names(selc_fit_loc), wi = 1/nrow(pool.df_output), pool.df_output)
  
  return(list(dfs_expert, plts_pool))
}
