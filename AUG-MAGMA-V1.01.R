### OPTIONS ##############################################################################################################################################################

options(scipen = 999)

### PACKAGES #############################################################################################################################################################

if (!(require(data.table))){
  install.packages("data.table")
  library(data.table)
}

if (!(require(stringi))){
  install.packages("stringi")
  library(stringr)
}

if (!(require(igraph))){
  install.packages("igraph")
  library(igraph)
}

if (!(require(foreach))){
  install.packages("foreach")
  library(foreach)
}

if (!(require(parallel))){
  install.packages("parallel")
  library(parallel)
}

if (!(require(doParallel))){
  install.packages("doParallel")
  library(doParallel)
}

### ARGUMENTS ############################################################################################################################################################

command <- commandArgs(trailingOnly = TRUE)

# command <- "--magma /home/elkon/davidgroen/tapuz/tools/magma/magma
# --sumstat /home/elkon/davidgroen/tapuz/input/sumstats/Insomnia_sumstats_Jansenetal.txt
# --sumstat-id 1
# --sumstat-pval 10
# --sumstat-nsample 11
# --binaries /home/elkon/davidgroen/tapuz/input/binaries/g1000_eur
# --baseline-model /home/elkon/davidgroen/tapuz/input/annotations/genes_u0d0.genes.annot
# --augmented-model /home/elkon/davidgroen/tapuz/input/annotations/genes_u35d10.genes.annot
# --gene-set-file /home/elkon/davidgroen/tapuz/input/sets/c5.go.bp.v7.4.entrez.gmt
# --output /home/elkon/davidgroen/tapuz/output_testing_insomnia/"

# --gene-scoring-model top (changes the way gene scores are calculated by MAGMA - this is not recommended; defaults to SNP-Wise Mean)
# --gene-set-format col=A,B (allows a column-format gene-set file, where A is the column number of the column with genes and B is the column number of the column with sets; defaults to assume row-based)
# --ignore-genes <path/to/file> (a list of genes that will be ignored when adjusting gene scores and running gene-set analysis; defaults to none)
# --permutations A (A is the number of permutations; defaults to 20)
# --cores A (A is the number of cores; defaults to a third of available cores)

### FUNCTIONS ############################################################################################################################################################

# 0. general commands

makeDir <- function(path.dir, overwrite = FALSE){
  if (!dir.exists(path.dir)){
    dir.create(path = path.dir, recursive = TRUE)
  } else {
    if (overwrite == TRUE){
      dir.create(path = path.dir, recursive = TRUE)
    }
  }
  return(invisible())
}

# i. command-line parsing

breakCmd <- function(cmd.raw){
  fun.cmd <- stri_split_regex(str = cmd.raw, pattern = "--")
  fun.cmd <- unlist(x = fun.cmd, recursive = FALSE, use.names = FALSE)
  fun.cmd <- stri_trim(str = fun.cmd, side = "both", pattern = "\\P{Wspace}")
  return(fun.cmd[-1])
}

getArgs <- function(cmd.parts){
  fun.args <- stri_replace_all_regex(str = cmd.parts, pattern = "^([:alpha:]+[-]*[:alpha:]{0})*\\s*", replacement = "", vectorize_all = TRUE)
  fun.args <- stri_replace_all_regex(str = fun.args, pattern = "\\s*$", replacement = "", vectorize_all = TRUE)
  fun.path.sums <- getFileName(path.file = cmd.parts[which(startsWith(x = cmd.parts, prefix = "sumstat "))])
  fun.path.aug <- getFileName(path.file = cmd.parts[which(startsWith(x = cmd.parts, prefix = "augmented-model "))])
  return(c(fun.args, fun.path.sums, fun.path.aug))
}

getFileName <- function(path.file){
  fun.path <- stri_replace_first_regex(str = path.file, pattern = ".*\\/.*\\/", replacement = "")
  fun.path <- stri_replace_all_regex(str = fun.path, pattern = "\\..*$", replacement = "")
  return(fun.path)
}

getFl <- function(cmd.parts){
  fun.flags <- stri_extract_first_regex(str = cmd.parts, pattern = "^([:alpha:]+[-]*[:alpha:]{0})*")
  return(c(fun.flags, "study", "augmentation"))
}

makeCmd <- function(cmd.raw){
  fun.cmd <- breakCmd(cmd.raw = cmd.raw)
  return(setNames(object = getArgs(cmd.parts = fun.cmd), nm = getFl(cmd.parts = fun.cmd)))
}

rmEmptyFl <- function(cmd.made){
  return(cmd.made[which(cmd.made != "" | !is.na(cmd.made))])
}

# ii. command-line argument validation

checkAugModelFl <- function(cmd.ready, flag = "augmented-model"){
  if (!(flag %in% names(cmd.ready)) || !file.exists(cmd.ready[flag])){
    stop(paste0("<< File containing augmented-model gene annotations not found (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " file (--", flag, "): ", cmd.ready[flag])
  }
  return(invisible())
} 

checkBslModelFl <- function(cmd.ready, flag = "baseline-model"){
  if (!(flag %in% names(cmd.ready)) || !file.exists(cmd.ready[flag])){
    stop(paste0("<< File containing baseline-model gene annotations not found (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " file (--", flag, "): ", cmd.ready[flag])
  }
  return(invisible())
}

checkBinFl <- function(cmd.ready, flag = "binaries"){
  if (!(flag %in% names(cmd.ready)) || (FALSE %in% file.exists(paste0(cmd.ready[flag], c(".bim", ".bed", ".fam"))))){
    stop(paste0("<< At least one of the binary files (.bim, .bam or .fam) is missing (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " (--", flag, "): ", cmd.ready[flag], ".(bim|bed|fam)")
  }
  return(invisible())
}

checkCoresFl <- function(cmd.ready, flag = "cores"){
  if (!(flag %in% names(cmd.ready))){
    fun.cores <- setNames(object = c(as.character(as.integer(detectCores()/4))), nm = "cores")
    message(">> Number of ", toupper(flag), " not specified (--", flag, "). Defaulting to ", fun.cores, " cores. Continuing.")
    return(fun.cores)
  } else {
    if (as.integer(cmd.ready[flag]) <= 0 || is.na(as.integer(cmd.ready[flag]))){
      stop(paste0("<< Number of ", toupper(flag), " specified is invalid (--", flag, "). Terminating."))
    } else if (as.integer(cmd.ready[flag]) > 35){
      message(">> Number of ", toupper(flag), " requested is large, which may result in unexpected crashing", " (--", flag, "): ", cmd.ready[flag])
    } else {
      message(">> Number of ", toupper(flag), " requested (--", flag, "): ", cmd.ready[flag])
    }
    return(c())
  }
}

checkMagmaFl <- function(cmd.ready, flag = "magma"){
  if (!(flag %in% names(cmd.ready)) || !file.exists(cmd.ready[flag])){
    stop(paste0("<< Path to MAGMA software is invalid (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " executable (--", flag, "): ", cmd.ready[flag])
  }
  return(invisible())
}

checkOutFl <- function(cmd.ready, flag = "output"){
  if (!(flag %in% names(cmd.ready)) || !dir.exists(cmd.ready[flag])){
    stop(paste0("<< ", toupper(flag), " directory does not exist (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " file (--", flag, "): ", cmd.ready[flag])
  }
  return(invisible())
}

checkPermsFl <- function(cmd.ready, flag = "permutations"){
  if (!(flag %in% names(cmd.ready))){
    fun.perms <- setNames(object = c("20"), nm = "permutations")
    message(">> Number of ", toupper(flag), " not specified (--", flag, "). Defaulting to 20 permutations. Continuing.")
    return(fun.perms)
  } else {
    if (as.integer(cmd.ready[flag])<=0 || is.na(as.integer(cmd.ready[flag]))){
      stop(paste0("<< Number of ", toupper(flag), " specified is invalid (--", flag, "). Terminating."))
    } else if (as.integer(cmd.ready[flag]) > 50){
      message(">> Number of ", toupper(flag), " requested is large, which may result in long run times", " (--", flag, "): ", cmd.ready[flag])
    } else {
      message(">> Number of ", toupper(flag), " requested (--", flag, "): ", cmd.ready[flag])
    }
    return(c())
  }
}

checkRmGenesFl <- function(cmd.ready, flag = "ignore-genes"){
  if (!(flag %in% names(cmd.ready)) || !file.exists(cmd.ready[flag])){
    message(">> ", toupper(flag), " file (genes excluded from gene-score adjustment and gene-set analysis) not specified or found (--", flag, "). Continuing.")
    return(setNames(object = c(""), nm = "ignore-genes"))
  } else {
    message(">> ", toupper(flag), " file (--", flag, "): ", cmd.ready[flag])
    return(c(""))
  }
}

checkSetFileFl <- function(cmd.ready, flag = "gene-set-file"){
  if (!(flag %in% names(cmd.ready)) || !file.exists(cmd.ready[flag])){
    stop(paste0("<< File containing gene sets not found (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " file (--", flag, "): ", cmd.ready[flag])
  }
  return(invisible())
}

checkSetFormFl <- function(cmd.ready, flag = "gene-set-format"){
  if (!(flag %in% names(cmd.ready)) || !(stri_detect_regex(str = cmd.ready[flag], pattern = "^col=[0-9]+,[0-9]+$"))){
    message(">> ", toupper(flag), " not specified or specified incorrectly (--", flag, "). Assuming row-based format. Continuing.")
    return(setNames(object = c(""), nm = "gene-set-format"))
  } else {
    fun.info <- sub(pattern = "^col=", replacement = "", x = cmd.ready[flag]) 
    fun.info <- stri_split_regex(str = fun.info, pattern = ",")
    fun.info <- unlist(x = fun.info, use.names = FALSE, recursive = FALSE)
    message(">> ", toupper(flag), " (--", flag, "): ", fun.info[1], " (genes) and ", fun.info[2], " (sets).")
    return(c(""))
  }
}

checkScoringFl <- function(cmd.ready, flag = "gene-scoring-model"){
  if (!(flag %in% names(cmd.ready)) || cmd.ready[flag] != "top"){
    message(">> ", toupper(flag), " not specified or recognized (--", flag, "). Defaulting to SNP-Wise Mean. Continuing.")
    return(setNames(object = c(""), nm = "gene-scoring-model"))
  } else {
    message(">> ", toupper(flag), " (--", flag, "): ", cmd.ready[flag])
    return(c(""))
  }
}

checkSumsFl <- function(cmd.ready, flag = "sumstat"){
  if (!(flag %in% names(cmd.ready)) || !file.exists(cmd.ready[flag])){
    stop(paste0("<< File containing summary statistics not found (--", flag, "). Terminating."))
  } else {
    message(">> ", toupper(flag), " file (--", flag, "): ", cmd.ready[flag])
  }
  return(invisible())
}

checkSumsColFl <- function(cmd.ready, flags = c("sumstat-id", "sumstat-pval", "sumstat-nsample")){
  if (FALSE %in% (flags %in% names(cmd.ready))){
    stop(paste0("<< At least one column index for the summary statistics was not specified (--sumstat-id, --sumstat-pval, --sumstat-nsample). Terminating."))
  } else {
    for (e in paste0(">> ", toupper(flags), " column (--", flags , "): ", cmd.ready[flags])){
      message(e)
    }
  }
  return(invisible())
}

checkCmd <- function(cmd.ready){
  fun.cmd <- cmd.ready
  checkMagmaFl(cmd.ready = fun.cmd)
  checkSumsFl(cmd.ready = fun.cmd)
  checkSumsColFl(cmd.ready = fun.cmd)
  checkBinFl(cmd.ready = fun.cmd)
  checkBslModelFl(cmd.ready = fun.cmd)
  checkAugModelFl(cmd.ready = fun.cmd)
  checkSetFileFl(cmd.ready = fun.cmd)
  checkOutFl(cmd.ready = fun.cmd)
  return(c(fun.cmd, 
           checkScoringFl(cmd.ready = fun.cmd), 
           checkSetFormFl(cmd.ready = fun.cmd), 
           checkRmGenesFl(cmd.ready = fun.cmd), 
           checkPermsFl(cmd.ready = fun.cmd), 
           checkCoresFl(cmd.ready = fun.cmd)))
}

# iii. preprocessing summary statistics

readSums <- function(path.sums, col.id, col.p, col.n){
  message("Reading-in summary statistics for processing...")
  fun.sums <- fread(file = path.sums, header = TRUE, fill = TRUE, select = as.integer(c(col.id, col.p, col.n)), col.names = c("id", "p", "n"), data.table = TRUE)
  fun.sums[, id:=as.character(id)][, p:=as.numeric(p)][, n:=as.integer(n)]
  message("Detected ", as.character(nrow(fun.sums)), " entries in summary statistics.")
  return(fun.sums)
}

readBim <- function(path.bins, col.id = 2, col.chr = 1, col.pos = 4){
  message("Reading-in reference population .bim file to annotate variants in summary statistics with coordinates...")
  fun.bim <- fread(file = paste0(path.bins, ".bim"), header = FALSE, fill = TRUE, select = as.integer(c(col.id, col.chr, col.pos)), col.names = c("id", "chr", "pos"), data.table = TRUE)
  fun.bim[, id:=as.character(id)][, chr:=as.character(chr)][, pos:=as.integer(pos)]
  message("Detected ", as.character(nrow(fun.bim)), " entries in reference population .bim file.")
  return(fun.bim)
}

prepSums <- function(data.sums, data.bim, path.main, name.sums, name.aug, type.sums.out = "original", prefix.sums.out = "sumstat", count.perms){
  fun.bim  <- copy(data.bim)
  setkey(fun.bim, id)
  fun.sums <- copy(data.sums)
  setkey(fun.sums, id)
  message("Matching variants in summary statistics to those in .bim file (by rs-id) for annotation...")
  fun.sums <- merge(x = fun.sums, y = fun.bim, by = "id")
  message("After matching, ", as.character(nrow(fun.sums)), " entries remain in summary statistics.")
  message("Applying some quality-control steps and filtering to summary statistics...")
  fun.sums[, chr:=prepChr(chr)][getDupId(fun.sums[, id]), id:=NA][getDupCoord(data.chr = fun.sums[, chr], data.pos = fun.sums[, pos]), pos:=NA]
  fun.sums[p < (1*(10^(-323))), p:=(1*(10^(-323)))][p > 0.9999999999999999, p:=0.9999999999999999]
  fun.sums[pos < 1, pos:=NA][chr == 6 & pos > 28477796 & pos < 33448355, pos:=NA][n < 50, n:=NA]
  fun.sums[chr %in% as.integer(fun.sums[, .N, by = chr][N < (as.integer(count.perms) + 1), chr]), chr:=NA]
  fun.sums <- na.omit(fun.sums)[order(chr, pos), .(chr, pos, id, p, n)]
  checkSums(data.sums = fun.sums, count.perms = count.perms) 
  writeSums(data.sums = fun.sums, path.main = path.main, name.sums = name.sums, name.aug = name.aug, type.sums.out = type.sums.out, prefix.sums.out = prefix.sums.out)
  return(fun.sums[, newp := as.numeric(0)][, newn := as.integer(0)][, .(chr, id, p, n, newp, newn)])
}

prepChr <- function(data.chr){
  fun.chr <- stri_replace_all(str = data.chr, replacement = "", regex = "(c|C)+(h|H)*(r|R)*")
  fun.chr <- as.integer(stri_replace_all(str = fun.chr, replacement = "23", regex = "(x|X)"))
  fun.chr <- ifelse(test = fun.chr %in% 1:23, yes = fun.chr, no = NA)
  return(fun.chr)
}

getDupId <- function(data.id){
  fun.dup <- data.id[which(duplicated(x = data.id))]
  return(which(data.id %in% fun.dup))
}

getDupCoord <- function(data.chr, data.pos){
  fun.chr.pos <- paste0(as.character(data.chr), "_", as.character(data.pos))
  fun.dup <- fun.chr.pos[which(duplicated(x = fun.chr.pos))]
  return(which(fun.chr.pos %in% fun.dup))
}

writeSums <- function(data.sums, path.main, name.sums, name.aug, type.sums.out, prefix.sums.out, log.verbose = "on"){
  fun.path <- paste0(path.main, name.sums, "/", name.aug, "/sumstat/", type.sums.out, "/")
  makeDir(path.dir = fun.path)
  fwrite(x = data.sums, file = paste0(fun.path, prefix.sums.out, ".txt"), sep = "\t", col.names = TRUE, quote = FALSE)
  if (log.verbose == "on"){
    message("Writing cleaned-up summary statistics for ", as.character(nrow(data.sums)), " variants to ", fun.path, prefix.sums.out, ".txt")
  }
  return(invisible())
}

checkSums <- function(data.sums, count.perms){
  if (nrow(data.sums) == 0){
    stop("<< No variants in summary statistics after matching and filtering. Terminating.")
  } else {
    fun.chr <- sort(unique(data.sums[, chr]))
    for (j in fun.chr){
      if (nrow(data.sums[chr == j]) > as.integer(count.perms)){
        message(">> Summary statistics are suitable for EPVP. Continuing.")
        break()
      } else {
        if (j == max(fun.chr)){
          stop("<< Not one chromosome contains sufficient entries for the number of permutations specified. Terminating.")
        }
      }
    }
    if (nrow(data.sums) < 3000000){
      warning("Fewer than 3 million variants remain after matching and filtering. Preprocessing of summary statistics may be required.")
    }
    if (FALSE %in% (1:22 %in% fun.chr)){
      fun.missing <- paste(c(1:22)[which(!(1:22 %in% fun.chr))], collapse = ", ")
      warning("Not all autosomes (", fun.missing, ") are present in summary statistics after matching and filtering. Preprocessing of the summary statistics may be required.")
    }
  }
  return(invisible())
}

# iv. preparing variant-to-gene annotation table

prepAnnotTbl <- function(data.annot.bsl, data.annot.aug, data.sums, path.main, name.sums, name.aug){
  message("Preparing annotation table...")
  fun.annot.aug <- copy(data.annot.aug)
  setkey(fun.annot.aug[, type := "a"], gene, id)
  fun.annot.bsl <- copy(data.annot.bsl)
  setkey(fun.annot.bsl[, type := "b"], gene, id)
  fun.annot <- merge(x = fun.annot.aug, y = fun.annot.bsl, by = c("gene", "id"), all = TRUE)
  if ((nrow(fun.annot) == 0) || (length(intersect(data.sums[, id], fun.annot.bsl[, id])) == 0)){
    stop("<< There is not one variant in the annotation table that also occurs in the summary statistics. Terminating.")
  } else {
    if (nrow(fun.annot) != nrow(fun.annot.aug)){
      stop("<< Issue detected after merging annotation files. Baseline model is probably not a subset of augmented model. Terminating.")
    } else {
      fun.annot[!is.na(type.y), type.x := "b"][, type.y := NULL]
      setnames(x = setcolorder(x = fun.annot, neworder = c("gene", "type.x", "id")), old = c("gene", "type.x", "id"), new = c("gene", "type", "id"))
      if (!("a" %in% fun.annot[, type])){
        stop("<< Nothing to permute. All variants are part of baseline model. Terminating.")
      } else {
        writeAnnotTbl(data.annot.tbl = fun.annot, path.main = path.main, name.sums = name.sums, name.aug = name.aug)
      }
    }
  }
  return(fun.annot)
}

readAnnot <- function(path.annot){
  fun.annot <- readLines(con = path.annot)
  return(fun.annot[which(!startsWith(x = fun.annot, prefix = "#"))])
}

reformatAnnot <- function(data.annot){
  fun.annot.split <- stri_split_fixed(str = data.annot, pattern = "\t", simplify = FALSE)
  fun.genes <- unlist(x = lapply(X = fun.annot.split, FUN = function(x){x = x[1]}), use.names = FALSE, recursive = FALSE)
  if (TRUE %in% duplicated(fun.genes)){
    stop("<< After reading-in annotation file, detected the same gene symbol on multiple rows. Terminating.")
  } else {
    fun.annot.whole <- matrix(nrow = sum(sapply(X = fun.annot.split, FUN = function(x){length(x)-2})), ncol = 3)
    fun.annot.whole <- setNames(object = as.data.table(x = fun.annot.whole, stringsAsFactors = FALSE), nm = c("gene", "type", "id"))
    fun.annot.whole[, c("gene", "type", "id"):=lapply(.SD, as.character), .SDcols = c("gene", "type", "id")]
    fun.annot.whole[, gene:=unlist(x = lapply(X = fun.annot.split, FUN = function(x){rep(x = x[1], length(x)-2)}), recursive = FALSE, use.names = FALSE)]
    fun.annot.whole[, id:=unlist(x = lapply(X = fun.annot.split, FUN = function(x){tail(x = x, -2)}), recursive = FALSE, use.names = FALSE)]
    fun.annot.whole <- unique(x = fun.annot.whole)
  }
  return(fun.annot.whole)
}

writeAnnotTbl <- function(data.annot.tbl, path.main, name.sums, name.aug){
  fun.path <- paste0(path.main, name.sums, "/", name.aug, "/annotation/")
  makeDir(fun.path)
  fun.rows.count <- nrow(data.annot.tbl)
  fun.genes.count <- length(unique(data.annot.tbl[, gene]))
  fun.id.count <- length(unique(data.annot.tbl[, id]))
  fun.a.count <- length(data.annot.tbl[type == "a", type])
  fun.b.count <- length(data.annot.tbl[type == "b", type])
  message("Completed annotation table. Writing a copy to ", fun.path, "annotation.table for reference.")
  message("Table contains ", as.character(fun.rows.count), " variant-to-gene annotations for ", as.character(fun.genes.count), " unique genes and ", as.character(fun.id.count), " unique variants.")
  message("Across all entries, ", as.character(fun.b.count), " are labelled as type b (part of baseline) and ", as.character(fun.a.count), " as type a (unique to augmentation)")
  fwrite(x = data.annot.tbl, file = paste0(fun.path, "annotation.table"), sep = "\t", col.names = TRUE, quote = FALSE)
  return(invisible())
}

# v. gene-scoring and gene-set analysis

delaySetsAnalysis <- function(path.file.in, timer = 1000000){
  stopifnot(file.exists(path.file.in))
  fun.size.start <- file.size(path.file.in)
  for (i in 1:timer){
    fun.size.end <- file.size(path.file.in)
    if (i == timer){
      fun.size.test <- fun.size.end > fun.size.start
      if (fun.size.test == TRUE){
        delayProgram(path.file.in = path.file.in, timer = 1000000)
      }
    }
  }
  return(invisible())
}

scoreGenes <- function(path.magma, path.bins, path.annot, path.main, name.sums, name.aug, type.sums.in, prefix.sums.in, type.score.out, prefix.score.out, score.mdl, genes.only = "no"){
  fun.path.sums  <- paste0(path.main, name.sums, "/", name.aug, "/sumstat/", type.sums.in, "/", prefix.sums.in, ".txt")
  fun.path.score <- paste0(path.main, name.sums, "/", name.aug, "/scores/", type.score.out, "/")
  makeDir(fun.path.score)
  fun.score.mdl  <- ifelse(test = (score.mdl  != "top"), yes = "", no = " --gene-model snp-wise=top")
  fun.genes.only <- ifelse(test = (genes.only != "yes"), yes = "", no = " --genes-only")
  fun.cmd <- paste0(path.magma, " --bfile ", path.bins, " --pval ", fun.path.sums, " use=id,p ncol=n --gene-annot ", path.annot, fun.score.mdl, fun.genes.only, " --out ", fun.path.score, prefix.score.out)
  message("Executing gene-scoring for ", type.score.out, " model...")
  try(expr = system(command = fun.cmd, ignore.stdout = TRUE))
  if (genes.only == "no"){delaySetsAnalysis(path.file.in = paste0(fun.path.score, prefix.score.out, ".genes.raw"))}
  message("Completed gene-scoring for ", type.score.out, " model. All output is located in the ", fun.path.score, " directory.")
  return(invisible())
}

scoreGenesPar <- function(path.magma, path.bins, path.main, name.sums, name.aug, type.sums.in, type.score.out, score.mdl, genes.only = "yes", count.cores){
  fun.path <- paste0(path.main, name.sums, "/", name.aug, "/sumstat/", type.sums.in, "/")
  fun.run.sums  <- list.files(path = fun.path, pattern = "^C-[0-9]+_P")
  fun.run.sums  <- sub(pattern = "\\.txt$", replacement = "", x = fun.run.sums)
  fun.run.annot <- paste0(stri_extract_first(str = fun.run.sums, regex = "^C-[0-9]+"), ".genes.annot")
  fun.run.info  <- cbind.data.frame(annot = fun.run.annot, sumstat = fun.run.sums, stringsAsFactors = FALSE)
  fun.cores <- makeForkCluster(nnodes = as.integer(count.cores))
  registerDoParallel(fun.cores)
  on.exit(stopCluster(fun.cores))
  fun.run.time  <- ceiling((1.5*(((nrow(fun.run.info)/as.integer(count.cores))*35)/60)/60))
  fun.run.time  <- ifelse(test = fun.run.time == 1, yes = "1 hour.", no = paste0(as.character(fun.run.time), " hours."))
  message("Parallel execution using ", as.character(count.cores), " cores. Procedure is estimated to take ", fun.run.time)
  foreach(r = 1:nrow(fun.run.info), .inorder = FALSE) %dopar%
    scoreGenes(path.magma       = path.magma,
               path.bins        = path.bins,
               path.annot       = paste0(path.main, name.sums, "/", name.aug, "/sumstat/", type.sums.in, "/", fun.run.info[r, "annot"]),
               path.main        = path.main,
               name.sums        = name.sums,
               name.aug         = name.aug, 
               type.sums.in     = type.sums.in,
               prefix.sums.in   = fun.run.info[r, "sumstat"],
               type.score.out   = type.score.out,
               prefix.score.out = fun.run.info[r, "sumstat"],
               score.mdl        = score.mdl,
               genes.only       = genes.only)
  return(invisible())
}

scoreSets <- function(path.magma, path.sets, path.main, name.sums, name.aug, type.score, prefix.score.in, prefix.score.out, format.sets, gene.exclude){
  fun.path <- paste0(path.main, name.sums, "/", name.aug, "/scores/", type.score, "/")
  fun.format.sets  <- ifelse(test = (format.sets   != ""), yes = paste0(" ", format.sets), no = "")
  fun.gene.exclude <- ifelse(test = (gene.exclude != ""), yes = paste0(" gene-exclude=", gene.exclude), no = "")
  fun.cmd <- paste0(path.magma, " --gene-results ", fun.path, prefix.score.in, ".genes.raw --set-annot ", path.sets, fun.format.sets, " --settings gene-info ", fun.gene.exclude, " --out ", paste0(fun.path, prefix.score.out))
  message("Executing gene-set analysis for ", type.score, " model...")
  try(expr = system(command = fun.cmd, ignore.stdout = TRUE))
  message("Completed gene-set analysis for ", type.score, " model. All output is located in the ", fun.path, " directory.")
  return(invisible())
}

scoreSetsPar <- function(path.magma, path.sets, path.main, name.sums, name.aug, type.score, format.sets, gene.exclude, count.cores){
  fun.run.info <- list.files(path = paste0(path.main, name.sums, "/", name.aug, "/scores/", type.score, "/"), pattern = ".*.genes.raw$")
  fun.run.info <- sub(pattern = "\\.genes\\.raw$", replacement = "", x = fun.run.info)
  fun.cores <- makeForkCluster(nnodes = as.integer(count.cores))
  registerDoParallel(fun.cores)
  on.exit(stopCluster(fun.cores))
  print(paste0("Parallel execution using ", as.character(count.cores), " cores. This won't take long."))
  foreach(r = 1:length(fun.run.info)) %dopar%
    scoreSets(path.magma       = path.magma,
              path.sets        = path.sets,
              path.main        = path.main,
              name.sums        = name.sums,
              name.aug         = name.aug, 
              type.score       = type.score,
              prefix.score.in  = fun.run.info[r],
              prefix.score.out = paste0("scores.permutation-", r, ".adjusted"),
              format.sets      = format.sets,
              gene.exclude     = gene.exclude)
  return(invisible())
}

# vi. permutation of summary statistics

permuteSums <- function(data.sums, data.annot.tbl, count.perms, path.main, name.sums, name.aug, type.sums.out = "permuted", count.cores){
  fun.track <- c()
  fun.progr <- c(". Current chromosome: ")
  fun.annot.tbl <- copy(data.annot.tbl)
  setkey(fun.annot.tbl, id)
  fun.sums.full <- copy(data.sums)
  fun.chr   <- sort(unique(fun.sums.full[, chr]))
  fun.perms <- as.integer(count.perms)
  fun.cores <- makeForkCluster(nnodes = as.integer(count.cores))
  registerDoParallel(fun.cores)
  on.exit(stopCluster(fun.cores))
  message("Permuting summary statistics (chromosome-by-chromosome for ", as.character(length(fun.chr)), " chromosomes). Requested ", count.perms, " permutations.")
  for (j in fun.chr){
    if (TRUE %in% duplicated(fun.track)){
      stop("<< At least one gene is annotated with variants from multiple chromosomes. This is incompatible with EPVP. Terminating.")
    } else {
      fun.count <- 0
      fun.progr <- paste0(fun.progr, as.character(j), "...")
      message(Sys.time(), fun.progr)
      fun.sums.chr <- fun.sums.full[chr == j]
      fun.shifts   <- sample(x = 1:(nrow(fun.sums.chr)-1), size = fun.perms)
      for (i in fun.shifts){
        fun.count <- fun.count + 1
        fun.sums  <- copy(fun.sums.chr)
        fun.rows  <- nrow(fun.sums)
        fun.sums[1:i, newp:=tail(x = fun.sums[, p], n = i)][1:i, newn:=tail(x = fun.sums[, n], n = i)]
        fun.sums[(i + 1):fun.rows, newp:=head(x = fun.sums[, p], n = (-1 * i))][(i + 1):fun.rows, newn:=head(x = fun.sums[, n], n = (-1 * i))]
        setkey(fun.sums, id)
        fun.sums  <- merge(x = fun.sums, y = fun.annot.tbl, by = "id")
        setcolorder(x = fun.sums[type == "a", p:=newp][type == "a", n:=newn][, chr:=NULL][, type:=NULL][, newp:=NULL][, newn:=NULL], neworder = c("gene", "id", "p", "n"))
        fun.genes <- unique(fun.sums[, gene])
        foreach(h=1:length(fun.genes), .inorder = FALSE) %dopar% 
          writeSums(data.sums = fun.sums[gene == fun.genes[h]], path.main = path.main, name.sums = name.sums, name.aug = name.aug, type.sums.out = type.sums.out,  prefix.sums.out = paste0(fun.genes[h], "_", as.character(fun.count)), log.verbose = "off")
        if (fun.count == fun.perms){
          fun.track <- c(fun.track, fun.genes)
        }
      }
    }
  }
  message("Completed all permutations. All files were written to the ", path.main, name.sums, "/", name.aug, "/sumstat/", type.sums.out, "/ directory.")
  return(invisible())
}

checkPermSums <- function(path.main, name.sums, name.aug, type.score.in = "augmented", prefix.score.in, type.sums.in = "permuted", count.perms){
  message("Checking for the existence of a summary statistics file for each gene and each permutation...")
  fun.path <- paste0(path.main, name.sums, "/", name.aug)
  fun.genes.orig <- extractRawCol(data.raw = readRaw(path.raw = paste0(fun.path, "/scores/", type.score.in, "/", prefix.score.in, ".genes.raw")), index = 1)
  for (j in 1:as.integer(count.perms)){
    fun.files <- list.files(path = paste0(fun.path, "/sumstat/", type.sums.in, "/"), pattern = paste0("_", as.character(j), ".txt$"))
    fun.files <- stri_replace(str = fun.files, replacement = "", regex = "_([0-9]+).txt$") 
    if (j == 1){
      fun.genes.perm <- fun.files 
    } else {
      fun.genes.perm <- intersect(x = fun.genes.perm, y = fun.files)
    }
  }
  if (!identical(sort(fun.genes.orig), sort(fun.genes.perm))){
    stop("<< At least one gene is missing either (1) permuted summary statistics or (2) a gene score (following gene-scoring with the original summary statistics). Terminating.")
  } else {
    message(">> Everything looks OK. Continuing.")
  }
  return(fun.genes.orig)
}

extractRawCol <- function(data.raw, index){
  return(sapply(X = stri_extract_all_regex(str = data.raw, pattern = "\\S+"), FUN = function(x){x[index]}))
}

readRaw <- function(path.raw, keep.comment = "off"){
  fun.raw <- readLines(con = path.raw)
  if (keep.comment != "on"){
    fun.raw <- fun.raw[which(!(startsWith(x = fun.raw, prefix = "#")))]
  }
  return(fun.raw)
}

# vii. aggregation of permuted summary statistics

clusterGenes <- function(data.sums, data.annot.tbl, list.genes){
  message("Defining clusters of genes for which summary statistics should not be aggregated...")
  fun.annot.tbl <- copy(data.annot.tbl)
  fun.clusters  <- fun.annot.tbl[id %in% data.sums[, id], .(gene, id)]
  fun.clusters  <- graph_from_data_frame(d = as.data.frame(x = fun.clusters, stringsAsFactors = FALSE), directed = FALSE)
  fun.clusters  <- na.omit(clusters(graph = fun.clusters)[["membership"]][list.genes])
  fun.clusters  <- cbind.data.frame(gene = names(fun.clusters), cluster = unname(fun.clusters), stringsAsFactors = FALSE)
  return(fun.clusters)
}

orderClusterGenes <- function(data.cluster.genes, col.to.order, vec.of.order){
  message("Ordering genes according to the cluster to which they belong, from largest cluster to smallest cluster...")
  fun.cluster.genes <- data.cluster.genes
  fun.cluster.genes[, col.to.order] <- factor(x = fun.cluster.genes[, col.to.order], levels = vec.of.order, ordered = TRUE)
  fun.cluster.genes <- fun.cluster.genes[order(fun.cluster.genes[, col.to.order]), ]
  fun.cluster.genes[, col.to.order] <- as.character(fun.cluster.genes[, col.to.order])
  row.names(fun.cluster.genes) <- NULL
  return(fun.cluster.genes)
} 

combineSums <- function(path.annot, data.sums, data.annot.tbl, list.genes, path.main, name.sums, name.aug, type.sums.in = "permuted", type.sums.out = "aggregated", count.perms, count.cores){
  fun.annot <- readAnnot(path.annot = path.annot)
  fun.genes <- clusterGenes(data.sums = data.sums, data.annot.tbl = data.annot.tbl, list.genes = list.genes)
  fun.sizes <- sort(x = tapply(X = fun.genes[["gene"]], INDEX = fun.genes[["cluster"]], FUN = length), decreasing = TRUE)
  fun.order <- as.integer(names(fun.sizes))
  fun.genes <- orderClusterGenes(data.cluster.genes = fun.genes, col.to.order = "cluster", vec.of.order = fun.order)
  fun.wsize <- max(fun.sizes)
  fun.path  <- paste0(path.main, name.sums, "/", name.aug)
  fun.cores <- makeForkCluster(nnodes = as.integer(count.cores))
  registerDoParallel(fun.cores)
  on.exit(stopCluster(fun.cores))
  message("Largest cluster contains ", as.character(fun.wsize), " genes.")
  message("Using a sliding window to define ", as.character(fun.wsize), " batches of independent genes, and aggregating (for each permutation) their permuted summary statistics")
  for (j in 1:fun.wsize){
    fun.group <- fun.genes[seq(from = j, to = nrow(fun.genes), by = fun.wsize), "gene"]
    fun.index <- rep(NA, length(fun.group))
    for (i in 1:as.integer(count.perms)){
      fun.batch <-  foreach(h = 1:length(fun.group), .combine = rbind, .inorder = FALSE) %dopar%
        readSums(path.sums = paste0(fun.path, "/sumstat/", type.sums.in, "/", fun.group[h], "_", as.character(i), ".txt"), col.id = 2, col.p = 3, col.n = 4)
      fun.prefix <- paste0("C-", as.character(j), "_P-", as.character(i))
      writeSums(data.sums = fun.batch, path.main = path.main, name.sums = name.sums, name.aug = name.aug, type.sums.out = type.sums.out, prefix.sums.out = fun.prefix, log.verbose = "off")
      if (i == as.integer(count.perms)){
        if (TRUE %in% duplicated(fun.batch[, id])){
          stop("<< Summary statistics were aggregated for at least one pair of genes sharing a variant. Terminating.")
        } else {
          message(">> Completed batch: ", as.character(j), "/", as.character(fun.wsize))
        }
        for (h in 1:length(fun.group)){
          fun.index[h] <- which(stri_startswith(str = fun.annot, fixed = paste0(fun.group[h], "\t")))
        }
      }
    }
    fwrite(x = list(fun.group), file = paste0(fun.path, "/sumstat/", type.sums.out, "/", "C-", as.character(j), ".list"), quote = FALSE, sep = "\t", col.names = FALSE)
    writeLines(text = fun.annot[fun.index], con = paste0(fun.path, "/sumstat/", type.sums.out, "/", "C-", as.character(j), ".genes.annot"), sep = "\n")
  }
  return(fun.wsize)
} 

# viii. building new MAGMA .raw files for permuted data

readGenesOut <- function(path.genes.out){
  return(fread(path.genes.out, sep = " ", fill = TRUE, header = TRUE))
}

readGenesList <- function(path.genes.list){
  return(fread(path.genes.list, sep = "\t", header = FALSE, col.names = "GENE", fill = TRUE))
}

checkBatches <- function(path.batch.scores, count.batches, count.perms){
  fun.files <- list.files(path = path.batch.scores, pattern = "^C-[0-9]+_P-[0-9]+\\.genes\\.out$")
  fun.files <- gsub(pattern = "_P-[0-9]+\\.genes\\.out$", replacement = "", x = fun.files)
  fun.files <- gsub(pattern = "^C-", replacement = "", x = fun.files)
  fun.files <- sort(x = as.integer(fun.files))
  if (FALSE %in% (sort(x = rep(x = 1:count.batches, as.integer(count.perms))) == fun.files)){
    stop("<< Expect ", count.perms, " files for each batch of independent genes. At least one file is missing. Terminating.")
  } else {
    fun.batches <- 1:count.batches
    return(fun.batches)
  }
}

makeRepTblEmpty <- function(count.row, count.col){
  fun.tbl <- as.data.table(matrix(nrow = count.row, ncol = count.col))
  setnames(x = fun.tbl, old = paste0("V", as.character(1:count.col)), new = c("GENE", paste0("P", as.character(1:(count.col-1)))))
  fun.tbl[, GENE:=as.character(GENE)]
  fun.tbl <- list(fun.tbl, copy(fun.tbl))
  names(fun.tbl) <- c("ZSTAT", "N")
  for (index in 2:(as.integer(count.col))){
    set(x = fun.tbl[["ZSTAT"]], j = index, value = as.numeric(fun.tbl[["ZSTAT"]][[index]]))
    set(x = fun.tbl[["N"]], j = index, value = as.integer(fun.tbl[["N"]][[index]]))
  }
  return(fun.tbl)
}

makeRepTbl <- function(path.main, name.sums, name.aug, count.batches, count.perms, type.score.orig.in, prefix.score.orig.in, type.score.perm.in, type.sums.in, count.cores){
  message("Building reference tables of gene scores and sample sizes for rapid construction of .raw files.")
  fun.path <- paste0(path.main, name.sums, "/", name.aug)
  fun.orig <- readGenesOut(path.genes.out = paste0(fun.path, "/scores/", type.score.orig.in, "/", prefix.score.orig.in, ".genes.out"))
  fun.batches <- checkBatches(path.batch.scores = paste0(fun.path, "/scores/", type.score.perm.in), count.batches = count.batches, count.perms = count.perms) 
  fun.rep.tbl <- makeRepTblEmpty(count.row = nrow(fun.orig), count.col = (as.integer(count.perms) + 1))
  fun.start.i <- 1 
  fun.cores <- makeForkCluster(nnodes = as.integer(count.cores))
  registerDoParallel(fun.cores)
  on.exit(stopCluster(fun.cores))
  for (j in fun.batches){
    fun.genes <- readGenesList(path.genes.list = paste0(fun.path, "/sumstat/", type.sums.in, "/C-", as.character(j), ".list"))
    fun.shift <- nrow(fun.genes)
    fun.rep.tbl <- lapply(X = fun.rep.tbl, function(x){x[fun.start.i:(fun.start.i + fun.shift - 1), GENE:=as.character(fun.genes[, GENE])]})
    fun.input <- foreach(i = 1:as.integer(count.perms), .inorder = TRUE) %dopar% 
      readGenesOut(path.genes.out = paste0(fun.path, "/scores/", type.score.perm.in, "/C-", as.character(j), "_P-", as.character(i), ".genes.out"))
    fun.input <- lapply(X = fun.input, FUN = function(x){x[match(x = fun.genes[, GENE], table = x[, GENE]), ]})
    for (g in 1:as.integer(count.perms)){
      fun.rep.tbl[["ZSTAT"]][fun.start.i:(fun.start.i + fun.shift -1), (g + 1):=fun.input[[g]][, ZSTAT]]
      fun.rep.tbl[["N"]][fun.start.i:(fun.start.i + fun.shift -1), (g + 1):=fun.input[[g]][, N]] 
    }
    fun.start.i <- fun.start.i + fun.shift
  }
  fun.order <- fun.orig[, GENE]
  fun.rep.tbl[["ZSTAT"]] <- fun.rep.tbl[["ZSTAT"]][match(x = fun.order, table = fun.rep.tbl[["ZSTAT"]][, GENE]), ]
  fun.rep.tbl[["N"]] <- fun.rep.tbl[["N"]][match(x = fun.order, table = fun.rep.tbl[["N"]][, GENE]), ]
  if(!(identical(fun.rep.tbl[["ZSTAT"]][, GENE], fun.rep.tbl[["N"]][, GENE]))){
    stop("<< Order of genes in reference table for gene scores does not match that in reference table for sample sizes. Terminating.")
  } else {
    message("Reference tables generated. Using reference tables to generate one .raw file for each permutation...")
    return(fun.rep.tbl)
  }
}


makeRawMtrcs <- function(data.rep.tbl, path.main, name.sums, name.aug, type.score.in, prefix.score.in, type.score.out, count.perms, count.cores){
  fun.orig <- readRaw(path.raw = paste0(path.main, name.sums, "/", name.aug, "/scores/", type.score.in, "/", prefix.score.in, ".genes.raw"), keep.comment = "on")
  fun.comm <- fun.orig[stri_startswith_fixed(str  = fun.orig, pattern = "#")]
  fun.orig <- fun.orig[!stri_startswith_fixed(str = fun.orig, pattern = "#")]
  fun.orig <- stri_split_fixed(str = fun.orig, pattern = " ")
  fun.cores <- makeForkCluster(nnodes = as.integer(count.cores))
  registerDoParallel(fun.cores)
  on.exit(stopCluster(fun.cores))
  if (FALSE %in% (sapply(fun.orig, function(x){x[1]}) == data.rep.tbl[["ZSTAT"]][, GENE])){
    stop("Order of genes between original raw-file and replacement table does not match. Terminating.")
  } else {
    foreach(j = 1:as.integer(count.perms)) %dopar%
      lapply(X   = c(fun.comm, updRawMtrx(data.raw = fun.orig, data.rep.tbl = data.rep.tbl, data.perm = j)),
             FUN = function(x){fwrite(x     = list(paste(x, collapse = " ")), 
                                      file  = paste0(path.main, name.sums, "/", name.aug, "/scores/", type.score.out, "/scores.permutation-", as.character(j), ".unadjusted.genes.raw"), 
                                      quote = FALSE, sep = " ", col.names = FALSE, append = TRUE)})
  }
  message("Done. One .raw file for each permutation was written to the ", path.main, name.sums, "/", name.aug, "/scores/", type.score.out, "/ directory.")
  return(invisible())
}

updRawMtrx <- function(data.raw, data.rep.tbl, data.perm){
  fun.raw <- data.raw
  for (i in 1:length(fun.raw)){
    fun.raw[[i]][7] <- data.rep.tbl[["N"]][[data.perm + 1]][i]
    fun.raw[[i]][9] <- data.rep.tbl[["ZSTAT"]][[data.perm + 1]][i]
  }
  return(fun.raw)
}

# ix. making gene-set results summary (convenience)

summarizeSets <- function(path.main, name.sums, name.aug, prefix.bsl, prefix.aug, prefix.rnd, count.perms){
  fun.path  <- paste0(path.main, name.sums, "/", name.aug, "/")
  fun.bsl   <- fread(file = paste0(fun.path, "scores/baseline/",  prefix.bsl, ".gsa.out"), header = TRUE, fill = TRUE, skip = "VARIABLE", select = c("FULL_NAME", "BETA", "SE"), col.names = c("geneset", "beta", "se"))
  setkey(setnames(x = fun.bsl[, z:=getZ(vct.beta = fun.bsl[, beta], vct.se = fun.bsl[, se])][, beta:=NULL][, se:=NULL], c("geneset", "z_bsl")), geneset)
  fun.aug   <- fread(file = paste0(fun.path, "scores/augmented/", prefix.aug, ".gsa.out"), header = TRUE, fill = TRUE, skip = "VARIABLE", select = c("FULL_NAME", "BETA", "SE"), col.names = c("geneset", "beta", "se"))
  setkey(setnames(x = fun.aug[, z:=getZ(vct.beta = fun.aug[, beta], vct.se = fun.aug[, se])][, beta:=NULL][, se:=NULL], c("geneset", "z_aug")), geneset)
  fun.out   <- merge(x = fun.bsl, y = fun.aug, by = "geneset")
  for (j in 1:as.integer(count.perms)){
    if (j == 1){
      fun.rnd <- fread(file = paste0(fun.path, "scores/random/",  prefix.rnd, as.character(j), ".adjusted.gsa.out"), header = TRUE, fill = TRUE, skip = "VARIABLE", select = c("FULL_NAME", "BETA", "SE"), col.names = c("geneset", "beta", "se"))
      setkey(setnames(x = fun.rnd[, z:=getZ(vct.beta = fun.rnd[, beta], vct.se = fun.rnd[, se])][, beta:=NULL][, se:=NULL], c("geneset", paste0("z_per-", as.character(j)))), geneset)
    } else {
      fun.tmp <- fread(file = paste0(fun.path, "scores/random/",  prefix.rnd, as.character(j), ".adjusted.gsa.out"), header = TRUE, fill = TRUE, skip = "VARIABLE", select = c("FULL_NAME", "BETA", "SE"), col.names = c("geneset", "beta", "se"))
      setkey(setnames(x = fun.tmp[, z:=getZ(vct.beta = fun.tmp[, beta], vct.se = fun.tmp[, se])][, beta:=NULL][, se:=NULL], c("geneset", paste0("z_per-", as.character(j)))), geneset)
      fun.rnd <- merge(x = fun.rnd, y = fun.tmp, by = "geneset") 
    }
  }
  fun.rnd <- cbind.data.frame(geneset = fun.rnd[, geneset], z_per = apply(X = fun.rnd[, 2:(j+1)], MARGIN = 1, FUN = mean), sd_per = apply(X = fun.rnd[, 2:(j+1)], MARGIN = 1, FUN = sd))
  fun.out <- merge(x = fun.out, y = fun.rnd, by = "geneset")
  fun.msg <- copy(fun.out)
  fun.msg <- fun.msg[z_aug > (qnorm(p = 0.05, lower.tail = FALSE))]
  message("Number of gene sets that are signficantly enriched for phenotype association with the augmented model: ", nrow(fun.msg))
  if (nrow(fun.msg) > 0){
    fun.msg <- fun.msg[z_aug > z_bsl]
    message("Amongst these, ", nrow(fun.msg), " gene sets are more enriched for phenotype association with the augmented model compared to with the baseline model")
    if (nrow(fun.msg) > 0){
      fun.msg[, z_val:=((z_aug-z_per)/sd_per)][, category:=sapply(X = fun.msg[, z_val], FUN = assignConfidence)]
      message("Amongst these gene sets, ", as.character(nrow(fun.msg[category == "Invalidated"])), " are invalidated, ", as.character(nrow(fun.msg[category == "Mildy Validated"])), " are mildy validated, and ", as.character(nrow(fun.msg[category == "Strongly Validated"])), " are strongly validated.")
    }
  }
  setnames(x = fun.out[, remarkable:=ifelse(test = fun.out[, z_aug] > (qnorm(p = 0.05, lower.tail = FALSE)), yes = "yes", no = "no")
                       ][, gain:=ifelse(test = fun.out[, z_aug] > fun.out[, z_bsl], yes = "gain", no = "no gain")
                         ][, z_val:=((z_aug-z_per)/sd_per)
                           ][, category:=sapply(X = fun.out[, z_val], FUN = assignConfidence)
                             ][, z_val:=NULL], c("geneset", "z_baseline", "z_augmented", "z_random-mean", "z_random-sd", "remarkable", "gain-or-not", "validation-category"))
  message("Writing summary of gene set results to: ", fun.path, "genesets.summary")
  fwrite(x = fun.out, file = paste0(fun.path, "genesets.summary"), quote = FALSE, sep = "\t", col.names = TRUE)
  return(invisible())
}


getZ <- function(vct.beta, vct.se){
  fun.p <- pnorm(q = (vct.beta/vct.se), lower.tail = FALSE)
  fun.q <- p.adjust(p = fun.p, method = "fdr")
  fun.q <- ifelse(fun.q < (1*(10^(-323))), yes = (1*(10^(-323))), no = fun.q)
  fun.q <- ifelse(fun.q > 0.9999999999999999, yes = 0.9999999999999999, no = fun.q)
  fun.z <- qnorm(p = fun.q, lower.tail = FALSE)
  return(fun.z)
} 

assignConfidence <- function(x){
  if (x >= 2){
    "Strongly Validated"
  } else if (x >= 1 & x < 2){
    "Mildy Validated"
  } else {
    "Invalidated"
  }
}


### PROCESS ##############################################################################################################################################################

# 1. parse and check command-line arguments
message("1. ", Sys.time(), ". Start. EPVP is a controlled way for augmenting a baseline annotation-model in MAGMA with larger flanks or regulatory interactions.")
command <- checkCmd(cmd.ready = rmEmptyFl(cmd.made = makeCmd(cmd.raw = paste(command, collapse = " ")))) 

message("")
# 2. read-in summary statistics and perform some simple filtering steps and quality controls
message("2. ", Sys.time(), ". Preparing summary statistics for analyses.")
sumstat <- prepSums(data.sums   = readSums(path.sums = command["sumstat"], col.id = command["sumstat-id"], col.p = command["sumstat-pval"], col.n = command["sumstat-nsample"]),
                    data.bim    = readBim(path.bins = command["binaries"]), 
                    path.main   = command["output"], 
                    name.sums   = command["study"], 
                    name.aug    = command["augmentation"],
                    count.perms = command["permutations"])

message("")
# 3. preparation of variant-to-gene annotation table (performed at this stage to check for potential issues in the files)
message("3. ", Sys.time(), ". Creating table of variant-to-gene annotation according to baseline and augmented model.")
annotation <- prepAnnotTbl(data.annot.bsl = reformatAnnot(data.annot = readAnnot(path.annot = command["baseline-model"])), 
                           data.annot.aug = reformatAnnot(data.annot = readAnnot(path.annot = command["augmented-model"])), 
                           data.sums      = sumstat,
                           path.main      = command["output"], 
                           name.sums      = command["study"], 
                           name.aug       = command["augmentation"])

message("")
# 4. execute magma gene-scoring and gene-set analysis for baseline model (non-permuted summary statistics)
message("4. ", Sys.time(), ". Executing MAGMA gene-scoring and gene-set analysis for baseline model using the prepared (non-permuted) summary statistics.")
scoreGenes(path.magma         = command["magma"],
           score.mdl          = command["gene-scoring-model"],
           path.bins          = command["binaries"],
           path.annot         = command["baseline-model"],
           path.main          = command["output"],
           name.sums          = command["study"],
           name.aug           = command["augmentation"],
           type.sums.in       = "original",
           prefix.sums.in     = "sumstat",
           type.score.out     = "baseline",
           prefix.score.out   = "scores.original.unadjusted")

scoreSets(path.magma       = command["magma"],
          path.sets        = command["gene-set-file"],
          format.sets      = command["gene-set-format"],
          path.main        = command["output"],
          name.sums        = command["study"],
          name.aug         = command["augmentation"],
          gene.exclude     = command["ignore-genes"],
          type.score       = "baseline",
          prefix.score.in  = "scores.original.unadjusted",
          prefix.score.out = "scores.original.adjusted")

message("")
# 5. execute magma gene-scoring and gene-set analysis for augmented model (non-permuted summary statistics)
message("5. ", Sys.time(), ". Executing MAGMA gene-scoring and gene-set analysis for augmented model using the prepared (non-permuted) summary statistics.")
scoreGenes(path.magma         = command["magma"],
           score.mdl          = command["gene-scoring-model"],
           path.bins          = command["binaries"],
           path.annot         = command["augmented-model"],
           path.main          = command["output"],
           name.sums          = command["study"],
           name.aug           = command["augmentation"],
           type.sums.in       = "original",
           prefix.sums.in     = "sumstat",
           type.score.out     = "augmented",
           prefix.score.out   = "scores.original.unadjusted")

scoreSets(path.magma       = command["magma"],
          path.sets        = command["gene-set-file"],
          format.sets      = command["gene-set-format"],
          path.main        = command["output"],
          name.sums        = command["study"],
          name.aug         = command["augmentation"],
          gene.exclude     = command["ignore-genes"],
          type.score       = "augmented",
          prefix.score.in  = "scores.original.unadjusted",
          prefix.score.out = "scores.original.adjusted")

message("")
# 6. executing EPVP permutation of summary statistics
message("6. ", Sys.time(), ". Permuting summary statistics (EPVP).")
permuteSums(data.sums      = sumstat,
            data.annot.tbl = annotation,
            count.perms    = command["permutations"],
            count.cores    = command["cores"],
            path.main      = command["output"],
            name.sums      = command["study"],
            name.aug       = command["augmentation"])

genes <- checkPermSums(path.main       = command["output"], 
                       name.sums       = command["study"], 
                       name.aug        = command["augmentation"], 
                       count.perms     = command["permutations"],
                       prefix.score.in = "scores.original.unadjusted")

message("")
# aggregating permuted summary statistics for efficiency
message("7. ", Sys.time(), ". Aggregating permuted summary statistics for independent collections of genes.")
batches <- combineSums(path.annot     = command["augmented-model"],
                       data.sums      = sumstat,
                       data.annot.tbl = annotation,
                       list.genes     = genes,
                       path.main      = command["output"],
                       name.sums      = command["study"],
                       name.aug       = command["augmentation"],
                       count.perms    = command["permutations"],
                       count.cores    = command["cores"],
                       type.sums.in   = "permuted",
                       type.sums.out  = "aggregated")

message("")
# 8. executing magma gene-scoring for aggregated files of summary statistics (parallel computing; specify cores with --nodes flag)
message("8. ", Sys.time(), ". Executing MAGMA gene scoring for the augmented model using permuted summary statistics (aggregated files).")
scoreGenesPar(path.magma      = command["magma"],
              path.bins       = command["binaries"],
              path.main       = command["output"],
              name.sums       = command["study"],
              name.aug        = command["augmentation"],
              count.cores     = command["cores"],
              score.mdl       = command["gene-scoring-model"],
              type.sums.in    = "aggregated",
              type.score.out  = "random/batches")

message("")
# 9. preparing a magma .raw file for each permutation
message("9. ", Sys.time(), ". Preparing a MAGMA .raw file for each permutation. These files are required for computing adjusted gene scores and gene set scores.")
reptable <- makeRepTbl(count.batches = batches,
                       path.main = command["output"],
                       name.sums = command["study"],
                       name.aug = command["augmentation"],
                       count.perms = command["permutations"],
                       count.cores = command["cores"],
                       type.score.orig.in = "augmented",
                       prefix.score.orig.in = "scores.original.unadjusted",
                       type.score.perm.in = "random/batches/",
                       type.sums.in = "aggregated")

makeRawMtrcs(data.rep.tbl    = reptable,
             path.main       = command["output"],
             name.sums       = command["study"],
             name.aug        = command["augmentation"],
             type.score.in   = "augmented",
             prefix.score.in = "scores.original.unadjusted",
             type.score.out  = "random",
             count.perms     = command["permutations"],
             count.cores     = command["cores"])

message("")
# 10. executing gene-set analysis for each permutation
message("10. ", Sys.time(), ". Executing gene-set analysis for each permuted .raw file.")
scoreSetsPar(path.magma   = command["magma"],
             path.sets    = command["gene-set-file"],
             path.main    = command["output"],
             name.sums    = command["study"],
             name.aug     = command["augmentation"],
             format.sets  = command["gene-set-format"],
             gene.exclude = command["ignore-genes"],
             count.cores  = command["cores"],
             type.score   = "random")

message("")
# 11. completion announcement
message("11. ", Sys.time(), ". Completed all analyses. Generating a simple results summary for reference...")
summarizeSets(path.main   = command["output"], 
              name.sums   = command["study"], 
              name.aug    = command["augmentation"], 
              count.perms = command["permutations"],
              prefix.bsl  = "scores.original.adjusted", 
              prefix.aug  = "scores.original.adjusted", 
              prefix.rnd  = "scores.permutation-")

















