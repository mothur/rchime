# Detect and remove chimeras from your [strollur](https://mothur.org/strollur/) object or data.frame using a denovo approach or alternatively a reference model.

The `rchime()` function allows you to detect and remove chimeras from
your data using a denovo approach or alternatively a reference model.

Our preferred way of doing this is to use the abundant sequences as our
reference (denovo).

This function uses code from the
[vsearch](https://github.com/torognes/vsearch) tools.

## Usage

``` r
rchime(
  data,
  reference = NULL,
  dereplicate = TRUE,
  verbose = TRUE,
  remove_chimeras = TRUE,
  rchime_options = NULL,
  table_names = list(sequence_name = "sequence_name", sequence = "sequence", abundance =
    "abundance", sample = "sample")
)
```

## Arguments

- data:

  a [strollur](https://mothur.org/strollur/) dataset object or a
  data.frame containing your sequence data.

- reference:

  a [strollur](https://mothur.org/strollur/) dataset object or a
  data.frame containing reference sequence data.

- dereplicate:

  logical. The dereplicate option allows you to remove chimeras by
  sample. When `dereplicate=FALSE`, if a sequence is flagged as chimeric
  in one sample, it is removed from all samples. Our experience suggests
  that this is a bit aggressive since we’ve seen rare sequences get
  flagged as chimeric when they’re the most abundant sequence in another
  sample. For a more conservative approach, we recommend using the
  default `dereplicate=TRUE` which will only remove sequences from the
  samples in which they are flagged as chimeric.

- verbose:

  logical, allow console outputs. Default = `TRUE`.

- remove_chimeras:

  Boolean, remove chimeras from dataset. Default = `TRUE`. Only used
  when `data` is a strollur object.

- rchime_options:

  List, You can fine tune the vsearch specific options using the
  \[[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)\]
  function. Default = `NULL`.

- table_names, :

  named list used to indicate the names of the columns in the
  data.frame. Only used when `data` or `reference` are a data.frames.

## Value

list() containing a chimera report, and vector of the chimeric
sequence's names.

The [strollur](https://mothur.org/strollur/) dataset object will also be
updated. The sequences flagged as chimeric are removed and the chimera
report is added.

## References

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a
versatile open source tool for metagenomics. PeerJ 4:e2584. doi:
10.7717/peerj.2584

Edgar,R.C., Haas,B.J., Clemente,J.C., Quince,C. and Knight,R. (2011),
UCHIME improves sensitivity and speed of chimera detection.
Bioinformatics 27:2194.

## See also

[`rchime_options()`](http://mothur.org/rchime/reference/rchime_options.md)
to set vsearch specific parameters.

## Author

Sarah Westcott, <swestcot@umich.edu>

## Examples

``` r

# Let's use a strollur object with 100 sequences

data <- strollur::load_dataset(
  rchime_example("strollur_multi_sample_tiny.rds")
)
#> Error in system.file("extdata", file, package = "rchime", mustWork = TRUE): no file found

chimera_report <- rchime(data)
#> Error in rchime.default(data): data must be a strollur object or a data.frame
data
#> function (..., list = character(), package = NULL, lib.loc = NULL, 
#>     verbose = getOption("verbose"), envir = .GlobalEnv, overwrite = TRUE) 
#> {
#>     fileExt <- function(x) {
#>         db <- grepl("\\.[^.]+\\.(gz|bz2|xz)$", x)
#>         ans <- sub(".*\\.", "", x)
#>         ans[db] <- sub(".*\\.([^.]+\\.)(gz|bz2|xz)$", "\\1\\2", 
#>             x[db])
#>         ans
#>     }
#>     my_read_table <- function(...) {
#>         lcc <- Sys.getlocale("LC_COLLATE")
#>         on.exit(Sys.setlocale("LC_COLLATE", lcc))
#>         Sys.setlocale("LC_COLLATE", "C")
#>         read.table(...)
#>     }
#>     stopifnot(is.character(list))
#>     names <- c(as.character(substitute(list(...))[-1L]), list)
#>     if (!is.null(package)) {
#>         if (!is.character(package)) 
#>             stop("'package' must be a character vector or NULL")
#>     }
#>     paths <- find.package(package, lib.loc, verbose = verbose)
#>     if (is.null(lib.loc)) 
#>         paths <- c(path.package(package, TRUE), if (!length(package)) getwd(), 
#>             paths)
#>     paths <- unique(normalizePath(paths[file.exists(paths)]))
#>     paths <- paths[dir.exists(file.path(paths, "data"))]
#>     dataExts <- tools:::.make_file_exts("data")
#>     if (length(names) == 0L) {
#>         db <- matrix(character(), nrow = 0L, ncol = 4L)
#>         for (path in paths) {
#>             entries <- NULL
#>             packageName <- if (file_test("-f", file.path(path, 
#>                 "DESCRIPTION"))) 
#>                 basename(path)
#>             else "."
#>             if (file_test("-f", INDEX <- file.path(path, "Meta", 
#>                 "data.rds"))) {
#>                 entries <- readRDS(INDEX)
#>             }
#>             else {
#>                 dataDir <- file.path(path, "data")
#>                 entries <- tools::list_files_with_type(dataDir, 
#>                   "data")
#>                 if (length(entries)) {
#>                   entries <- unique(tools::file_path_sans_ext(basename(entries)))
#>                   entries <- cbind(entries, "")
#>                 }
#>             }
#>             if (NROW(entries)) {
#>                 if (is.matrix(entries) && ncol(entries) == 2L) 
#>                   db <- rbind(db, cbind(packageName, dirname(path), 
#>                     entries))
#>                 else warning(gettextf("data index for package %s is invalid and will be ignored", 
#>                   sQuote(packageName)), domain = NA, call. = FALSE)
#>             }
#>         }
#>         colnames(db) <- c("Package", "LibPath", "Item", "Title")
#>         footer <- if (missing(package)) 
#>             paste0("Use ", sQuote(paste("data(package =", ".packages(all.available = TRUE))")), 
#>                 "\n", "to list the data sets in all *available* packages.")
#>         else NULL
#>         y <- list(title = "Data sets", header = NULL, results = db, 
#>             footer = footer)
#>         class(y) <- "packageIQR"
#>         return(y)
#>     }
#>     paths <- file.path(paths, "data")
#>     for (name in names) {
#>         found <- FALSE
#>         for (p in paths) {
#>             tmp_env <- if (overwrite) 
#>                 envir
#>             else new.env()
#>             if (file_test("-f", file.path(p, "Rdata.rds"))) {
#>                 rds <- readRDS(file.path(p, "Rdata.rds"))
#>                 if (name %in% names(rds)) {
#>                   found <- TRUE
#>                   if (verbose) 
#>                     message(sprintf("name=%s:\t found in Rdata.rds", 
#>                       name), domain = NA)
#>                   objs <- rds[[name]]
#>                   lazyLoad(file.path(p, "Rdata"), envir = tmp_env, 
#>                     filter = function(x) x %in% objs)
#>                   break
#>                 }
#>                 else if (verbose) 
#>                   message(sprintf("name=%s:\t NOT found in names() of Rdata.rds, i.e.,\n\t%s\n", 
#>                     name, paste(names(rds), collapse = ",")), 
#>                     domain = NA)
#>             }
#>             files <- list.files(p, full.names = TRUE)
#>             files <- files[grep(name, files, fixed = TRUE)]
#>             if (length(files) > 1L) {
#>                 o <- match(fileExt(files), dataExts, nomatch = 100L)
#>                 paths0 <- dirname(files)
#>                 paths0 <- factor(paths0, levels = unique(paths0))
#>                 files <- files[order(paths0, o)]
#>             }
#>             if (length(files)) {
#>                 for (file in files) {
#>                   if (verbose) 
#>                     message("name=", name, ":\t file= ...", .Platform$file.sep, 
#>                       basename(file), "::\t", appendLF = FALSE, 
#>                       domain = NA)
#>                   ext <- fileExt(file)
#>                   if (basename(file) != paste0(name, ".", ext)) 
#>                     found <- FALSE
#>                   else {
#>                     found <- TRUE
#>                     switch(ext, R = , r = {
#>                       library("utils")
#>                       sys.source(file, chdir = TRUE, envir = tmp_env)
#>                     }, RData = , rdata = , rda = load(file, envir = tmp_env), 
#>                       TXT = , txt = , tab = , tab.gz = , tab.bz2 = , 
#>                       tab.xz = , txt.gz = , txt.bz2 = , txt.xz = assign(name, 
#>                         my_read_table(file, header = TRUE, as.is = FALSE), 
#>                         envir = tmp_env), CSV = , csv = , csv.gz = , 
#>                       csv.bz2 = , csv.xz = assign(name, my_read_table(file, 
#>                         header = TRUE, sep = ";", as.is = FALSE), 
#>                         envir = tmp_env), found <- FALSE)
#>                   }
#>                   if (found) 
#>                     break
#>                 }
#>                 if (verbose) 
#>                   message(if (!found) 
#>                     "*NOT* ", "found", domain = NA)
#>             }
#>             if (found) 
#>                 break
#>         }
#>         if (!found) {
#>             warning(gettextf("data set %s not found", sQuote(name)), 
#>                 domain = NA)
#>         }
#>         else if (!overwrite) {
#>             for (o in ls(envir = tmp_env, all.names = TRUE)) {
#>                 if (exists(o, envir = envir, inherits = FALSE)) 
#>                   warning(gettextf("an object named %s already exists and will not be overwritten", 
#>                     sQuote(o)))
#>                 else assign(o, get(o, envir = tmp_env, inherits = FALSE), 
#>                   envir = envir)
#>             }
#>             rm(tmp_env)
#>         }
#>     }
#>     invisible(names)
#> }
#> <bytecode: 0x55b1014a3808>
#> <environment: namespace:utils>
```
