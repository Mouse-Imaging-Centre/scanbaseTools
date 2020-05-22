#' @param db The google sheet title.
#' @export
read_scanbase <- function(db = "test_scanbase_40um"){
  tfile <- tempfile(fileext = ".xlsx")
  googledrive::drive_download(db, path = tfile)
  list(scans = readxl::read_excel(tfile, "Scans")
     , studies = readxl::read_excel(tfile, "Studies")
     , genotypes = readxl::read_excel(tfile, "Genotypes")
     , treatments = readxl::read_excel(tfile, "Treatments")
     , global = readxl::read_excel(tfile, "Global")
       )
}

#' 

#' @export
load_pydpiper_results <-
  function(ppd
         , common = "/hpf/largeprojects/MICe/scanbase/pipeline-40um/scanbase_second_level_nlin/scanbase_second_level-nlin-3.mnc"
         , mask = sub("\\.mnc$", "_mask.mnc", common)
         , fwhm = 0.2
         , clobber = FALSE
         , dry = FALSE){

    ppd <- normalizePath(ppd, mustWork = TRUE)
    
    if(!file.exists(common))
      stop("Unable to find your common space file: ", common)
    
    if(!file.exists(mask))
      stop("Unable to find your common space mask: ", mask)
    
    transforms <- read.csv(file.path(ppd, "transforms.csv")
                         , stringsAsFactors = FALSE)

    if(is.null(transforms[["overall_xfm_to_common"]])){
      if(is.null(transforms[["overall_xfm_to_common_inv"]])){
        stop("Neither overall_xfm_to_common nor overall_xfm_to_common_inv found in your transforms.csv, can't compute determinants.\\n Aborting.")
      } else {
        forwards <-
          transforms$overall_xfm_to_common_inv %>%
          sub("_inverted\\.xfm", ".xfm", .)

        forwards_full <-
          forwards %>%
          ifelse(grepl("^/", .), ., file.path(ppd, .))
        
        if(!all(file.exists(file.path(forwards_full))))
          stop("Only overall_xfm_to_common_inv found in transforms, removing `_inverted` suffix does not produce valid files.\\n Aborting.")
        
        transforms$overall_xfm_to_common <- forwards
      }
    }

    determinants <-
      read.csv(file.path(ppd, "determinants.csv")
             , stringsAsFactors = FALSE) %>%
      filter(fwhm == 0.2)

    column_mapping <-
      c(
        Distortion_Corrected_Scan = "native_file"
      , Scan_To_Study_Absolute_Jacobians = "log_full_det"
      , Scan_To_Study_Relative_Jacobians = "log_nlin_det"
      , Scan_To_Global_Absolute_Jacobians = "log_full_overall_det_common"
      , Scan_To_Global_Relative_Jacobians = "log_nlin_overall_det_common"
      , Scan_To_Study_Global_Space_Resampled_Absolute_Jacobians = "log_full_det_common"
      , Scan_To_Study_Global_Space_Resampled_Relative_Jacobians = "log_nlin_det_common"
      , Rigid_Transform = "rigid_xfm"
      , Rigid_Filepath = "lsq6_file"
      , Scan_To_Study_Transform = "lsq12_nlin_xfm"
      , Labels = "label_file"
      )

    known_columns <- keep(column_mapping, ~ . %in% names(transforms) || . %in% names(determinants))

    full_data <-
      inner_join(transforms, determinants
               , by = c("lsq12_nlin_xfm" = "inv_xfm")) %>%
      rename(
        !!! map(known_columns, as.symbol)
      ) %>%
      mutate_at(
        .vars = vars(!!!names(known_columns), overall_xfm_to_common)
      , .funs = funs(
          ifelse(grepl("^/", .), ., file.path(ppd, .)))
      ) %>%
      mutate(Processed_dir = dirname(dirname(Scan_To_Study_Relative_Jacobians))
           , Labels =
                 `if`(is.null(.$Labels)
                    , map_chr(file.path(Processed_dir, "voted.mnc")
                            , function(f){
                              if(file.exists(f)){
                                return(f)
                              } else {
                                stop("these subjects have not been processed with MAGeT")
                              }
                            })
                    , Labels
                      ))

    scans <-
      full_data %>%
      mutate(
        xfm_nlin3_to_global =
          map_chr(overall_xfm_to_common
                , ~ global_from_concat(., "nlin3-to-global.xfm"
                                     , clobber = clobber
                                    ,  dry = dry))) 
    
    
    if(is.null(scans$Scan_To_Study_Global_Space_Resampled_Absolute_Jacobians)){
      scans <- scans %>%
        mutate(Scan_To_Study_Global_Space_Resampled_Absolute_Jacobians =
                 future_map2_chr(overall_xfm_to_common
                                 , basename(Processed_dir)
                                 , ~ compute_inverse_determinants(.x
                                                                  , "abs"
                                                                  , like = common
                                                                  , output = .y
                                                                  , mask = mask
                                                                  , clobber = clobber
                                                                  , dry = dry))
               )
    }
        
    if(is.null(scans$Scan_To_Study_Global_Space_Resampled_Relative_Jacobians)){
      scans <- scans %>%
        mutate(Scan_To_Study_Global_Space_Resampled_Relative_Jacobians =
                 future_map2_chr(overall_xfm_to_common
                                 , basename(Processed_dir)
                                 , ~ compute_inverse_determinants(.x
                                                                  , "rel"
                                                                  , like = common
                                                                  , output = .y
                                                                  , mask = mask
                                                                  , clobber = clobber
                                                                  , dry = dry))
        )
    } 
    
    scans <- scans %>%
      select(!!! map(names(column_mapping), as.symbol))

    
    command_and_version_to_date <-
      function(cav){
        cav %>%
          basename() %>%
          sub(".*-([^-]*-[^-]*-[^-]*-at-[^-]*-[^-]*-[^-]*)\\.sh", "\\1", ., perl = TRUE) %>%
          as.POSIXct(format = "%d-%m-%Y-at-%H-%M-%OS")
      }

    get_first_command <-
      function(){
        quiet_readLines <- function(x){
          withCallingHandlers(
            readLines(x)
          , warning = function(w){
            if(!grepl("incomplete final line", w)){
              warning(w)
            } else {
              invokeRestart("muffleWarning")
            }
          })
        }
          
        Sys.glob(file.path(ppd, "*command-and-version*")) %>%
          setNames(command_and_version_to_date(.), .) %>%
          sort(decreasing = TRUE) %>%
          head(1) %>%
          names %>%
          readLines %>%
          (function(lines){
            vers <- lines[6] %>% sub(".*is: *", "", .)
            cmd <- lines[8] %>% strsplit(" ") %>% .[[1]] %>% .[1]
            paste0(cmd, " v", vers)
          })
      }
    
    get_study_to_global_xfm <- function(){
      xfm <- full_data$overall_xfm_to_common[1]
      xfm_dir <- dirname(xfm)
      
      files <-
        xfm %>%
        scan(what = character(), quiet = TRUE) %>%
        grep("*.mnc", ., value = TRUE) %>%
        file.path(xfm_dir, .) %>%
        c(xfm, .)

      print(files)
      new_files <- file.path(ppd, basename(files))
      
      walk2(files, new_files, ~ {
        if(!file.exists(.y) || clobber)
          file.copy(.x, .y)
      })

      study_to_global <- "study_to_global.xfm"
      global_from_concat(new_files[1], study_to_global, clobber = clobber, dry = dry)

      file.path(ppd, study_to_global)
    }

    study <-
      data_frame(Study_Average =
                   Sys.glob(file.path(ppd, "*nlin", "/*nlin-3.mnc"))
               , Study_To_Global_Transform =
                   get_study_to_global_xfm()
               , Study_Mask =
                   Sys.glob(file.path(ppd, "*nlin", "*nlin-3_mask.mnc"))
               , Pydpiper_Path =
                   get_first_command()
                 )

    if(nrow(study) != 1)
      stop("Error identifying study mask or average \n"
         , "Should match <pydpiper_dir>/*nlin/*nlin-3(_mask)?.mnc")        
    
    list(scans = scans, study = study)
  }



#' @export
upload_study <- function(ppd, scan, study, treatments, genotypes
                       , db = "test_scanbase_40um"
                       , common = "/hpf/largeprojects/MICe/scanbase/pipeline-40um/scanbase_second_level_nlin/scanbase_second_level-nlin-3.mnc"
                       , mask = sub("\\.mnc$", "_mask.mnc", common)
                       , clobber = FALSE ){

  sb <- read_scanbase(db)
  pp <- load_pydpiper_results(ppd, common, mask = mask, clobber = clobber)
  
  scanf <- read.csv(scan, stringsAsFactors = FALSE) %>% mutate(Mouse_ID = as.character(Mouse_ID))
  studf <- read.csv(study, stringsAsFactors = FALSE) %>% bind_cols(pp$study)
  treatf <- read.csv(treatments, stringsAsFactors = FALSE)
  genef <- read.csv(genotypes, stringsAsFactors = FALSE)

  set_diffs <- function(x,y){
    list(missing = setdiff(x, y)
       , extra = setdiff(y, x))
  }

  column_diff_test <- function(f1, f2, name){
  column_diffs <- set_diffs(names(f2), names(f1))
  if(length(c(column_diffs$missing, column_diffs$extra)) != 0)
    stop("Columns of ", name, " sheet update differ from those present in scanbase\n"
       , "MISSING: ", paste0(column_diffs$missing, collapse = ", "), "\n"
       , "EXTRA: ", paste0(column_diffs$extra, collapse = ", "))
  }

  # Check the study doesn't already exist in scanbase
  if(studf$Study_Name %in% sb$studies$Study_Name)
    stop("This study name already appears in scanbase, has it been uploaded already?")

  # Merge the pydpiper results with the scan frame
  scanf_merged <- inner_join(scanf, pp$scans, by = "Distortion_Corrected_Scan")

  if(nrow(scanf_merged) != nrow(scanf)){
    cat("Failed to merge pydpiper results with your scan sheet, trying resolving symlinks")

    scanf_merged <-
      pp$scans %>%
      mutate(Distortion_Corrected_Scan =
               map_chr(Distortion_Corrected_Scan, ~ system(paste("readlink -f", .), intern = TRUE))) %>%
      inner_join(scanf, . , by = "Distortion_Corrected_Scan")

    if(nrow(scanf_merged) != nrow(scanf))
      stop("Scans present in the pydpiper directory were not merged successfully, "
         , "this potentially indicates a problem with you scan metadata")    
  }

  if(nrow(scanf_merged) != nrow(pp$scans)){
    choice <- readline("There were scans found in your pydpiper directory that aren't in your scan sheet. Proceed? [Y/n]")
    while(! choice %in% c("", "Y", "y")){      
      if(choice == "n") stop("Aborting")
      choice <- readline("Please answer y or n")
    }
  }
   

  if(is.null(scanf_merged$POND_Mouse_ID))
    scanf_merged$POND_Mouse_ID <-
      seq_len(nrow(scanf_merged))

  column_diff_test(sb$scans, scanf_merged, "scans")
  column_diff_test(sb$studies, studf, "studies")
  column_diff_test(sb$treatments, treatf, "treatments")
  column_diff_test(sb$genotypes, genef, "genotypes")

  # Will error if there are type mismatches
  new_scans <- bind_rows(sb$scans, scanf_merged) %>% tail(nrow(scanf_merged))
  new_studies <- bind_rows(sb$studies, studf) %>% tail(nrow(studf))
  new_treatments <- bind_rows(sb$treatments, treatf) %>% tail(nrow(treatf))
  new_genotypes <- bind_rows(sb$genotypes, genef) %>% tail(nrow(genef))

  sheet <- googlesheets4::sheets_find(db) %>% filter(name == db)
  if(nrow(sheet) != 1) stop("Sheet ", db, " not found")
  id <- sheet$id[1]
  token <- googlesheets4::sheets_token()

  upload_sheet <- function(new, old, name){
    col <- cellranger::num_to_letter(ncol(old))
    old_max <- nrow(old) + 1
    new_max <- old_max + nrow(new)

    range_ref <- glue("{name}!A{old_max + 1}:{col}{new_max}")

    body <-
      list(range = range_ref
         , majorDimension = "ROWS"
         , values = as.matrix(new))

    httr::PUT(
            glue("https://sheets.googleapis.com/v4/spreadsheets/"
               , "{id}/"
               , "values/{range_ref}?valueInputOption=USER_ENTERED")
          , token
          , body = jsonlite::toJSON(body, auto_unbox = TRUE))
            
  }

  list(upload_sheet(new_scans, sb$scans, "Scans")
     , upload_sheet(new_studies, sb$studies, "Studies")
     , upload_sheet(new_treatments, sb$treatments, "Treatments")
     , upload_sheet(new_genotypes, sb$genotypes, "Genotypes"))
}

#' @export
unconcat_xfm <- function(xfm){
  lines <- readLines(xfm)
  xfm_starts <- grep("Transform_Type", lines)

  header <- lines[1:(xfm_starts[1] - 1)]
  xfms <-
    xfm_starts %>%
    { map2(., lead(., default = length(lines) + 1) - 1
         , function(b,e) lines[b:e]) }
                 
  list(header = header, xfms = xfms)
}

#' @export
reconstruct_xfm <- function(unc_xfms, which, file = NULL
                           , clobber = TRUE, dry = FALSE){
  header <- c(unc_xfms$header %>% { .[!grepl("^ *$", .)] }
            , glue("%{Sys.Date()}>>> Unconcat in R and reassembled from pieces:",
                   " {paste0(which, collapse = ', ')}"))

  xfm <- c(header, "", unlist(unc_xfms$xfms[which]))
  if(!is.null(file)){
    if(!file.exists(file) || clobber){
      if(dry){
        cat("Generating reconstructed xfm ", file, "\n")
      } else {
        cat(xfm, file = file, sep = "\n")
      }
      return(invisible(xfm))
    }
  }

  xfm
}

#' ((scan->nlin)->global) -> (nlin->global)
#'
#' Deconstruct a scan to global transform into an nlin to global transform
#' The output file is created in the same directory as the input file.
#' @param scan_to_global path to an xfm file composed of a nonlinear scan->nlin and a nonlinear
#' nlin->global. The xfm file should have four component transforms linear->grid->linear->grid.
#' The final two are extracted and written out.
#' @param nlin_to_global A filename for the resultant transform, to be written to the directory
#' of `scan_to_global`
#' @param dry whether to do a dry run, meaning not run the command.
#' @return the path to the output file invisibly.
#' @export
global_from_concat <- function(scan_to_global, nlin_to_global
                             , clobber = TRUE, dry = FALSE){
  dir <- dirname(scan_to_global)
  if(grepl("/", nlin_to_global)) stop("nlin_to_global cannot be a path")
  nlin_to_global <- file.path(dir, nlin_to_global)
  
  s2g <- unconcat_xfm(scan_to_global)
  reconstruct_xfm(s2g, 3:4, file = nlin_to_global, clobber = clobber
                , dry = dry)
  invisible(nlin_to_global)
}

#' Apply an xfm to a minc file
#'
#' @param xfm The transform
#' @param input The minc to transform
#' @param output an output file
#' @param like a like file
#' @return the transformed file invisibly
#' @export
apply_xfm <- function(xfm, input, output, like, clobber = TRUE
                    , dry = FALSE){
  
  runner <- `if`(dry, function(x) cat(x, "\n"), system)
  
  if(!file.exists(output) || clobber)
    runner(
      glue("mincresample -transform {xfm} {if(clobber) '-clobber' else ''} "
      , " {input} -like {like} {output}"))

  invisible(output)
}

#' Compute the determinants of nonlinear transform
#' @export
compute_inverse_determinants <-
  function(xfm, type, like, mask
         , output = NULL
         , clobber = TRUE
         , fwhm = 0.2
         , dry = FALSE
         , tempdir = tempfile(pattern = "compute_determinants")
           ){

    if(! type %in% c("abs", "rel"))
      stop("You must specify \"abs\"-olute or \"rel\"-ative jacobians with the `type` "
         , "argument. You specified ", type)
    
    if(is.null(output))
      output <- sub("\\.[^.]$", "", like)

    if(!grepl("^/", output))
      output <- file.path(dirname(xfm), output)

    output_path <- glue("{output}_log_{type}_fwhm{fwhm}.mnc")
    
    if(file.exists(output_path) && !clobber)
      return(invisible(output_path))
    
    cmd <- 
      glue("compute_determinant.py --transform {xfm} "
         , "--inverse --smooth {fwhm} "
         , "--like {like} --mask {mask} "
         , "{if(type == 'rel') '--non-linear-only' else ''} "
         , "--temp-dir {tempdir} "
         , "--log "
         , "--determinant {output_path}"
           )

    runner <-
      `if`(dry
         , function(cmd) cat(cmd, "\n")
         , system)

    runner(cmd)

    invisible(output_path)
  }
