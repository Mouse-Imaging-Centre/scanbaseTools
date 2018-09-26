#' @param db The google sheet title.
#' @export
read_scanbase <- function(db = "test_scanbase_40um"){
  sbase <- googlesheets::gs_title(db)
  tfile <- tempfile(fileext = ".xlsx")
  gs_download(sbase, to = tfile)
  list(scans = readxl::read_excel(tfile, "Scans")
     , studies = readxl::read_excel(tfile, "Studies")
     , genotypes = readxl::read_excel(tfile, "Genotypes")
     , treatments = readxl::read_excel(tfile, "Treatments")
     , global = readxl::read_excel(tfile, "Global")
       )
}

#' @export
load_pydpiper_results <-
  function(ppd
         , common = "/hpf/largeprojects/MICe/scanbase/pipeline-40um/scanbase_second_level_nlin/scanbase_second_level-nlin-3.mnc"
         , clobber = FALSE){
    
    transforms <- read.csv(file.path(ppd, "transforms.csv")
                         , stringsAsFactors = FALSE)

    determinants <-
      read.csv(file.path(ppd, "determinants.csv")
             , stringsAsFactors = FALSE) %>%
      filter(fwhm == 0.2)

    column_mapping <-
      c(
        Distortion_Corrected_Scan = "native_file"
      , Scan_To_Study_Absolute_Jacobians = "log_full_det"
      , Scan_To_Study_Relative_Jacobians = "log_nlin_det"
      , Scan_To_Global_Absolute_Jacobians = "log_full_det_common"
      , Scan_To_Global_Relative_Jacobians = "log_nlin_det_common"
      , Scan_To_Study_Global_Space_Resampled_Absolute_Jacobians = "unknown"
      , Scan_To_Study_Global_Space_Resampled_Relative_Jacobians = "unknown"
      , Rigid_Transform = "rigid_xfm"
      , Rigid_Filepath = "lsq6_file"
      , Scan_To_Study_Transform = "lsq12_nlin_xfm"
      , Labels = "unknown"
      )

    known_columns <- discard(column_mapping, ~ . == "unknown")

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
           , Labels = map_chr(file.path(Processed_dir, "voted.mnc")
                            , function(f){
                              if(file.exists(f)){
                                return(f)
                              } else {
                                stop("these subjects have not been processed with MAGeT")
                              }
                            }))
    

    apply_xfm_wrapper <- function(xfm, input, suffix, like, clobber = TRUE){
      dir <- dirname(input)
      subject <- basename(dirname(dir)) ## Subject processed dir
      output <- file.path(dir, paste0(subject, "_", suffix))
      
      apply_xfm(xfm, input, output, like, clobber = clobber)
    }

    ## FIX
    scans <-
      full_data %>%
      mutate(
        xfm_nlin3_to_global =
          map_chr(overall_xfm_to_common, ~ global_from_concat(., "nlin3-to-global.xfm"
                                                            , clobber = clobber))
      , Scan_To_Study_Global_Space_Resampled_Absolute_Jacobians =
          future_map2_chr(xfm_nlin3_to_global, Scan_To_Study_Absolute_Jacobians
                 , ~ apply_xfm_wrapper(.x, .y, "scan_to_study_abs_jacobians_resampled_global.mnc"
                                     , like = common
                                     , clobber = clobber))
      , Scan_To_Study_Global_Space_Resampled_Relative_Jacobians =
          future_map2_chr(xfm_nlin3_to_global, Scan_To_Study_Relative_Jacobians
                 , ~ apply_xfm_wrapper(.x, .y, "scan_to_study_rel_jacobians_resampled_global.mnc"
                                     , like = common
                                     , clobber = clobber))
      ) %>%
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

      files <-
        xfm %>%
        scan(what = character()) %>%
        grep("*.mnc", ., value = TRUE) %>%
        file.path(ppd, .) %>%
        c(xfm, .)

      print(files)
      new_files <- file.path(ppd, basename(files))
      
      walk2(files, new_files, ~ {
        if(!file.exists(.y) || clobber)
          file.copy(.x, .y)
      })

      study_to_global <- "study_to_global.xfm"
      global_from_concat(new_files[1], study_to_global)

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
                       , clobber = FALSE ){

  sb <- read_scanbase(db)
  pp <- load_pydpiper_results(ppd, common, clobber)
  
  scanf <- read.csv(scan, stringsAsFactors = FALSE)
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

  if(nrow(scanf_merged) != nrow(pp$scans))
    stop("Scans present in the pydpiper directory were not merged successfully, "
       , "this potentially indicates a problem with you scan metadata")

  if(nrow(scanf_merged) != nrow(scanf))
    stop("Scans present in the scan metadata were not in your pydpiper results, "
       , "this potentially indicates a problem with you scan metadata")

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
       
  sbase <- googlesheets::gs_title(db)
  token <- googlesheets::gs_token()

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
               , "{sbase$sheet_key}/"
               , "values/{range_ref}?valueInputOption=USER_ENTERED")
          , httr::config(token = token)
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
                           , clobber = TRUE){
  header <- c(unc_xfms$header %>% { .[!grepl("^ *$", .)] }
            , glue("%{Sys.Date()}>>> Unconcat in R and reassembled from pieces:",
                   " {paste0(which, collapse = ', ')}"))

  xfm <- c(header, "", unlist(unc_xfms$xfms[which]))
  if(!is.null(file)){
    if(!file.exists(file) || clobber){
      cat(xfm, file = file, sep = "\n")
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
#' @return the path to the output file invisibly.
#' @export
global_from_concat <- function(scan_to_global, nlin_to_global
                             , clobber = TRUE){
  dir <- dirname(scan_to_global)
  if(grepl("/", nlin_to_global)) stop("nlin_to_global cannot be a path")
  nlin_to_global <- file.path(dir, nlin_to_global)
  
  s2g <- unconcat_xfm(scan_to_global)
  reconstruct_xfm(s2g, 3:4, file = nlin_to_global, clobber = clobber)
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
apply_xfm <- function(xfm, input, output, like, clobber = TRUE){
  if(!file.exists(output) || clobber)
    system(
      glue("mincresample -transform {xfm} -clobber {input} -like {like} {output}"))

  invisible(output)
}

compute_determinants <-
  function(xfm, output, mask, clobber = TRUE){

  }
