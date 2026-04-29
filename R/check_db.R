#' Check neuroimaGene database downloaded
#'
#' Check if the NeuroimaGene database exists in the proper location prior to running the query and prompt user to download if not.
#' @keywords installation
#' @param timeout time to spend downloading the NeuroimaGene database in seconds (default = 900)
#' @param localdb path to local copy of NeuroimaGenefast.db file (default = NA)
#' @export
#' @importFrom utils download.file
#' @importFrom tools md5sum
#' @returns no return value, called to give information on status of neuroimaGene
#' database and prompt user the user to download if resource file is missing.
#' @examples
#' check_db(timeout = 600)



# Define the function to check and download the database
check_db <- function(timeout=900, localdb = NA) {
  
  pkg_dir <- system.file(package = "neuroimaGene")
  db_path <- file.path(pkg_dir, "extdata", "NeuroimaGenefast.db")
  db_url <- "https://zenodo.org/records/10994978/files/NeuroimaGenefast.db"
  
  
  #Download local copy of neuroimaGene SQL database as specified by user
  if (!is.na(localdb)){    
    if (tools::md5sum(localdb) == '60dab0d523b046f242c0c0f34d918b76') { #check md5 of local file
      message("md5 of local db file passes internal check. integrating file into package data")
      file.rename(from = localdb, to = db_path)
      if (tools::md5sum(db_path) == '60dab0d523b046f242c0c0f34d918b76') {
        message("neuroimaGene SQL database integrated into package data")
      } else {
        message("neuroimaGene SQL database integration was unsuccessful. See documentation for troubleshooting.")
      }
    } else { 
      message("md5 of local db file does not pass internal check. neuroimaGene SQL file may be corrupted or incorrect file")
      response <- readline(prompt = "Would you like to proceed using specified file? (yes/no): ")
      if (tolower(response) == "yes") {
        file.rename(from = localdb, to = db_path)
      } else {
        message("neuroimaGene SQL database integration was unsuccessful. See documentation for troubleshooting.")
      }
    }
  } else if (!file.exists(db_path)) { # If the file does not exist, ask about proceeding with download
    # Prompt the user for permission to download the database
    response <- readline(prompt = "NeuroimaGene database not found. It will require ~1.9Gb of space. Do you want to download it from Zenodo? (yes/no): ")
    if (tolower(response) == "yes") {
      message("Downloading database from Zenodo...\n(This step may take >10-15 minutes depending on download speed)")
      message(paste("user determined timeout:", timeout, "seconds"))
      options(timeout = timeout)  # Set timeout for download
      download_success <- FALSE
      tryCatch({
        if (Sys.info()["sysname"] == "Windows") {
          download.file(db_url, db_path, mode = "wb", method = "wininet", quiet = FALSE)
        } else {
          download.file(db_url, db_path, mode = "wb", method = "libcurl", quiet = FALSE)
        }
        if (tools::md5sum(db_path) == '60dab0d523b046f242c0c0f34d918b76') {
          download_success <- TRUE
        }else { file.remove(db_path)
          message("neuroimaGene SQL file was corrupted.")
        }
      }, error = function(e) {
        file.remove(db_path)
        message("Download failed: ", e$message)
      })
      if (!download_success) {
        stop("Database download failed. Please see vignette/check_db() documentation for further options.")
      } else {
        message("Download complete.")
      }
    } else {
      message("Database download was not permitted by the user.")
    }
  } else {
    # Prompt the user for permission to re-download the database
    response <- readline(prompt = "Database already exists locally. Would you like to delete current version and redownload from Zenodo? It will require ~1.9Gb of space. (yes/no): ")
    if (tolower(response) == "yes") {
      message("Deleting current database")
      file.remove(db_path)
      message("Downloading database from Zenodo...\n(This step may take >10-15 minutes depending on download speed)")
      message(paste("user determined timeout:", timeout, "seconds"))
      options(timeout = timeout)  # Set timeout for download
      download_success <- FALSE
      tryCatch({
        if (Sys.info()["sysname"] == "Windows") {
          download.file(db_url, db_path, mode = "wb", method = "wininet", quiet = FALSE)
        } else {
          download.file(db_url, db_path, mode = "wb", method = "libcurl", quiet = FALSE)
        }
        if (tools::md5sum(db_path) == '60dab0d523b046f242c0c0f34d918b76') {
          download_success <- TRUE
        } else { message("neuroimaGene SQL file was corrupted")
        }
      }, error = function(e) {
        file.remove(db_path)
        message("Download failed: ", e$message)
      })
      if (!download_success) {
        stop("Database download failed. Please see vignette/check_db() documentation for further options.")
      } else {
        message("Download complete.")
      }
    } else {
      message("Database download was not permitted by the user.")
    }
  }
}
