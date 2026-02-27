#Download of GZs

safe_download <- function(url, destfile) {
  if (file.exists(destfile)) {
    return(destfile)
  }
  download.file(url, destfile, mode = "wb")
  return(destfile)
}

library(GEOquery)
library(dplyr)

sub_data <- pData(getGEO("GSE97693")$GSE97693_series_matrix.txt.gz) %>%
  filter(stringr::str_detect(supplementary_file_1, "TPM"))

# all(list.files("GSE97693/TPMs/") == basename(sub_data$supplementary_file_1))

if(all(list.files("GSE97693/TPMs/") == basename(sub_data$supplementary_file_1))){
  print("Todos los archivos han sido descargados.")
}else{
  mapply(
    FUN = safe_download,
    url = sub_data$supplementary_file_1,
    destfile = file.path("GSE97693/TPMs", basename(sub_data$supplementary_file_1))
  )
}


# Extraction of TXTs from .GZs

library(R.utils)

input_dir  <- "../Datos/GSE97693/TPMs"
output_dir <- "GSE97693/TPMs/TXTs"

# Crear la subcarpeta si no existe
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

files <- list.files(input_dir, pattern = "\\.gz$", full.names = TRUE)

for (f in files) {
  # Nombre del archivo sin .gz
  out_file <- file.path(output_dir, sub("\\.gz$", "", basename(f)))
  
  gunzip(f, destname = out_file, remove = FALSE)
}
