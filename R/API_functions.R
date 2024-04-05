#' get canonical SwissProt ID
#'
#' @param gene_name gene name to transform
#' @param organism taxonomy ID
#' @param reviewed (True/False)
#'
#' @return Canonical SwissProtID
#' @export
#' @import httr
#' @examples get_swissprot_id("A2M")
get_swissprot_id <- function(gene_name, organism=9606, reviewed="true") {
  base_url <- "https://rest.uniprot.org/uniprotkb/"
  query <- paste0("search?query=(gene:", gene_name,")+(organism_id:",organism,")+(reviewed:",reviewed,")")
  full_url <- paste0(base_url, query)
  response <- httr::GET(full_url, httr::add_headers("Content-Type" = "application/json", "Accept" = "application/json"))
  if (response$status_code == 200) {
    content <- httr::content(response, "parsed")
    if (length(content$results) != 0) {
      protID <- content$results[[1]]$primaryAccession
    } else {protID <- NA}
    return(protID)
  } else {return(NA)}
}

#' get isoelectric point
#'
#' @param SwissProtID
#'
#' @return Isoelectric Point
#' @export
#' @import httr
#' @importFrom jsonlite fromJSON
#' @examples get_nextprot_PI("P01023")
get_nextprot_PI <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/biophysicochemical-property.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    canonical_isoform <- subset(parsed_data$entry$isoforms, canonicalIsoform == TRUE)
    return(as.numeric(canonical_isoform$isoelectricPointAsString))
  } else {return(NA)}
}

#' get protein mass
#'
#' @param protname UniProtID
#'
#' @return Mass (Da)
#' @export
#' @import httr
#' @importFrom jsonlite fromJSON
#' @examples get_nextprot_mass("P01023")
get_nextprot_mass <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/biophysicochemical-property.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    canonical_isoform <- subset(parsed_data$entry$isoforms, canonicalIsoform == TRUE)
    return(as.integer(canonical_isoform$massAsString))
  } else {return(NA)}
}

#' get number of TM domains
#'
#' @param protein name
#'
#' @return number of TM domains
#' @export
#' @import httr
#' @importFrom jsonlite fromJSON
#' @examples get_nextprot_transmembrane("Q9NRK6")
get_nextprot_transmembrane <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/transmembrane-region.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory$`transmembrane-region`
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    if (is.null(data2mine)) {
      print(0)
      return(0)
    } else {return(nrow(data2mine))}
  } else {return(NA)}
}


#' Get subcellular locations
#'
#' @param protname UniProtID
#'
#' @return subcellular location(s) in strings
#' @export
#' @import httr
#' @importFrom jsonlite fromJSON
#' @examples get_nextprot_subcell("P01023")
get_nextprot_subcell <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/subcellular-location.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory$`subcellular-location`
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    result <- paste(data2mine$cvTermName, collapse=", ")
    return(result)
  } else {return(NA)}
}

#' Transform list of subcell locations into a wide logistic table.
#'
#' @param csv_column the specific column with the comma separated string
#' @param df facultative data.frame to retrieve rownames if needed
#'
#' @return logistic table of the mined strings
#' @export
#'
#' @examples table_maker_from_strings(newdata$Subcell, newdata)
table_maker_from_strings <- function(csv_column, df=NULL) {
  locations_list <- strsplit(csv_column, ", ")
  unique_locations <- unique(unlist(locations_list))
  logical_matrix <- matrix(FALSE, nrow = length(csv_column), ncol = length(unique_locations), dimnames = list(NULL, unique_locations))
  for (i in seq_along(locations_list)) {
    logical_matrix[i, unique_locations %in% locations_list[[i]]] <- TRUE
  }
  logical_df <- as.data.frame(logical_matrix)
  if (!missing(df)) {rownames(logical_df) <- rownames(df)}
  return(logical_df)
}

#' get glycosylation for a given protein
#'
#' @param protname UniProtID
#'
#' @return the type of glycosylation in strings or vector if few
#' @export
#' @import httr
#' @importFrom jsonlite fromJSON
#' @examples get_nextprot_glycosite("Q8N2K0")
get_nextprot_glycosite <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/glycosylation-site.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory$`glycosylation-site`
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    result <- paste(data2mine$cvTermName, collapse=", ")
    return(result)
  } else {return(NA)}
}

#' Get number of coiled-coil of the protein
#'
#' @param protname UniProtID
#'
#' @return number of coiled-coil
#' @export
#'
#' @examples get_nextprot_coiledcoil("Q99996")
get_nextprot_coiledcoil <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/coiled-coil-region.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory$`coiled-coil`
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    if (is.null(data2mine)) {
      return(0)
    } else return(nrow(data2mine))
    #return(result)
  } else {return(NA)}
}

#' Title Get number of lipidation site of the protein
#'
#' @param protname UniProtID
#'
#' @return number of lipidation site
#' @export
#'
#' @examples get_nextprot_lipidationsite("Q14254")
get_nextprot_lipidationsite <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/lipidation-site.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory$`lipidation-site`
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    if (is.null(data2mine)) {
      return(0)
    } else {return(nrow(data2mine))}
    #return(result)
  } else {return(NA)}
}

#' Title Get number of Zinc-Finger region of the protein
#'
#' @param protname UniProtID
#'
#' @return number of zinc-finger region
#' @export
#'
#' @examples get_nextprot_zincfinger("O14686")
get_nextprot_zincfinger <- function(protname) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- "/zinc-finger-region.json"
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory$`zinc-finger-region`
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    if (is.null(data2mine)) {
      return(0)
    } else {return(nrow(data2mine))}
  } else {return(NA)}
}

#' Title General function to mine integer annotations (number of TM domains, catalytic activity...)
#'
#' @param protname UniProtID
#' @param name_of_annotation name of the string (separated with "-")
#'
#' @return the number of annotations
#' @export
#'
#' @examples get_nextprot_annotationsnumber("O15118", "transmembrane-region")
get_nextprot_annotationsnumber <- function(protname, name_of_annotation) {
  # Nextprot API endpoint and ProteinID
  api_endpoint <- "https://api.nextprot.org/entry/"
  api_end <- paste0("/", name_of_annotation, ".json")
  protein_id <- paste0("NX_", protname)
  # Send the API request and retrieve the JSON data
  response <- httr::GET(paste0(api_endpoint, protein_id, api_end))
  if (response$status_code == 200) {
    json_data <- httr::content(response, as = "text")
    # Parse JSON data into an R object
    parsed_data <- (jsonlite::fromJSON(json_data))
    data2mine <- parsed_data$entry$annotationsByCategory[[name_of_annotation]]
    data2mine <- data2mine[data2mine$qualityQualifier == "GOLD",]
    if (is.null(data2mine)) {
      return(0)
    } else {return(nrow(data2mine))}
  } else {return(NA)}
}

#' Create a table traducing gene name into canonical UniProtID
#'
#' @param df containing gene names in first columnx§
#'
#' @return df containing UniProtID
#' @export
#' @import pbapply
#' @examples uniprot_assembling(head(datamined))
uniprot_assembling <- function(df) {
  # first please choose a dataframe containing the gene.names in first column
  data2assemble <- data.frame(name=df[,1], UniProtID=rep(NA, length(df[,1])), stringsAsFactors=FALSE)
  data2assemble$UniProtID <- unlist(pbapply::pblapply(data2assemble$name, function(gene) get_swissprot_id(gene)))
  return(na.omit(data2assemble))
}

#' Incorporate a new column from the datamining functions above
#'
#' @param df containing "UniProtID" as colname
#' @param colname in strings of the datamined column
#' @param fun name of the function to use
#'
#' @return df containing the new column
#' @export
#' @import pbapply
#' @examples datamining_assembling(head(datamined), "PI", get_nextprot_PI)
datamining_assembling <- function(df, colname, fun) {
  # first please choose a dataframe containing "UniProtID" as colname
  data2assemble <- df
  data2assemble$column <- unlist(pbapply::pblapply(data2assemble$UniProtID, function(prot) fun(prot)))
  names(data2assemble)[names(data2assemble) == "column"] <- colname
  return(na.omit(data2assemble))
}

#' Incorporate a new column from the datamining function "get_nextprot_annotationsnumber")
#'
#' @param df containing "UniProtID" as colname
#' @param colname in strings of the datamined column
#' @param fun name of the function to use
#' @param annotation name of the annotation parameter of the function to use
#'
#' @return df containing the new column
#' @export
#' @import pbapply
#' @examples datamining_annotationsnumber(head(datamined), "catalytic_activity")
datamining_annotationsnumber <- function(df, annotation) {
  # first please choose a dataframe containing "UniProtID" as colname
  data2assemble <- df
  data2assemble$column <- unlist(pbapply::pblapply(data2assemble$UniProtID, function(prot) get_nextprot_annotationsnumber(prot, annotation)))
  names(data2assemble)[names(data2assemble) == "column"] <- sub("^([a-z])", "\\U\\1", gsub("-", "_", annotation), perl = TRUE)
  return(na.omit(data2assemble))
}

#' Choose which proteins to analyze from MaxQuant format
#'
#' @param df a df with a dataframe containing one or more protein per row
#' @param colname the column of interest in string
#' @param sep the separator for the proteins
#' @param KeepAllProteins logical. If True, all the proteins will be added to the dataset.
#'
#' @return a dataframe with only one protein per row
#' @export
#'
#' @examples proteinid_parser(Log2_Extraction_comparison_Sep22, "Majority.protein.IDs")
proteinid_parser <- function(df, colname, sep=";", KeepAllProteins=FALSE) {
  if (KeepAllProteins) {
    newdata <- as.data.frame(tidyr::separate_rows(df, colname, sep=sep))
    return(as.data.frame(newdata))
  } else {
    df[[colname]] <- sub(";.*", "", df[[colname]])
    return(as.data.frame(df))
  }
}

#' column parser for analysis
#'
#' @param df the dataframe of interest
#' @param colname2check the column of interest
#' @param newcolname the returned column
#'
#' @return a new dataframe with a new column
#' @export
#'
#' @examples column_maker(datamined, "Glycosite", "Glycological")
column_maker <- function(df, colname2check, newcolname) {
  if (is.character(df[[colname2check]])) {
    df[[newcolname]] <- ifelse(df[[colname2check]]=="", FALSE, TRUE)
  }
  if (is.numeric(df[[colname2check]]) || is.integer(df[[colname2check]])) {
    if (0 %in% df[[colname2check]]) {
      df[[newcolname]] <- as.logical(df[[colname2check]])
    } else {
      mean <- mean(df[[colname2check]])
      df[[newcolname]] <- ifelse(df[[colname2check]]>mean, "High", "Low")
      df[[newcolname]] <- as.factor(df[[newcolname]])
    }
  }
  return(df)
}

#' create a new factor column based on the number of transmembrane domains.
#'
#' @param df with a column containing transmembrane domains
#' @param colname string name of the new column created by the function
#' @param newcolname string of the new colname
#' @param treshold number of transmembrane to seperate. It will create 3 groups, 0, 1-treshold, and >treshold.
#'
#' @return a new dataframe with a new column
#' @export
#'
#' @examples transmembrane_factor(datamined, "TMregions", "TMgroups")
transmembrane_factor <- function(df, colname, newcolname, treshold = 3) {
  middle <- paste0("1-", treshold)
  high <- paste0(">", treshold)
  df[[newcolname]] <- ifelse(df[[colname]] > 0, ifelse(df[[colname]] <= treshold, middle, high), "0")
  df[[newcolname]] <- as.factor(df[[newcolname]])
  return(df)
}

#' Transform glycosylation enrichment as a logisitc table enumrating O & N glycosylation
#'
#' @param table the dataframe containing the glycosylation site from neXtProt
#' @param colname the colname of interest
#'
#' @return a new dataframe, with 2 logistic column O_linked and N_linked
#' @export
#'
#' @examples extract_glyco_string(datamined, "Glycosite")
extract_glyco_string <- function(table, colname) {
  extract_glycosylation <- function(glyco_string) {
    if (grepl("N-linked", glyco_string)) {
      return("N-linked")
    } else if (grepl("O-linked", glyco_string)) {
      return("O-linked")
    } else {
      return(NA)
    }
  }
  glyco_types <- sapply(table[[colname]], extract_glycosylation)
  table$O_linked <- grepl("O-linked", glyco_types)
  table$N_linked <- grepl("N-linked", glyco_types)
  return(table)
}

#' Combine all the membrane terms from the Subcell Table
#'
#' @param df output of table_maker_from_strings(data$Subcell, data)
#'
#' @return the same table but with all the Membrane columns grouped in one
#' @export
#'
#' @examples combineMembraneColumns(subcell)
combineMembraneColumns <- function(df) {
  # Ensure there's a "Membrane" column, initialize as FALSE if it doesn't exist
  if (!"Membrane" %in% names(df)) {
    df$Membrane <- FALSE
  }

  # Identify all columns that contain the word "membrane", excluding the exact "Membrane"
  membrane_cols <- grep("membrane", names(df), ignore.case = TRUE, value = TRUE)
  membrane_cols <- membrane_cols[membrane_cols != "Membrane"]

  # Combine all "membrane" related columns into the "Membrane" column
  for(col in membrane_cols) {
    df$Membrane <- df$Membrane | df[[col]]
  }

  # Drop the original "membrane" columns, keeping the updated "Membrane"
  df <- df[, !(names(df) %in% membrane_cols)]

  return(df)
}
