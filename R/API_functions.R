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
  } else {return(paste0("Error: ", response$status_code))}
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
    canonical_isoform$isoelectricPointAsString
  } else {return(paste0("Error: ", response$status_code))}
}

#' get protein mass
#'
#' @param protname
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
    canonical_isoform$massAsString
  } else {return(paste0("Error: ", response$status_code))}
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
    annotationgold <- parsed_data$entry$annotationsByCategory
    sum(annotationgold$`transmembrane-region`$qualityQualifier=="GOLD")
  } else {return(paste0("Error: ", response$status_code))}
}

#' Get subcellular locations
#'
#' @param protname
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
    data2mine$cvTermName
  } else {return(paste0("Error: ", response$status_code))}
}

#' Transform list of subcell locations into a wide logistic table.
#'
#' @param df
#'
#' @return logistic table
#' @export
#'
#' @examples subcell_table_maker(dataf)
subcell_table_maker <- function(df) {
  table <- df
  rownames(table) <- table$Gene
  table <- table["SubcellStrings"]
  SubCellLoc <- unique(unlist(table$SubcellStrings))
  for (locations in SubCellLoc) {
    table[[locations]] <- NA
    table[[locations]] <- sapply(table$SubcellStrings, function(x) locations %in% x)
  }
  table <- table[,-1]
  return(table)
}

#' get glycosylation for a given protein
#'
#' @param protname
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
    data2mine$cvTermName
  } else {return(paste0("Error: ", response$status_code))}
}


#' Create a table traducing gene name into canonical UniProtID
#'
#' @param df containing gene names
#'
#' @return df containing UniProtID
#' @export
#' @import pbapply
#' @examples uniprot_assembling(big)
uniprot_assembling <- function(df) {
  # first please choose a dataframe containing the gene.names
  data2assemble <- data.frame(Gene=df[,1], UniProtID=rep(NA, length(df[,1])), stringsAsFactors=FALSE)
  data2assemble$UniProtID <- pbapply::pblapply(data2assemble$Gene, function(gene) get_swissprot_id(gene))
  return(data2assemble)
}

#' Incorporate a new column from the datamining functions above
#'
#' @param df containing UniProtID
#' @param colname in strings of the new column
#' @param fun name of the function to use
#'
#' @return df containing the new column
#' @export
#' @import pbapply
#' @examples datamining_assembling(test, "PI", get_nextprot_PI)
datamining_assembling <- function(df, colname, fun) {
  # first please choose a dataframe containing "UniProtID"
  data2assemble <- df
  data2assemble$column <- pbapply::pblapply(data2assemble$UniProtID, function(prot) fun(prot))
  names(data2assemble)[names(data2assemble) == "column"] <- colname
  return(data2assemble)
}
