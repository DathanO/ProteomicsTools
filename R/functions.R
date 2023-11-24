#' Creates MDS Plot from Log2 Protein Abundances to show how the different methods are close to each other
#'
#' @param log2 your first table with log2 normalized data of protein abundances
#' @param rows If True, rownames will take the first column and delete it from the data frame
#' @return A MDS Plot of your methods
#' @export
#' @import edgeR
#' @import limma
#' @examples data("log2")
#' myMDSplot(log2)
myMDSplot <- function(log2, rows = TRUE) {
  if (rows) {
    rownames(log2) <- log2$name
    log2 <- log2[,-1]
  }
  explog <- trunc(2^log2)
  groups <- gsub("_.*", "", colnames(log2))
  y <- edgeR::DGEList(counts=explog, group=groups)
  y <- edgeR::normLibSizes(y)
  return(limma::plotMDS(y))
}

#' Automatize pairwise t.test of your data.
#'
#' @param df A data.frame with protein abundances
#' @param padj A string parameter defining the correction methods. See ?p.adjust
#' @export
#' @return A table of pairwise t.test of 2 times your amount of methods
#' @import stats
#'
#' @examples data("log2")
#' mypairwise_t(log2)
mypairwise_t <- function(df, padj="fdr") {
  df <- replicates_factor(df)
  method <- sort(unique(df$Method))
  methods <- expand.grid(method, method)
  methods <- methods[as.character(methods$Var1) > as.character(methods$Var2),]
  colnames <- unlist(Map(paste, methods$Var1, methods$Var2, sep="_vs_"))
  vectorgene <- unique(df$name)
  ptt <- function(gene) {
    aux <- df[df$name == gene,]
    temp <- na.omit(melt(stats::pairwise.t.test(aux$value, aux$Method, p.adjust.method = padj)$p.value))
    temp$value
  }
  data <- data.frame(Reduce(rbind, Map(ptt, vectorgene)))
  rownames(data) <- vectorgene
  colnames(data) <- colnames
  data
}


#' Automatize and Save VennDigram of the result of the mypairwise_t() above
#'
#' @param output of mypairwise_t()
#' @param path will define directory to store output. by default, files are saved within wd
#' @export
#' @return Nothing is returned, but x plots will be stored, x being the amount of methods in your dataset
#' @import ggvenn
#' @import ggplot2
#'
#' @examples data("ttest")
#' myvenndiagram(ttest)
myvenndiagram <- function(output_ttest, path="examples/") {
  methods <- unique(unlist(strsplit(colnames(output_ttest), "_vs_")))
  pairs <- colnames(output_ttest)
  for (method in methods) {
    listmethod <- list()
    for (pair in pairs) {
      if (grepl(method, pair)) {
        listmethod[[pair]] <- rownames(output_ttest[output_ttest[, pair] < 0.05, ])
      }
    }
    plot <- ggvenn::ggvenn(listmethod, set_name_size = 2)
    plot
    filename <- paste("VennDiagram", method, sep="_")
    filename <- paste(path, filename, ".png", sep="")
    ggplot2::ggsave(filename, bg="white", limitsize = FALSE)
  }
  return(plot)
}

#' Used melt() and separate() to reformat your table into a long one, which will seperate your replicates as a new variable.
#'
#' @param df Unmelted df (log2)
#' @export
#' @return A long and reformated data.frame
#'
#' @examples data("log2")
#' replicates_factor(log2)
replicates_factor <- function(df) {
  longdf <- reshape2::melt(df)
  final <- magrittr::"%>%"(longdf, tidyr::separate(., variable, c("Method", "Replicates"), sep = "_"))
  return(final)
}

#' This function ensure to get the right input for fishertest()
#'
#' @param df A Dataframe that needs to be logical. Numeric value like p-values will be transformed as logical depending to a certain treshold. Integer values will be taken as default (0=false, else=true)
#' @param padj Correction method. See p.adjust()
#' @param treshold The p-value treshold.
#'
#' @return a data.frame with only logical variable instead of integer and numeric variables.
#' @export
#'
#' @examples data("ttest")
#' reformat_table(ttest)
reformat_table <- function(df, padj="fdr", treshold=0.05) {
  for (ncol in 1:length(df)) {
    if (is.integer(df[,ncol]) == TRUE) {
      df[,ncol] <- as.logical(df[,ncol])
    }
    if (is.numeric(df[,ncol]) == TRUE) {
      if (padj != "none") {df[,ncol] <- as.logical(p.adjust(df[,ncol] < treshold, padj))}
      else {df[,ncol] <- as.logical(df[,ncol] < treshold)}
    }
  }
  output <- df
  return(output)
}

#' Makes a Fisher Exact Test to get enrichment of your logical columns
#'
#' @param df_mined_logical Is a dataframe of your data.set but with biophysical properties
#' @param ttest_logical A logical data-frame of the ouput of mypairwise_t().
#' @param padj Correction method for your p.value. See ?p.adjust
#'
#' @return A fisher table with the numbers that were used to make the f.test, and the p-value and p-adjust.
#' @export
#'
#' @examples data("ttest")
#' # taille du points: nombre de gene T_T
#' # couleur du point: p.value
#' data("big")
#' fishertest(reformat_table(big), reformat_table(ttest), "fdr")
fishertest <- function(df_mined, ttest_logical, padj = "none") {
  results <- data.frame(Table = character(), False_False = integer(), False_True = integer(), True_False = integer(), True_True = integer(), p_val_fisher = numeric(), p_adj = numeric(), stringsAsFactors = FALSE)
  for (ncol in (1:length(df_mined))) {
    if (is.logical(df_mined[,ncol]) == TRUE) {
      for (colsign in 1:length(ttest_logical)) {
        conting <- table(ttest_logical[,colsign], df_mined[,ncol])
        ftest <- stats::fisher.test(conting, simulate.p.value = TRUE)
        varname <- paste(list(colnames(df_mined[ncol])), list(colnames(ttest_logical[colsign])), sep="_X_")
        pval <- ftest$p.value
        vec <- cbind(varname, conting[1,1], conting[1,2], conting[2,1], conting[2,2], pval)
        colnames(vec) <- c("Table", "False_False", "False_True", "True_False", "True_True", "pval_fisher")
        results <- rbind(results, vec)
      }
    }
  }
  results$True_True <- as.numeric(results$True_True)
  results$False_True <- as.numeric(results$False_True)
  results$gene.ratio <- as.numeric(results$True_True) / as.numeric((results$True_True + results$False_True))
  results$padj <- p.adjust(results$pval_fisher, padj)
  return(results)
}


#' This function makes the mean of your replicates for each methods
#'
#' @param df with protein abundances
#'
#' @return A data.frame of the mean of each method
#' @export
#'
#' @examples data("log2")
#' mean_function(log2)
mean_function <- function(df) {
  df <- replicates_factor(df)
  meanagg <- aggregate(value ~ name + Method, data = df, FUN = mean)
  meandf <- reshape(meanagg, idvar="name", timevar= "Method", direction="wide")
  names(meandf) <- gsub("value\\.", "", names(meandf))
  rownames(meandf) <- meandf$name
  meandf <- meandf[,-1]
  meandf
}

#' Makes foldchange around your methods
#'
#' @param df with protein abundances in mean. removing replicates can be done with mean_function(). gene name should be named "name" if you want to keep it as rowname.
#'
#' @return A data.frame of the different fold.change values.
#' @export
#'
#' @examples data("log2")
#' myfoldchange(mean_function(log2))
myfoldchange <- function(meandflog2) {
  log2fc <- function(x,y) {
    log2(x/y)
  }
  meandf <- reshape2::melt(meandflog2)
  method <- sort(unique(meandf$variable))
  methods <- expand.grid(method, method)
  methods <- methods[as.character(methods$Var1) > as.character(methods$Var2),]
  colnames <- unlist(Map(paste, methods$Var1, methods$Var2, sep="_vs_"))
  data <- Map(function(x,y) log2fc(meandf[meandf$variable == x, 'value'],
                                   meandf[meandf$variable == y, 'value']),
              methods$Var1, methods$Var2)
  names(data) <- colnames
  data <- as.data.frame(data)
  rownames(data) <- rownames(meandflog2)
  data
}

#' Automatize VolcanoPlots of your dataset
#'
#' @param dfpval A data.frame of p.values (Output of mypairwise_t())
#' @param dffoldchange A data.frame of fold changes (Output of myfoldchange())
#' @param binarydf A reformated data.frame of mined data (reformat_table(big))
#' @param binarycol A String parameter of the column you wish to highlight
#' @param path will define directory to store output. by default, files are saved within working directory
#' @param threshold should depend on the p.value adjustment, you can adapt it to your data
#' @export
#' @return Saves X plots of your pairs of methods, X being 2 times the amount of methods
#' @import ggplot2
#' @importFrom rlang sym
#'
#' @examples data("ttest")
#' data("foldc")
#' data("big")
#' volcano_maker(ttest, foldc, reformat_table(big), "Groups.of.TM")
volcano_maker <- function(dfpval, dffoldchange, binarydf, binarycol, threshold=5, path="examples/") {
  for (i in colnames(dfpval)) {
    new <- as.list(sort(unlist(strsplit(i, "_vs_"))))
    for (l in colnames(dffoldchange)) {
      cut <- as.list(sort(unlist(strsplit(l, "_vs_"))))
      if (setequal(new, cut)) {
        name <- paste(new, collapse="_vs_")
        table2plot <- cbind(dfpval[,i], dffoldchange[,l], binarydf[,binarycol])
        table2plot <- data.frame(table2plot)
        rownames(table2plot) <- rownames(dfpval)
        colnames(table2plot) <- c("log10pval", "log2foldchange", binarycol)
        table2plot[,binarycol] <- as.logical(table2plot[,binarycol])
        plot <- ggplot2::ggplot(data=table2plot, ggplot2::aes(x=log2foldchange, y=-log10(log10pval), color=!!rlang::sym(binarycol))) +
          ggplot2::geom_point() + ggplot2::ylab("-log10pval") +
          ggplot2::geom_vline(xintercept=c(-0.1, 0.1), col="black") +
          ggplot2::geom_hline(yintercept=threshold, col="black") +
          ggplot2::ylim(c(0, 11)) +
          ggplot2::xlim(c(-0.6, 0.6)) +
          ggplot2::ggtitle(rlang::sym(l)) +
          ggplot2::theme_classic()
        filename <- paste("Volcano", binarycol, l, sep="_")
        filename <- paste(path, filename, ".png", sep="")
        plot
        ggplot2::ggsave(filename)
      }
    }
  }
  return(plot)
}


#' Title
#'
#' @param num_df A df with numeric values
#' @param scale_data A logical parameter scaling your data or not. Note that data should be scaled.
#'
#' @return Heatmap of your methods
#' @import ComplexHeatmap
#' @export
#' @examples data("log2")
#' heatmap_maker(log2[,-1])
heatmap_maker <- function(num_df, scale_data=TRUE) {
  mat <- as.matrix(num_df)
  if (scale_data) {
    mat <- t(scale(t(mat)))
    ComplexHeatmap::Heatmap(mat)
  } else {
    ComplexHeatmap::Heatmap(mat)
  }
}

#' Title
#'
#' @param ftest The output of fishertest()
#' @param path will define directory to store output. by default, files are saved within wd
#' @export
#' @return A plot of the fisher test enrichment
#' @import ggplot2
#' @import tidyr
#'
#' @examples data("ftest")
#' fisher_dotplot(ftest)
fisher_dotplot <- function(ftest, path="examples/") {
  results <- ftest[ftest$padj<0.05,]
  result <- magrittr::"%>%"(results, tidyr::separate(., Table, c("Mined", "Pair"), sep = "_X_"))
  plot <- ggplot2::ggplot(result, ggplot2::aes(x=Pair, y=Mined, color=padj, size=gene.ratio)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low="red", high="blue") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust = 1)) +
    ggplot2::xlab("Pairs of Methods") +
    ggplot2::ylab("Mined Properties") +
    ggplot2::ggtitle("DotPlot [Enrichment Results]")
  plot
  filename <- paste(path, "FisherDotPlot_Enrichment.png", sep="")
  ggplot2::ggsave(filename, width = 25, height = 6)
  return(plot)
}

#' Format your table for others functions
#'
#' @param data.frame to be reformated
#' @import reshape2
#' @return a long data.frame (melted)
#' @export
#'
#' @examples data("ttest")
#' col_sorted(ttest)
col_sorted <- function(data) {
  data$name <- rownames(data)
  dataL <- reshape2::melt(data)
  dataL$variable <- as.character(dataL$variable)
  dataL$variable <- strsplit(dataL$variable, "_vs_")
  dataL$variable <- lapply(dataL$variable, sort)
  dataL$variable <- sapply(dataL$variable, paste, collapse="_vs_")
  return(dataL)
}

#' Compare your volcano plots together
#'
#' @param ttest output of mypairwise_t
#' @param foldc output of myfoldchange
#' @param big dataframe with a logical variable to highlight
#' @param binarycol the variable to highlight
#' @param path will define directory to store output. by default, files are saved within wd
#' @param threshold should depend on the p.value adjustment, you can adapt it to your data
#' @import ggplot2
#' @importFrom rlang sym
#' @return Volcanoplots for each pair of methods
#' @export
#'
#' @examples data("ttest")
#' data("foldc")
#' volcano_all(ttest, foldc, big, "Groups.of.TM")
#' # to visualize each variable in a for loop, you can use colnames():
#' # for (colname in colnames(big)[-1]) {volcano_all(ttest, foldc, big, colname)}
volcano_all <- function(ttest, foldc, big, binarycol, threshold=5, path="examples/") {
  ttestL <- col_sorted(ttest)
  foldcL <- col_sorted(foldc)
  table <- merge(ttestL, foldcL, by=c("name", "variable"))
  table2highlight <- subset(big, select=c("name", binarycol))
  table2plot <- merge(table, table2highlight, by="name")
  table2plot[,binarycol] <- as.logical(table2plot[,binarycol])
  plot <- ggplot2::ggplot(data=table2plot, ggplot2::aes(x=value.y, y=-log10(value.x), color=!!rlang::sym(binarycol))) +
    ggplot2::geom_point() + ggplot2::facet_grid(~variable) +
    ggplot2::ylab("-log10pval") +
    ggplot2::xlab("log2fc") +
    ggplot2::geom_vline(xintercept=c(-0.1, 0.1), col="black") +
    ggplot2::geom_hline(yintercept=threshold, col="black") +
    ggplot2::theme_classic()
  plot
  filename <- paste(path, "Volcano_", binarycol, ".png", sep="")
  ggplot2::ggsave(filename, width = 40, height = 6)
  return(plot)
}

