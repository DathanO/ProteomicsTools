#' Creates MDS Plot from Log2 Protein Abundances to show how the different methods are close to each other
#'
#' @param log2 your first table with log2 normalized data of protein abundances, with a "name" column of your identifiers
#' @param rows If True, rownames will take the first column and delete it from the data frame
#' @return A MDS Plot of your methods
#' @export
#' @import edgeR
#' @import limma
#' @examples myMDSplot(log2)
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
#' @param df A data.frame with protein abundances, with a "name" column of your identifiers
#' @param padj A string parameter defining the correction methods. See ?p.adjust
#' @export
#' @return A table of pairwise t.test of 2 times your amount of methods
#' @import stats
#'
#' @examples data("log2")
#' mypairwise_t(head(log2))
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
myvenndiagram <- function(output_ttest, path="", saveplot=FALSE) {
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
    if (saveplot) {
      filename <- paste("VennDiagram", method, sep="_")
      filename <- paste(path, filename, ".png", sep="")
      ggplot2::ggsave(filename, bg="white", limitsize = FALSE)
    }
  }
  return(plot)
}

#' Used melt() and separate() to reformat your table into a long one, which will seperate your replicates as a new variable. Default separator is "_".
#'
#' @param df Unmelted df (log2)
#' @export
#' @return A long and reformated data.frame
#'
#' @examples data("log2")
#' replicates_factor(log2)
replicates_factor <- function(df, sep="_") {
  longdf <- reshape2::melt(df)
  final <- magrittr::"%>%"(longdf, tidyr::separate(., variable, c("Method", "Replicates"), sep = sep))
  final$Method <- as.factor(final$Method)
  final$Replicates <- as.factor(final$Replicates)
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
#' data("datamined")
#' fishertest(reformat_table(datamined), reformat_table(ttest), "fdr")
fishertest <- function(df_mined, ttest_logical, padj = "none") {
  results <- data.frame(Table = character(), False_False = integer(), False_True = integer(), True_False = integer(), True_True = integer(), p_val_fisher = numeric(), p_adj = numeric(), stringsAsFactors = FALSE)
  diff <- setdiff(ttest_logical$name, df_mined$name)
  for (na in diff) {ttest_logical <- ttest_logical[ttest_logical$name!=na,]}
  for (ncol in (1:length(df_mined))) {
    if (is.logical(df_mined[,ncol]) == TRUE) {
      for (colsign in 1:length(ttest_logical)) {
        conting <- table(ttest_logical[,colsign], df_mined[,ncol])
        if (all(dim(conting) == c(2, 2))) {
          ftest <- stats::fisher.test(conting, simulate.p.value = TRUE)
          varname <- paste(list(colnames(df_mined[ncol])), list(colnames(ttest_logical[colsign])), sep="_X_")
          pval <- ftest$p.value
          vec <- cbind(varname, conting[1,1], conting[1,2], conting[2,1], conting[2,2], pval)
          colnames(vec) <- c("Table", "False_False", "False_True", "True_False", "True_True", "pval_fisher")
          results <- rbind(results, vec)
        }
      }
    }
  }
  results$True_True <- as.numeric(results$True_True)
  results$False_True <- as.numeric(results$False_True)
  results$Protein.Ratio <- as.numeric(results$True_True) / as.numeric((results$True_True + results$False_True))
  results$padj <- p.adjust(results$pval_fisher, padj)
  return(results)
}


#' This function makes the mean of your replicates for each methods
#'
#' @param df with protein abundances, with a "name" column of your identifiers
#'
#' @return A data.frame of the mean of each method
#' @export
#'
#' @examples data("log2")
#' mean_function(head(log2))
mean_function <- function(df) {
  dfrep <- replicates_factor(df)
  df$order <- seq_len(nrow(df))
  meanagg <- aggregate(value ~ name + Method, data = dfrep, FUN = mean)
  meandf <- reshape(meanagg, idvar="name", timevar= "Method", direction="wide")
  names(meandf) <- gsub("value\\.", "", names(meandf))
  columns <- subset(df, select=c("name", "order"))
  meandf <- merge(meandf, columns, by="name")
  meandf <- meandf[order(meandf$order),]
  rownames(meandf) <- meandf$name
  meandf <- subset(meandf, select=-(order))
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
#' myfoldchange(head(log2))
myfoldchange <- function(log2) {
  meandflog2 <- mean_function(log2)
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

#' Compute fold-change for each method vs all other methods
#'
#' @param log2 protein abundance table (log2 normalized)
#'
#' @return foldchange matrix for each method
#' @export
#'
#' @examples myfoldchange_unpaired(head(log2))
myfoldchange_unpaired <- function(log2) {
  meandflog2 <- mean_function(log2)
  meandf <- reshape2::melt(meandflog2, id.vars = "name")
  fold_changes_list <- list(ID = meandf$name[meandf$variable == unique(meandf$variable)[1]])
  methods <- unique(meandf$variable)
  for (method in methods) {
    current_method_data <- meandf[meandf$variable == method,]
    other_methods_data <- meandf[meandf$variable != method,]
    other_means <- aggregate(value ~ name, data = other_methods_data, mean)
    merged_data <- merge(current_method_data, other_means, by = "name")
    log2fc <- function(x, y) log2(x / y)
    fold_changes_list[[method]] <- sapply(seq_along(merged_data$value.x),
                                          function(i) log2fc(merged_data$value.x[i], merged_data$value.y[i]))
  }
  fold_changes_df <- data.frame(fold_changes_list)
  rownames(fold_changes_df) <- fold_changes_df$ID
  colnames(fold_changes_df) <- colnames(meandflog2)
  fold_changes_df <- fold_changes_df[,-1]
  return(fold_changes_df)
}

#' Automatize VolcanoPlots of your dataset
#'
#' @param ttest A data.frame of p.values (Output of mypairwise_t())
#' @param foldc A data.frame of fold changes (Output of myfoldchange())
#' @param datamined A reformated data.frame of mined data (reformat_table(datamined))
#' @param binarycol A String parameter of the column you wish to highlight
#' @param path will define directory to store output. by default, files are saved within working directory
#' @param threshold should depend on the p.value adjustment, you can adapt it to your data
#' @param saveplot boolean parameter to whether save or not the plot
#' @export
#' @return Saves X plots of your pairs of methods, X being 2 times the amount of methods
#' @import ggplot2
#' @importFrom rlang sym
#'
#' @examples data("ttest")
#' data("foldc")
#' data("datamined")
#' volcano_maker(ttest, foldc, reformat_table(datamined), "Groups.of.TM")
volcano_maker <- function(ttest, foldc, datamined, binarycol, threshold=5, namefile="Volcano", path="examples/", width=10, height=10, saveplot=FALSE) {
  for (i in colnames(ttest)) {
    new <- as.list(sort(unlist(strsplit(i, "_vs_"))))
    for (l in colnames(foldc)) {
      cut <- as.list(sort(unlist(strsplit(l, "_vs_"))))
      if (setequal(new, cut)) {
        name <- paste(new, collapse="_vs_")
        table2plot <- cbind(ttest[,i], foldc[,l], datamined[,binarycol])
        table2plot <- data.frame(table2plot)
        rownames(table2plot) <- rownames(ttest)
        colnames(table2plot) <- c("log10pval", "log2foldchange", binarycol)
        table2plot[,binarycol] <- as.factor(datamined[,binarycol])
        plot <- ggplot2::ggplot(data=table2plot, ggplot2::aes(x=log2foldchange, y=-log10(log10pval), color=!!rlang::sym(binarycol))) +
          ggplot2::geom_point() + ggplot2::ylab("-log10pval") +
          ggplot2::geom_vline(xintercept=c(-0.1, 0.1), col="black") +
          ggplot2::geom_hline(yintercept=threshold, col="black") +
          ggplot2::ylim(c(0, 11)) +
          ggplot2::xlim(c(-0.6, 0.6)) +
          ggplot2::ggtitle(rlang::sym(l)) +
          ggplot2::theme_classic()
        plot
        if (saveplot) {
          filename <- paste(namefile, binarycol, l, sep="_")
          filename <- paste(path, filename, ".png", sep="")
          ggplot2::ggsave(filename, width = width, height = height)
        }
      }
    }
  }
  return(plot)
}


#' Make Heatmap of your protein abundances table
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

#' Enrichment plot of your data
#'
#' @param ftest The output of fishertest()
#' @param path will define directory to store output. by default, files are saved within wd
#' @export
#' @return A plot of the fisher test enrichment
#' @import ggplot2
#' @import tidyr
#' @examples data("ftest")
#' dotplot_maker(fisher_subcell, width=25, height=10)
dotplot_maker <- function(ftest, namefile="FisherDotPlot_Enrichment.png", path="", width = 25, height = 6, saveplot=FALSE) {
  results <- ftest[ftest$padj<0.05,]
  result <- magrittr::"%>%"(results, tidyr::separate(., Table, c("Mined", "Pair"), sep = "_X_"))
  plot <- ggplot2::ggplot(result, ggplot2::aes(x=Pair, y=Mined, color=padj, size=Protein.Ratio)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(name="P-Value adjusted", low="red", high="blue") +
    ggplot2::scale_size(name="Protein Ratio\nsignificative\nto both T-Test\nand Mined Property", guide=guide_legend(title.theme=element_text(size=8))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust = 1)) +
    ggplot2::xlab("Methods") +
    ggplot2::ylab("Mined Properties") +
    ggplot2::ggtitle("DotPlot [Enrichment Results]") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#BFD5E3", colour = "#6D9EC1"),
                   panel.grid.major = ggplot2::element_line(linewidth = 0.15, linetype = 'solid', colour = 'white'))
  plot
  if (saveplot) {
    filename <- paste(path, namefile, sep="")
    ggplot2::ggsave(filename, width = width, height = height)
  }
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
#' @param ttest output of mypairwise_t, with a "name" column of your identifiers
#' @param foldc output of myfoldchange, with a "name" column of your identifiers
#' @param datamined dataframe with a logical variable to highlight, with a "name" column of your identifiers
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
#' volcano_all(ttest, foldc, datamined, "Groups.of.TM")
volcano_all <- function(ttest, foldc, datamined, binarycol, threshold=5, namefile="Volcano_", path="", width=40, height=6, saveplot=FALSE) {
  ttestL <- col_sorted(ttest)
  foldcL <- col_sorted(foldc)
  table <- merge(ttestL, foldcL, by=c("name", "variable"))
  table2highlight <- subset(datamined, select=c("name", binarycol))
  table2plot <- merge(table, table2highlight, by="name")
  plot <- ggplot2::ggplot(data=table2plot, ggplot2::aes(x=value.y, y=-log10(value.x), color=!!rlang::sym(binarycol))) +
    ggplot2::geom_point() + ggplot2::facet_grid(~variable) +
    ggplot2::ylab("-log10pval") +
    ggplot2::xlab("log2fc") +
    ggplot2::geom_vline(xintercept=c(-0.1, 0.1), col="black") +
    ggplot2::geom_hline(yintercept=threshold, col="black") +
    ggplot2::theme_classic()
  plot
  if (saveplot) {
    filename <- paste(path, namefile, binarycol, ".png", sep="")
    ggplot2::ggsave(filename, width = width, height = height)
  }
  return(plot)
}

#' ggpairs plot with extra highlight on a specific column
#'
#' @param log2 log2 protein abundances table; with "name" column which is the protein/gene identifier
#' @param datamined datamined table, with the exact same column "name"
#' @param column2highlight factor or logical column to highlight plots
#'
#' @return exploratory plots from ggpairs
#' @export
#'
#' @examples ggpairs_proteomics(head(log2), head(datamined), "Groups.of.TM")
ggpairs_proteomics <- function(log2, datamined, column2highlight) {
  meanmethods <- mean_function(log2)
  data2highlight <- subset(datamined, select=c("name", column2highlight))
  data2plot <- merge(meanmethods, data2highlight, by="name")
  data2plot$name <- NULL
  plottitle <- paste0("Exploratory graphs of the data, colored by: ", column2highlight)
  GGally::ggpairs(data2plot, ggplot2::aes(color=!!rlang::sym(column2highlight), alpha = 0.5), upper = list(continuous = GGally::wrap("cor", size = 4)),
                  title=plottitle) + ggplot2::theme(axis.text = element_text(size = 10))
}

#' create a dotplot colored by Zscore and sized by the percent of abundance
#'
#' @param logicaltable a dataframe with mined data contained as logical variables; with protein/gene name as "name"
#' @param log2 a dataframe with protein abundances; with protein/gene name as "name"
#' @param saveplot if True, the plot will be saved
#' @param namefile the name of the plot to save if saveplot=True
#' @param path the path of the output to save if saveplot=True
#' @param width the width of the plot to save if saveplot=True
#' @param height the height of the plot to save if saveplot=True
#'
#' @return dotplot; saved or not
#' @export
#'
#' @examples dotplot_by_methods(logical_mined_data, log2)
dotplot_by_methods <- function(logicaltable, log2, mean_computing=TRUE, saveplot=FALSE, namefile="DotPlot.png", path="", width = 25, height = 6) {
  table <- melt(mean_function(log2))
  table$Zscore <- ave(table$value, table$variable, FUN = function(x) ((x-mean(x))/sd(x)))
  table$percent <- ave(table$value, table$variable, FUN = function(x) ((x - min(x)) / (max(x) - min(x))) * 100)
  logical2plot <- melt(logicaltable, "name")
  tableforplot <- merge(logical2plot, table, by="name")
  tableforplot <- tableforplot[tableforplot$value.x==TRUE,]
  tableforplot <- subset(tableforplot, select = -value.y)
  table2plot <- unique(within(tableforplot, {
    Zscore <- ave(Zscore, variable.y, variable.x, FUN = mean)
    percent <- ave(percent, variable.y, variable.x, FUN = mean)
  }))
  if (mean_computing==FALSE) {table2plot <- tableforplot}
  plot <- ggplot2::ggplot(table2plot, ggplot2::aes(x=variable.x, y=variable.y, color=Zscore, size=percent)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low="red", high="blue") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust = 1)) +
    ggplot2::xlab("Variables") +
    ggplot2::ylab("Methods") +
    ggplot2::ggtitle("DotPlot by Zscore and abundance percent") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#BFD5E3", colour = "#6D9EC1"),
                   panel.grid.major = ggplot2::element_line(linewidth = 0.15, linetype = 'solid', colour = 'white'))
  if (saveplot) {
    filename <- paste(path, namefile, sep="")
    ggplot2::ggsave(filename, width = width, height = height, limitsize = FALSE)
  }
  plot
}

#' Compute T-Test of one method against all others per protein
#'
#' @param df log2 protein abundance matrix
#' @param padj p.adjustment method
#'
#' @return ttest matrix of each method
#' @export
#' @import stats
#' @examples compute_ttest(head(log2))
compute_ttest <- function(df, padj="fdr") {
  df <- replicates_factor(df)
  vectormethod <- sort(unique(df$Method))
  vectorgene <- unique(df$name)
  compare_one_method <- function(method) {
    pergene <- function(gene) {
      df <- df[df$name == gene,]
      t_tests <- t.test(data=df, value ~ Method == method)$p.value
      p_adjusted <- p.adjust(t_tests, method = padj)
      return(p_adjusted)
    }
    results <- Reduce(rbind, Map(pergene, vectorgene))
    return(results)
  }
  results_df <- data.frame(Reduce(cbind, Map(compare_one_method, vectormethod)))
  rownames(results_df) <- vectorgene
  colnames(results_df) <- vectormethod
  return(results_df)
}


#' Compare the mean of one method against all others per protein
#'
#' @param df log2 protein abundance matrix
#' @param padj p.adjustment method
#'
#' @return mean differences matrix
#' @export
#'
#' @examples compute_meandiff_methods(head(log2))
compute_meandiff_methods <- function(df) {
  df <- replicates_factor(df)
  vectormethod <- sort(unique(df$Method))
  vectorgene <- unique(df$name)
  compare_one_method <- function(method) {
    pergene <- function(gene) {
      df_subset <- df[df$name == gene,]
      current_method_data <- df_subset[df_subset$Method == method, ]
      other_methods_data <- df_subset[df_subset$Method != method, ]
      mean_difference <- mean(current_method_data$value) - mean(other_methods_data$value)
      return(mean_difference)
    }
    results <- Reduce(rbind, Map(pergene, vectorgene))
    names(results) <- c("mean_difference")
    return(results)
  }
  results_df <- data.frame(Reduce(cbind, Map(compare_one_method, vectormethod)))
  rownames(results_df) <- vectorgene
  colnames(results_df) <- vectormethod
  return(results_df)
}


#' Create a Plot of the difference of the means after using T-Test and MeanDiff
#'
#' @param tablettest output of compute_ttest()
#' @param tablemean output of mean_diff_methods()
#' @param minedproperties logical table of mined properties
#' @param saveplot True/False to save the plot in the desired path
#' @param namefile name of the output file if saved
#' @param path path of the output file if saved
#' @param width width of the output file if saved
#' @param height height of the output file if saved
#'
#' @return plot
#' @export
#' @import ggplot2
#' @import reshape2
#' @examples dotplot_meandiff(meantt, meandiff, glycolink, treshold=0.0001)
dotplot_meandiff <- function(tablettest, tablemean, minedproperties, treshold=0.05, savetable=TRUE, nametable="SignificativeProteins.csv", saveplot=FALSE, nameplot="MeanDifferenceDotPlot.png", path="", width=10, height=10) {
  tablettest$name <- rownames(tablettest)
  tablemean$name <- rownames(tablemean)
  tablettest <- reshape2::melt(tablettest, id.vars = "name")
  tablemean <- reshape2::melt(tablemean, id.vars = "name")
  minedproperties <- reshape2::melt(minedproperties[,sort(names(minedproperties))], id.vars = "name")
  firstmerge <- merge(tablettest, tablemean, by = c("name", "variable"))
  table2plot <- merge(firstmerge, minedproperties, by = "name")
  table2plot <- table2plot[table2plot$value==TRUE,]
  table2plot <- table2plot[table2plot$value.x < treshold,]
  table2plot <- within(table2plot, {
    MeanDiff <- ave(value.y, variable.x, variable.y)
    MeanDiffBinary <- as.factor(ifelse(MeanDiff<0, "Down", "Up"))
  })
  if (savetable) {write.csv(table2plot, file=nametable, row.names = FALSE)}
  plot <- ggplot2::ggplot(table2plot, ggplot2::aes(x=variable.x, y=variable.y)) +
    ggplot2::geom_point(ggplot2::aes(size=abs(MeanDiff), shape=MeanDiffBinary, fill=value.x), stroke=0.5) +
    ggplot2::scale_fill_continuous(name="Adjusted P-Value", low="red", high="yellow") +
    ggplot2::scale_size(name="Mean Differences\n(Absolute Value)") +
    ggplot2::scale_shape_manual(name="Mean Differences\nabove or below 0", values=c("Up"=24, "Down"=25)) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, vjust=1, hjust = 1)) +
    ggplot2::xlab("Variables") +
    ggplot2::ylab("Methods") +
    ggplot2::guides(shape=ggplot2::guide_legend(override.aes=list(size=3, stroke=1.5))) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#BFD5E3", colour = "#6D9EC1"),
                   legend.text = ggplot2::element_text(size=12),
                   panel.grid.major = ggplot2::element_line(linewidth = 0.15, linetype = 'solid', colour = 'white'))
  plot
  if (saveplot) {
    filename <- paste(path, namefile, ".png", sep="")
    ggplot2::ggsave(filename, width = width, height = height)
  }
  return(plot)
}

#' Create a simple dotplot around the max and min values per method
#'
#' @param logicaltable logical table of the mined properties
#' @param log2 log2 protein abundances matrix
#' @param mean_computing True/False; if True, you will get average values to get one color for the dot. If False, applies a gradient within each dot based on values.
#' @param saveplot True/False to save the plot in the desired path
#' @param namefile name of the output file if saved
#' @param path path of the output file if saved
#' @param width width of the output file if saved
#' @param height height of the output file if saved
#'
#' @return dotplot around the max and min values per method
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr summarise group_by
#'
#' @examples dotratio(subcell, log2)
dotratio <- function(logicaltable, log2, mean_computing=TRUE, saveplot=FALSE, namefile="DotPlot.png", path="", width = 25, height = 6) {
  table <- replicates_factor(log2)
  table$ratio <- ave(table$value, table$Method, FUN = function(x) ((x - min(x)) / (max(x) - min(x))) * 100)
  logical2plot <- reshape2::melt(logicaltable, "name")
  colnames(logical2plot) <- c("name", "annotations", "logical")
  tableforplot <- merge(logical2plot, table, by="name")
  tableforplot <- tableforplot[tableforplot$logical == TRUE,]
  range_values <- dplyr::summarise(
    dplyr::group_by(tableforplot, annotations),
    overall_max = max(value),
    overall_min = min(value),
    .groups = 'drop'
  )
  method_max <- dplyr::summarise(
    dplyr::group_by(tableforplot, annotations, Method),
    method_max = max(value),
    .groups = 'drop'
  )
  mean_ratios <- dplyr::summarise(
    dplyr::group_by(tableforplot, annotations, Method),
    mean_percent_ratio = mean(ratio),
    .groups = 'drop'
  )
  merged_data <- merge(method_max, range_values, by = "annotations")
  merged_data <- merge(merged_data, mean_ratios, by = c("annotations", "Method"))
  merged_data <- within(merged_data, {
    normalized_proportion <- (method_max - overall_min) / (overall_max - overall_min)
    Zscore <- ave(mean_percent_ratio, annotations, FUN = function(x) ((x - mean(x)) / sd(x)))
  })
  plot <- ggplot2::ggplot(merged_data, ggplot2::aes(x = annotations, y = Method, size = normalized_proportion, color = Zscore)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "yellow", high = "red") +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggplot2::xlab("Variables") +
    ggplot2::ylab("Methods") +
    ggplot2::ggtitle("DotPlot by Zscore and abundance ratio") +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "#BFD5E3", colour = "#6D9EC1"),
                   panel.grid.major = ggplot2::element_line(linewidth = 0.15, linetype = 'solid', colour = 'white'))
  if (saveplot) {
    filename <- paste(path, namefile, sep = "")
    ggplot2::ggsave(filename, plot = plot, width = width, height = height, limitsize = FALSE)
  }
  return(plot)
}


#' Compute anova analysis for each protein to compare each method altogether
#'
#' @param log2 protein abundance matrix (log2 normalized)
#'
#' @return anova output matrix
#' @export
#'
#' @examples anova_analysis(head(log2))
anova_analysis <- function(log2) {
  long <- replicates_factor(log2)
  proteins <- unique(long$name)
  output <- data.frame()
  for (protein in proteins) {
    df <- long[long$name == protein,]
    model <- aov(value ~ Method + Replicates, data = df)
    anova_table <- summary(model)
    output[protein, "Sum_of_Squares"] <- anova_table[[1]]["Method", "Sum Sq"]
    output[protein, "Mean_Squares"] <- anova_table[[1]]["Method", "Mean Sq"]
    output[protein, "Fvalue"] <- anova_table[[1]]["Method", "F value"]
    output[protein, "Pvalue"] <- anova_table[[1]]["Method", "Pr(>F)"]
  }
  return(output)
}

#' Compute linear models for each protein to compare each method
#'
#' @param log2 protein abundance table
#' @return a matrix summarizing results of lm
#' @export
#'
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#'
#' @examples lm_analysis(head(log2))
lm_analysis <- function(log2) {
  long <- replicates_factor(log2)
  proteins <- unique(long$name)
  results_list <- list()
  for (protein in proteins) {
    df <- long[long$name == protein,]
    model <- lm(value ~ Method + Replicates, data = df)
    output_table <- summary(model)
    method_coefficients <- output_table$coefficients[grep("Method", rownames(output_table$coefficients)), ]
    method_details <- data.frame(
      Method = rownames(method_coefficients),
      Estimate = method_coefficients[, "Estimate"],
      PValue = method_coefficients[, "Pr(>|t|)"],
      row.names = NULL
    )
    protein_summary <- data.frame(
      name = protein,
      ResidualStandardError = output_table$sigma,
      RSquared = output_table$r.squared,
      AdjustedRSquared = output_table$adj.r.squared,
      FStatistic = output_table$fstatistic[1],
      FPValue = pf(output_table$fstatistic[1], output_table$fstatistic[2], output_table$fstatistic[3], lower.tail = FALSE),
      row.names = NULL
    )
    combined_details <- cbind(protein_summary, method_details)
    results_list[[protein]] <- combined_details
  }
  output <- do.call(rbind, results_list)
  output_wide <- output %>%
    tidyr::pivot_wider(
      names_from = Method,
      values_from = c(Estimate, PValue),
      names_glue = "{Method}_{.value}"
    ) %>%
    dplyr::select(name, starts_with("Method"), everything())

  # Removing "Method" from column names
  output_wide <- as.data.frame(output_wide)
  colnames(output_wide) <- gsub("Method", "", colnames(output_wide))
  rownames(output_wide) <- output_wide$name
  return(output_wide[,-1])
}
