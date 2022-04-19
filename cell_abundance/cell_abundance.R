#generating merged files
library(ggplot2)
library(ggthemes)
library(reshape2)
library(plyr)
library(dplyr)
library(ggrepel)
library(gtools)
library(compositions)
library(gplots)
library(ggalt)
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(readxl)
library(stringr)

# load metadata and color palette
abbrev_color_scheme <- read_csv("metadata_cell_types.csv")

color_scheme <-
abbrev_color_scheme %>%
    arrange(Cell_type_level1, Cell_type_level_2) %>%
    mutate(Cell_type = factor(as.factor(Cell_type_level_2), levels = unique(Cell_type_level_2))) %>%
    select(c("Cell_type_abbreviation", "Cell_type_color")) %>%
    unique() %>%
    spread(key = Cell_type_abbreviation, value = Cell_type_color) %>%
    as.list()

# the SlopePlot uses Celltype abbrevations and therfore requires an own color scheme list
color_scheme_slopeplot <-
abbrev_color_scheme %>%
    arrange(Cell_type_level1, Cell_type_level_2) %>%
    mutate(Cell_type_level_2 = factor(as.factor(Cell_type_level_2), levels = unique(Cell_type_level_2))) %>%
    select(c("Cell_type_abbreviation", "Cell_type_color")) %>%
    unique() %>%
    spread(key = Cell_type_abbreviation, value = Cell_type_color) %>%
    as.list()


#produces DF with mean, CI, SD, SE values for all cell populations
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
conf.interval=.95, .drop=TRUE) {

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(! is.na(x))
        else length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop = .drop,
    .fun = function(xx, col) {
        c(N = length2(xx[[col]], na.rm = na.rm),
        mean = mean   (xx[[col]], na.rm = na.rm),
        sd = sd     (xx[[col]], na.rm = na.rm)
        )
    },
    measurevar
    )

    # Rename the "mean" column
    #datac <- rename(datac, c("mean" = measurevar))
    colnames(datac)[3] = measurevar
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval / 2 + .5, datac$N - 1)
    datac$ci <- datac$se * ciMult

    return(datac)
}



#produces DF for plotting of cell populations values (for bar plot, slope plot and circle plot)
prepareR = function(organ_combined_identity){
  plot_df_final = data.frame()
  plot_df_final_com = data.frame()
  tgc_sd = data.frame()
  plot_df_final_model = data.frame()
  plot_df_final_model_com = data.frame()
      for (i in levels(organ_combined_identity$condition)) {
        tissue_combined_identity = organ_combined_identity[grep(paste0("\\b", i, "\\b"), organ_combined_identity$condition),]
        tissue_combined_identity$batch = as.character(tissue_combined_identity$batch)
        tissue_combined_identity$batch = as.factor(tissue_combined_identity$batch)
        for (j in levels(tissue_combined_identity$batch)) {
            tissue_sample = tissue_combined_identity[grep(j, tissue_combined_identity$batch),]
            cell_types = names(table(tissue_sample$identity))
            counts_com = unname(table(tissue_sample$identity)) / sum(unname(table(tissue_sample$identity)))
            counts = unname(table(tissue_sample$identity))
            plot_df = data.frame("cell_types" = cell_types, "counts" = counts)
            plot_df_com = data.frame("cell_types" = cell_types, "counts" = counts_com)
            plot_df$sample = j
            plot_df$condition = i
            plot_df_final = rbind(plot_df, plot_df_final)
            plot_df_com$sample = j
            plot_df_com$condition = i
            plot_df_final_com = rbind(plot_df_com, plot_df_final_com)
        }
        tgc_temp = summarySE(plot_df_final, measurevar = "counts.Freq", groupvars = c("cell_types"))
        tgc_temp$counts.Freq = tgc_temp$counts.Freq / sum(tgc_temp$counts.Freq)
        tgc_temp$phase = i
        tgc_temp_sd = tgc_temp %>%
            mutate(SD = sd) %>%
            mutate(SDPos = cumsum(counts.Freq)) %>%
            mutate(pos = cumsum(counts.Freq) + (0.3 * counts.Freq))
        tgc_sd = rbind(tgc_temp_sd, tgc_sd)
        tgc_temp = data.frame()
        plot_df = data.frame()
        plot_df_final_model = rbind(plot_df_final_model, plot_df_final)
        plot_df_final = data.frame()

        plot_df_com = data.frame()
        plot_df_final_model_com = rbind(plot_df_final_model_com, plot_df_final_com)
        plot_df_final_com = data.frame()
        assign("tgc_sd", tgc_sd, envir = .GlobalEnv)
        assign("plot_df_final_model", plot_df_final_model, envir = .GlobalEnv)
        assign("plot_df_final_model_com", plot_df_final_model_com, envir = .GlobalEnv)
    }
}



x_limits <- c(1.25, NA)
x_limits_2 <- c(2.25, NA)
x_limits_3 <- c(3.25, NA)
x_limits_4 <- c(4.25, NA)


add_cell_type_lvl1_column <- function(df, abbrev_color_scheme) {
    # Add the level1 annotation column to the "tgc_sd" data.table in order
    # to allow re-ordering the barplot by level1 celltypes
    tmp <- abbrev_color_scheme %>%
        select(Cell_type_level1, Cell_type_abbreviation) %>%
        unique()
    colnames(tmp) <- c("Cell_type_level1", "cell_types")
    df.joined <- left_join(df, tmp, by = c("cell_types")) %>%
        arrange(Cell_type_level1, cell_types) %>%
        mutate(cell_types = factor(as.factor(cell_types), levels = unique(cell_types)))
    return(df.joined)
}


SlopePlot <- function(tgc_sd, tgc_sd_3, cell_abundance, org_str, output_dir, show_value_labels) {
    # This function creates the large SlopePlot as implemented by Ivanas.
    # The argument 'show_value_labels' (TRUE, FALSE) has currently no function, but is intended to trigger the display
    # of value labels similar to the barplot.
    SlopeGraphSettings <- list(
    # move the x axis labels up top
    scale_x_discrete(position = "top", labels = c("proestrus" = "proestrus", "estrus_1" = "estrus",
    "metestrus" = "metestrus", "diestrus" = "diestrus",
    "proestrus_2" = "proestrus", "estrus_1_2" = "estrus",
    "metestrus_2" = "metestrus", "diestrus_2" = "diestrus", "proestrus_3" = "proestrus")),
    theme_classic(),
    # Format tweaks
    # Remove the legend
    theme(legend.position = "none"),
    # Remove the panel border
    theme(panel.border = element_blank()),
    # Remove just about everything from the y axis
    theme(axis.title.y = element_blank()),
    #theme(axis.text.y      = element_blank()),
    theme(panel.grid.major.y = element_blank()),
    theme(panel.grid.minor.y = element_blank()),
    # Remove a few things from the x axis and increase font size
    theme(axis.title.x = element_blank()),
    theme(axis.line.x = element_blank()),
    theme(panel.grid.major.x = element_blank()),
    theme(axis.text.x.top = element_text(size = 36)),
    theme(axis.text.y = element_text(size = 36)),
    # Remove x & y tick marks
    #theme(axis.ticks       = element_blank()),
    # Format title & subtitle
    theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5)),
    theme(plot.subtitle = element_text(hjust = 0.5))
    )

    tgc_sd_3$phase_labels = factor(tgc_sd_3$phase_labels, levels = c("proestrus", "estrus_1", "metestrus", "diestrus", "proestrus_2",
    "estrus_1_2", "metestrus_2", "diestrus_2", "proestrus_3"))

    #tgc_vec = levels(tgc_sd_3$phase_labels)

    #idx <- sapply(tgc_vec, function(x) {
    #    which(tgc_sd_3$phase_labels == x)
    #})

    #tgc_sd_3 <- tgc_sd_3[idx,]
    rownames(tgc_sd_3) <- NULL

    #plotting slope-plot of cell abundance
    pdf(file = file.path(output_dir,
    paste0(org_str, "_slopeplot.pdf")), width = 23, height = 10, paper = 'special')
    slopeplot <-
    ggplot(data = tgc_sd_3, aes(x = phase_labels, y = counts.Freq, group = cell_types)) +
        stat_xspline(geom = "path", size = 0.5, aes(color = cell_types)) +
        SlopeGraphSettings +
        geom_point(aes(x = phase_labels, y = counts.Freq, color = cell_types, alpha = .1), size = 5) +
        geom_text_repel(data = tgc_sd_3 %>% filter(phase_labels == "proestrus"),
        aes(label = cell_types) ,
        hjust = "left",
        fontface = "bold",
        size = 8,
        nudge_x = - 10,
        direction = "y", force = 3) +
        geom_text_repel(data = tgc_sd_3 %>% filter(phase_labels == "proestrus_3"),
        aes(label = cell_types) ,
        hjust = "right",
        fontface = "bold",
        size = 8,
        nudge_x = 3.4,
        direction = "y") +
        scale_colour_manual(values = color_scheme_slopeplot %>% unlist())
    print(slopeplot)


    dev.off()
}



SlopePlotSmall <- function(tgc_sd, tgc_sd_3, cell_abundance, org_str, output_dir, show_value_labels, regression_df = NULL) {
    # This is a variation of the large SlopePlot function.
    # The argument 'regression_df' takes the dataframe containing the p-values from the regression analysis in order
    # to select cell_type for display based on the the presence of significant changes between phase changes.

    SlopeGraphSettings <- list(
    # move the x axis labels up top
    scale_x_discrete(position = "top", labels = c("proestrus" = "P", "estrus_1" = "E",
    "metestrus" = "M", "diestrus" = "D",
    "proestrus_2" = "P", "estrus_1_2" = "E",
    "metestrus_2" = "M", "diestrus_2" = "D", "proestrus_3" = "P")),
    theme_classic(),
    # Format tweaks
    # Remove the legend
    theme(legend.position = "none"),
    # Remove the panel border
    theme(panel.border = element_blank()),
    # Remove just about everything from the y axis
    theme(axis.title.y = element_blank()),
    #theme(axis.text.y      = element_blank()),
    theme(panel.grid.major.y = element_blank()),
    theme(panel.grid.minor.y = element_blank()),
    # Remove a few things from the x axis and increase font size
    theme(axis.title.x = element_blank()),
    theme(axis.line.x = element_blank()),
    theme(panel.grid.major.x = element_blank()),
    theme(axis.text.x.top = element_text(size = 26)),
    theme(axis.text.y = element_text(size = 20)),
    # Remove x & y tick marks
    #theme(axis.ticks       = element_blank()),
    # Format title & subtitle
    theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5)),
    theme(plot.subtitle = element_text(hjust = 0.5))
    )

    tgc_sd_3$phase_labels = factor(tgc_sd_3$phase_labels, levels = c("proestrus", "estrus_1", "metestrus", "diestrus", "proestrus_2",
    "estrus_1_2", "metestrus_2", "diestrus_2", "proestrus_3"))

    #tgc_vec = levels(tgc_sd_3$phase_labels)

    #idx <- sapply(tgc_vec, function(x) {
     #   which(tgc_sd_3$phase_labels == x)
    #})

    #tgc_sd_3 <- tgc_sd_3[idx,]
    rownames(tgc_sd_3) <- NULL

    # SELECTION OF CELL TYPES
    # a. Automatic selection based on p-values:
    high_pass_threshold <- 0.05
    p_val_threshold <- 0.05
    cell_types_selected <- character()
    # try to use p-values
    if (! is.null(regression_df)) {
        # ... using p-values
        regression_df %<>% pivot_longer(cols = -contrast, names_to = "cell_types", values_to = "p_val")
        print(regression_df)
        cell_types_selected <- regression_df %>%
            filter(as.numeric(p_val) < p_val_threshold) %>%
            group_by(cell_types) %$% cell_types 

        print(paste0("Selecting cell_types based on p_values of compositional regression:", as.character(cell_types_selected)))
    }
    # b. if no p-values provided to no significant changes was detected use high_pass_threshold
    if ((length(cell_types_selected) == 0) || is.null(regression_df)) {
        cell_types_selected <- tgc_sd_3 %>%
            as_tibble() %>%
            group_by(cell_types) %>%
            dplyr::summarise(amplitude = (max(counts.Freq) - min(counts.Freq))) %>%
            filter(amplitude > high_pass_threshold) %$% cell_types

        print(paste0("Selecting cell_types based on high_pass_threshold:", as.character(cell_types_selected)))

    }

    #plotting slope-plot of cell abundance
    pdf(file = file.path(output_dir,
    paste0(org_str, "_slopeplot_small.pdf")), width = 7, height = 2.5, paper = 'special')
    slopeplot <-
    ggplot(data = tgc_sd_3 %>% filter(cell_types %in% cell_types_selected),
    aes(x = phase_labels, y = counts.Freq, group = cell_types)) +
        stat_xspline(geom = "path", size = 2.5, spline_shape = -0.5 ,aes(color = cell_types)) +
        SlopeGraphSettings +
    #geom_point(aes(x = phase_labels, y = counts.Freq,color = cell_types, alpha = .1), size = 5)+
        scale_colour_manual(values = color_scheme_slopeplot %>% unlist()) +
        ylim(c(- 0.01, 1)) +
        geom_text_repel(data = tgc_sd_3 %>%
            filter(phase_labels == "proestrus") %>%
            filter(cell_types %in% cell_types_selected),
        aes(label = cell_types) ,
        hjust = "left",
        fontface = "bold",
        size = 6,
        segment.colour ="grey",
        nudge_x = - 10,
        direction = "y", force = 3)
    print(slopeplot)


    dev.off()
}

plotteR = function(tgc_sd, tgc_sd_3, cell_abundance, org_str, output_dir, show_value_labels, regression_df = NULL){

    # check if output directory exisits and create one if not
    if (! dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)}

    tgc_sd <- add_cell_type_lvl1_column(tgc_sd, abbrev_color_scheme )

    # check for consistency between data and color scheme
    if (length(setdiff(tgc_sd$cell_types, names(color_scheme))) > 0) {
        print(paste0("Cell type ",
        setdiff(tgc_sd$cell_types, names(color_scheme)),
        "is missing from the color scheme"))
    }

    if (show_value_labels) {
        output_filename <- file.path(output_dir,
        paste0(org_str, "_valuelabels_barplot.pdf"))
    } else {
        output_filename <- file.path(output_dir,
        paste0(org_str, "_barplot.pdf"))
    }
    # do the barplo
    
    pdf(file = output_filename, width = 10, height = 8, paper = 'special')
    #plotting barplot of cell abundance
    p <- ggplot(tgc_sd,
    aes(x = phase, y = counts.Freq, fill = cell_types)) +
        geom_bar(stat = "identity", colour = "black", width = 0.4) +
    #    guides(fill=guide_legend(reverse=F))+
        scale_fill_manual(values = color_scheme) +
    #scale_fill_brewer(palette="RdGy")+
        theme_classic() +
        theme(axis.text = element_text(size = 35), axis.title = element_text(size = 35)) +
        theme(legend.text = element_text(size = 35)) +
    #geom_errorbar(aes(ymin=SDPos-SD, ymax=SDPos+SD),position="identity",width=.2,size=1)+
        labs(y = "Cell proportions", fill = "") +
        scale_x_discrete(name = "", labels = c("proestrus" = "P", "estrus_1" = "E",
        "metestrus" = "M", "diestrus" = "D"))+guides(fill=guide_legend(ncol=2))


    #geom_text(size = 3, position = position_stack(vjust = 0.5))
    #geom_text(data=tgc_sd, aes(x = phase, y = pos, label = round(counts.Freq,2)),
    #           size=3)
    if (show_value_labels == TRUE) {
        p <- p +
            geom_text_repel(data = tgc_sd %>% filter(phase == "proestrus"),
            aes(label = round(counts.Freq, 3)) ,
            #hjust = 1,
            #vjust=1,
            fontface = "bold",
            size = 6,
            #nudge_x = 1,
            #direction = "y",
            force = 10,
            position = position_stack(vjust = 0.5),
            xlim = x_limits) +
            geom_text_repel(data = tgc_sd %>% filter(phase == "estrus_1"),
            aes(label = round(counts.Freq, 3)) ,
            #hjust = 1,
            #vjust=1,
            fontface = "bold",
            size = 6,
            #nudge_x = 1,
            #direction = "y",
            force = 10,
            position = position_stack(vjust = 0.5),
            xlim = x_limits_2) +
            geom_text_repel(data = tgc_sd %>% filter(phase == "metestrus"),
            aes(label = round(counts.Freq, 3)) ,
            #hjust = 1,
            #vjust=1,
            fontface = "bold",
            size = 6,
            #nudge_x = 1,
            #direction = "y",
            force = 10,
            position = position_stack(vjust = 0.5),
            xlim = x_limits_3) +
            geom_text_repel(data = tgc_sd %>% filter(phase == "diestrus"),
            aes(label = round(counts.Freq, 3)) ,
            #hjust = 1,
            #vjust=1,
            fontface = "bold",
            size = 6,
            #nudge_x = 1,
            #direction = "y",
            force = 10,
            position = position_stack(vjust = 0.5),
            xlim = x_limits_4)
    }
    print(p)
    dev.off()

    #plotting slopegraph of cell abundance
    SlopePlot(tgc_sd = tgc_sd, tgc_sd_3 = tgc_sd_3, cell_abundance = cell_abundance, org_str = org_str, output_dir = output_dir, show_value_labels = show_value_labels)

    SlopePlotSmall(tgc_sd = tgc_sd, tgc_sd_3 = tgc_sd_3,
        cell_abundance = cell_abundance,
        org_str = org_str,
        output_dir = output_dir,
        show_value_labels = show_value_labels,
        regression_df = regression_df)

    # do the circleplot
    pdf(file = file.path(output_dir,
    paste0(org_str, "_circleplot.pdf")),
    width = 8, height = 8, paper = 'special')
    print(ggplot(tgc_sd, aes(phase, cell_types, size = counts.Freq)) +
        geom_point(aes(color = cell_types)) +
    #scale_colour_gradientn(colours = cols) +
        theme(legend.position = "bottom") +
        scale_size(range = c(1, 15)) +
        guides(color = FALSE, size = guide_legend()) +
        scale_color_manual(values = color_scheme) +
    #  scale_colour_manual(values=color_vec)+
        theme_classic() +
        theme(
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5)
        ) +
        theme(legend.text = element_text(size = 20)) +
        labs(y = "Cell proportions", size = "")+scale_x_discrete(name = "", labels = c("proestrus" = "proestrus", "estrus_1" = "estrus",
                                                                                        "metestrus" = "metestrus", "diestrus" = "diestrus")))

    dev.off()
    # do the heatmap
    png(file = file.path(output_dir,
    paste0(org_str, "_heatmap.png")),
    width = 676, height = 571)
    col_b = brewer.pal(9, "Blues")
    col_r = brewer.pal(9, "Reds")
    colors <- colorRampPalette(c(col_b[9], col_b[7], col_b[6], "white", col_r[4], col_r[6], col_r[7]))(300)
    heatmap_matrix = as.matrix(cell_abundance[, 3 : ncol(cell_abundance) - 1])
    rownames(heatmap_matrix) = cell_abundance$condition
    print(heatmap.2(heatmap_matrix, trace = "none", col = colors, offsetRow = 0.001, Rowv = T,
    margins = c(12, 11), cexRow = 1, density.info = "none", cexCol = 1, distfun = function(x) dist(x, method = "euclidean"),
    hclustfun = function(x) hclust(x, method = "centroid")))
    dev.off()
}


