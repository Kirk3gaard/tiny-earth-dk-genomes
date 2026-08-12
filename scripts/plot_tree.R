# ── Logging ───────────────────────────────────────────────────────────────────
# See the note in compile_results.R: `script:` rules do not get automatic
# redirection to `log:`, so do it explicitly.
if (length(snakemake@log) > 0) {
  .log_path <- snakemake@log[[1]]
  dir.create(dirname(.log_path), showWarnings = FALSE, recursive = TRUE)
  .log_con <- file(.log_path, open = "wt")
  sink(.log_con, type = "output", split = TRUE)
  sink(.log_con, type = "message")
}

library(tidyverse)
library(ggtreeExtra)
library(ggtree)
library(ape)
library(glue)
library(cowplot)
library(ggtext)



# ── Snakemake I/O ─────────────────────────────────────────────────────────────
# Paths come from the rule rather than being hardcoded, so the tree source can
# be switched via `tree_method` in config.yaml without editing this script.
tree_path    <- snakemake@input$tree
meta_path    <- snakemake@input$metadata
mibig_path   <- snakemake@input$mibig
product_path <- snakemake@input$product
output_plot  <- snakemake@output$plot

tree <- read.tree(file = tree_path)
meta_data <- read_tsv(meta_path, col_types = cols())
mibig_summary <- read_tsv(mibig_path, col_types = cols())
product_summary <- read_tsv(product_path, col_types = cols())

message("Tree has ", length(tree$tip.label), " tips; metadata has ",
        nrow(meta_data), " rows.")

# The tree only carries genomes that survived assembly QC and bac120 alignment,
# so a shortfall here is expected. Fail loudly if the join finds nothing at all,
# which would mean the tip labels and SEQIDs have diverged.
matched <- sum(meta_data$SEQID %in% tree$tip.label)
message("Tips matched to metadata: ", matched, " of ", length(tree$tip.label))
if (matched == 0) {
  stop("No tree tip labels matched SEQID in the metadata. Example tips: ",
       paste(utils::head(tree$tip.label, 3), collapse = ", "))
}

#
gram_positive <- c(
  "Staphylococcus aureus",
  "Enterococcus faecium",
  "Staphylococcus epidermidis",
  "Enterococcus mundtii"
)

gram_negative <- c(
  "Acinetobacter baumannii",
  "Enterobacter mori",
  "Enterobacter cloacae",
  "Escherichia coli",
  "Klebsiella oxytoca",
  "Pseudomonas aeruginosa",
  "Pseudomonas putida",
  "Acinetobacter baylyi"
)


# Filter to match tree labels
meta_data_filter <- meta_data %>%
  filter(SEQID %in% tree[["tip.label"]]) %>%
  distinct(SEQID, .keep_all = TRUE)  %>% mutate(new_tip_labels=case_when(!is.na(ID)~paste0(Species," (",ID,")"),
                                                                     TRUE~paste0(Species," (",SEQID,")")))



p_tree<-ggplot(tree)+geom_tree()+geom_tiplab(align = T)+theme_classic()+expand_limits(x=5)
ordered_tips <- rev(get_taxa_name(p_tree))

# ── Row alignment across panels ───────────────────────────────────────────────
# cowplot::plot_grid stacks the panels side by side and aligns them by POSITION,
# not by genome. ggplot's discrete scales drop unused factor levels by default,
# so a panel whose data covers only some genomes gets fewer rows than the tree
# and its rows silently stop corresponding to tips. In the 2026-08-06 run the
# tree had 376 tips while the antiSMASH panels had data for 17 genomes, and
# those 17 rows were stretched across the full panel height.
#
# align_rows() forces every panel to carry exactly the tree's tips, in tree
# order, keeping empty rows for genomes with no data. Use it on every panel that
# sits beside the tree.
align_rows <- function(df) {
  df %>%
    filter(SEQID %in% ordered_tips) %>%
    mutate(SEQID = factor(SEQID, levels = ordered_tips))
}
same_rows <- function() {
  scale_y_discrete(limits = ordered_tips, drop = FALSE)
}

message("Tree tips: ", length(ordered_tips))
for (nm in c("meta_data", "mibig_summary", "product_summary")) {
  d <- get(nm)
  n <- length(intersect(unique(d$SEQID), ordered_tips))
  message("  ", nm, ": ", n, " of ", length(ordered_tips), " tips have data")
  if (n == 0) warning(nm, " has no genomes in common with the tree.")
}
p_tree_new_labels <- ggtree(tree, aes(label=new_tip_labels)) %<+% meta_data_filter + 
  geom_tiplab(align = T,aes(label = new_tip_labels), tibble = TRUE)+expand_limits(x=15)+
  theme_minimal()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p_inhib<-meta_data_filter %>% align_rows() %>% pivot_longer(cols = starts_with("inhibit"),names_to = "target",values_to = "inhibitionstatus") %>%
  mutate(target=gsub(x = target,pattern = "inhibit ([A-Za-z]+ [a-z]+) .*",replacement ="\\1" ),
         gram_status = case_when(
           target %in% gram_positive ~ "Gram\nPositive",
           target %in% gram_negative ~ "Gram\nNegative",
           TRUE                      ~ "Unknown" # Fallback just in case of typos
         ),
         inhibitionstatus = case_when(
           inhibitionstatus == 1 ~ "Inhibition",
           inhibitionstatus == 0 ~ "No inhibition",
           is.na(inhibitionstatus) ~ "Not determined"),
         inhibitionstatus = factor(inhibitionstatus, levels = c("Inhibition", "No inhibition", "Not determined"))) %>% 
  ggplot(aes(x=target,y=SEQID,fill=inhibitionstatus))+
  geom_tile()+
  same_rows()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")+
  facet_grid(. ~ gram_status, scales = "free_x", space = "free_x")+
  scale_fill_manual(
    values = c(
      "Inhibition" = "blue",  
      "No inhibition"  = "red",  
      "Not determined" = "#bdbdbd"  
    ),
    na.value = "#bdbdbd" # Just in case R treats actual NAs natively
  )+
  ggtitle("Inhibition\ntarget")


p_growthmedia<-meta_data_filter %>% align_rows() %>%
  mutate(growth_medium = str_match(ID, "-([^-]+)-[^-)]+\\)?$")[, 2]) %>%
  ggplot(aes(y=SEQID,fill=growth_medium,label=growth_medium,x="media"))+geom_tile()+geom_text(col="white")+
  same_rows()+
  theme_minimal()+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none")+
  ggtitle("Media")

p_mibig<-mibig_summary %>% align_rows() %>%
  ggplot(aes(y=SEQID,x=simplified_cluster_type,fill=Number_mibig_type,label=Number_mibig_type))+
  geom_tile()+
  same_rows()+
  scale_fill_gradientn(colours = c("white","red"))+
  geom_text()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  ggtitle("MiBig clusters")

p_product<-product_summary %>% align_rows() %>%
  ggplot(aes(y=SEQID,x=qualifiers_product_list,fill=number_product,label=number_product))+
  geom_tile()+
  same_rows()+
  scale_fill_gradientn(colours = c("white","red"))+
  geom_text()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  ggtitle("Predicted\nproducts")

title_box <- theme(
  plot.title = ggtext::element_textbox_simple(
    halign    = 0.5,                    # center the text
    fill      = "grey85",
    color     = "grey20",
    box.color = "grey70",
    linewidth = 0.3,
    width     = unit(1, "npc"),         # full width of the plot
    minheight = unit(2.4, "lines"),     # keeps all banners the same height
    padding   = margin(4, 4, 4, 4),
    margin    = margin(b = 6),
    lineheight = 1
  )
)
blank_strip <- function(p) {
  p +
    facet_grid(. ~ "dummy",
               labeller = as_labeller(c(dummy = " \n "))) +   # 2 blank lines
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(color = NA)  # reserves height, draws nothing
    )
}

p_growthmedia <- blank_strip(p_growthmedia)
p_mibig       <- blank_strip(p_mibig   + theme(panel.grid = element_blank()))
p_product     <- blank_strip(p_product + theme(panel.grid = element_blank()))

plot_grid(p_tree_new_labels,
          p_growthmedia+ labs(title = "Media<br>")+title_box,
          p_inhib+ labs(title = "Inhibition<br>target")+title_box,
          p_mibig+ labs(title = "MiBig<br>clusters")+title_box,
          p_product+ labs(title = "Predicted<br>products")+title_box,
          nrow=1,align = "h",axis = "tb",rel_widths = c(3,0.5,1, 1, 5))
p_height=nrow(meta_data_filter)/3+5
ggsave(filename = output_plot, height = p_height, width = 25, limitsize = FALSE)

