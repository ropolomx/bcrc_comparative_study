volcano_class <- read.csv('stats/AMR/TypeFixedLocationRandom/TypeFixedLocationRandom_AMR_Class_TypeFecal_Composite - TypeCatch_Basin_Model_Contrasts.csv')

volcano_class$threshold = as.factor(abs(volcano_class$logFC) > 2 & volcano_class$adj.P.Val < 0.05)

volcano_class_plot <- ggplot(volcano_class, aes(logFC,-log(adj.P.Val), label = Node.Name, colour = threshold)) +
  geom_point(alpha=0.5) + 
  geom_text(aes(label=ifelse(threshold == 1,as.character(Node.Name),'')),hjust=0,vjust=0) +
  facet_wrap(~ Contrast, nrow = 2)
  # geom_label(aes(Node.Name))

# ggplotly(volcano_class_plot, tooltip = c("Node.Name"))

  # geom_label()

ggplotly(volcano_class_plot)

volcano_gene <- read.csv('stats/AMR/TypeFixedLocationRandom/TypeFixedLocationRandom_AMR_Gene_TypeFecal_Composite - TypeCatch_Basin_Model_Contrasts.csv')

volcano_gene$threshold = as.factor(abs(volcano_gene$logFC) > 5 & volcano_gene$adj.P.Val < 0.05)

volcano_gene_plot <- ggplot(volcano_gene, aes(logFC,-log(adj.P.Val), label = Node.Name, colour = threshold)) +
  geom_point(alpha=0.5) + 
  geom_text(aes(label=ifelse(threshold == 1,as.character(Node.Name),'')),hjust=0,vjust=0) +
  facet_wrap(~ Contrast, nrow = 2)

ggplotly(volcano_gene_plot)

volcano_group <- read.csv('stats/AMR/TypeFixedLocationRandom/TypeFixedLocationRandom_AMR_Group_TypeFecal_Composite - TypeCatch_Basin_Model_Contrasts.csv')

volcano_group$threshold = as.factor(abs(volcano_group$logFC) > 4 & volcano_group$adj.P.Val < 0.05)

volcano_group_plot <- ggplot(volcano_group, aes(logFC,-log(adj.P.Val), label = Node.Name, colour = threshold)) +
  geom_point(alpha=0.5) + 
  geom_text(aes(label=ifelse(threshold == 1,as.character(Node.Name),'')),hjust=0,vjust=0) +
  facet_wrap(~ Contrast, nrow = 2)

ggplotly(volcano_group_plot)

volcano_mechanism <- read.csv('stats/AMR/TypeFixedLocationRandom/TypeFixedLocationRandom_AMR_Mechanism_TypeFecal_Composite - TypeCatch_Basin_Model_Contrasts.csv')

volcano_mechanism$threshold = as.factor(abs(volcano_mechanism$logFC) > 4 & volcano_mechanism$adj.P.Val < 0.05)

volcano_mechanism_plot <- ggplot(volcano_mechanism, aes(logFC,-log(adj.P.Val), label = Node.Name, colour = threshold)) +
  geom_point(alpha=0.5) + 
  geom_text(aes(label=ifelse(threshold == 1,as.character(Node.Name),'')),hjust=0,vjust=0) +
  facet_wrap(~ Contrast, nrow = 2)

ggplotly(volcano_mechanism_plot)