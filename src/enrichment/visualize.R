# ==== Enrichment Visualization ====

visualize.enrichment = function(results, gene.data, output.dir) {
    library(ggplot2)
    library(enrichplot)
    
    dotplot.total    = dotplot(results$total,    showCategory = 10)
    dotplot.positive = dotplot(results$positive, showCategory = 10)
    dotplot.negative = dotplot(results$negative, showCategory = 10)
    
    dotplot.combined = cowplot::plot_grid(dotplot.positive, dotplot.negative, dotplot.total,
                                          ncol = 2, labels = c("Positive", "Negative", "Total"))
    cowplot::save_plot(file.path(output.dir, "dotplot.png"), dotplot.combined,
                       ncol = 2, nrow = 2, base_height = 5, base_asp = 1.5)
    
    cnetplot.total = (
        cnetplot(results$total, foldChange = gene.data$total,
                 node_label = "category", showCategory = 5)
        + scale_color_gradient2(name = 'z-score'))
    cnetplot.positive = (
        cnetplot(results$positive, foldChange = gene.data$positive,
                 node_label = "category", showCategory = 5)
        + scale_color_gradient2(name = 'z-score'))
    cnetplot.negative = (
        cnetplot(results$negative, foldChange = gene.data$negative,
                 node_label = "category", showCategory = 5)
        + scale_color_gradient2(name = 'z-score'))
    
    cnetplot.combined = cowplot::plot_grid(
        cnetplot.positive, cnetplot.negative, cnetplot.total,
        ncol = 1, labels = c("Positive", "Negative", "Total"))
    cowplot::save_plot(file.path(output.dir, "cnetplot.png"), cnetplot.combined,
                       ncol = 1, nrow = 3, base_height = 5, base_asp = 1)
    
    heatplot.total = (
        heatplot(results$total, foldChange = gene.data$total,
                 showCategory = 10)
        + scale_color_gradient2(name = 'z-score'))
    
    heatplot.combined = cowplot::plot_grid(
        heatplot.total,
        ncol = 1, labels = c("Positive", "Negative", "Total"))
    cowplot::save_plot(file.path(output.dir, "heatplot.png"), heatplot.combined,
                       ncol = 1, nrow = 1, base_height = 5, base_asp = 5)
    
    emapplot.total    = emapplot(pairwise_termsim(results$total))
    emapplot.positive = emapplot(pairwise_termsim(results$positive))
    emapplot.negative = emapplot(pairwise_termsim(results$negative))
    
    emapplot.combined = cowplot::plot_grid(
        emapplot.positive, emapplot.negative, emapplot.total,
        ncol = 1, labels = c("Positive", "Negative", "Total"))
    cowplot::save_plot(file.path(output.dir, "emapplot.png"), emapplot.combined,
                       ncol = 1, nrow = 3, base_height = 6, base_asp = 1.2)
}

# ==== Pathway Visualization ====

custom.pathview = function(gene.data, pathway, organism, ...) {
    pathview(gene.data = gene.data,
             pathway.id = pathway,
             species = organism,
             gene.idtype = "KEGG",
             gene.idtype = "KEGG",
             low = list(gene = "red"),
             high = list(gene = "green"),
             ...)
}

visualize.pathways = function(gene.data, pathways, organism, output.dir) {
    initial.path = getwd()
    target.dir = file.path(output.dir, 'pathway')
    dir.create(target.dir, showWarnings = FALSE)
    setwd(target.dir)
    
    library(pathview)
    
    for (i in 1:length(pathways)) {
        
        tryCatch({
            custom.pathview(gene.data = gene.data,
                            pathway = pathways[i],
                            organism = organism)
        },
        error = function(e){
            message(paste(
                'An error has been raised during pathview, repeating with same.layer = TRUE.',
                'Pathway:', pathways[i]
                ))
            
            tryCatch({
                custom.pathview(gene.data = gene.data,
                                pathway = pathways[i],
                                organism = organism,
                                same.layer = FALSE)
            },
            error = function(e) {
                message(
                    paste('Cannot recover from error.', 'Pathway:', pathways[i])
                )
            })
        })
    }
    
    setwd(initial.path)
}
