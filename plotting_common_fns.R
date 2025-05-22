## Common plotting functions 

stat_box_data <- function(y) {
  return( 
    data.frame(
      y=quantile(y,probs=1)*1.4,  #may need to modify this depending on your data
      label = paste(
        'n:', length(y), '\n',
        'mu:', round(mean(y), 1)
      )
    )
  )
}

umap_theme <- theme(legend.position = "right",
                    axis.title=element_text(size=12,face="bold",hjust = 0),
                    # plot.title = element_blank(),
                    legend.title=element_text(size=14,face="bold"),
                    legend.text=element_text(size=14),
                    legend.spacing.y = unit(0.25, "cm"),
                    axis.line = element_line(arrow = arrow(angle = 20, length = unit(0.15, "inches"), ends = "last", type = "closed"))
)

trunc_axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2.5, "cm")
)