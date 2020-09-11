plot.simu <- function(simu, f_studied=2, f_view=1, solution_space, breaks=100,
                      choice = "all"){

    if (!choice %in% c("all", "sol_space"))
        stop("choice must be in c(\"all\", \"sol_space\")")

    res <- simu$res
    colnames(solution_space)[c(f_view, f_studied)] = c("x", "y")

    res_tot <- lapply(res, function(S){
        S$Y
    }) %>% bind_rows() %>%
        mutate(id_simu = rep(1:32, each = 50) %>% as.factor())

    colnames(res_tot)[f_studied] <- "f_studied"
    colnames(res_tot)[f_view] <- "f_view"

    front_variability <- res_tot %>%
        mutate(df = cut(f_view, breaks = breaks)) %>%
        group_by(df) %>% summarise(sigma = sd(f_studied), n = n()) %>%
        arrange(df)

    p1 <- ggplot() +
        geom_point(data = solution_space, aes(x = x, y=y)) +
        geom_line(data = res_tot, aes(x = f_view, y=f_studied, colour = id_simu)) +
        xlab(paste0("f", f_view)) +
        ylab(paste0("f", f_studied))

    p2 <- p1 + xlim(min(res_tot$f_view), max(res_tot$f_view)) +
        theme(legend.position='none')

    if (choice == "all"){
        p3 <- ggplot(data = front_variability) +
            geom_bar(aes(x = df, y=sigma), stat = "identity") +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  plot.margin = unit(c(1,1,1,1), "cm"))

        p4 <- ggplot(data = front_variability)  +
            geom_line(aes(x=df, y=n, group=1)) +
            theme(axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.title.x = element_blank(),
                  plot.margin = unit(c(1,1,1,1), "cm"))

        figure <- ggarrange(p2 + theme(plot.margin = unit(c(1,0,1,0), "cm")),
                                       p3, p4,
                            labels = c("Front Pareto", "Variabilt\u00e9 du front",
                                       "Nombre de points par unit\u00e9 de f1"),
                            ncol = 1, nrow = 3)

        evaluation <- paste0("Variabilite moyenne = ",
                             mean(front_variability$sigma, na.rm = T) %>% round(1),
                             "                        Variabilite ecart type = ",
                             sd(front_variability$sigma, na.rm = T) %>% round(1))

        plot_res <- annotate_figure(figure,
                        top = text_grob(paste0("B = ", simu$config$B), color = "red",
                                        face = "bold", size = 14),
                        bottom = text_grob(evaluation,
                                           color = "blue", size = 11))

        list(plot = plot_res,
             var_mean = mean(front_variability$sigma, na.rm = T) %>% round(3),
             var_sd = sd(front_variability$sigma, na.rm = T) %>% round(3))

    }else if (choice == "sol_space"){
        annotate_figure(p2,
                        top = text_grob(paste0("B = ", simu$config$B), color = "red",
                                        face = "bold", size = 14)
        )
    }

}
