

get_reps <- function(x)
{
    x %>%
        basename() %>%
        map_chr( ~ gsub(paste0(".",tools::file_ext(.)), "",.)) %>%
        map_chr( ~ gsub(paste0(".",tools::file_ext(.)),"",.))      
    
}

DESEQ2_full_model <- function(txdata,coldata,formula = "~interac")
{    
    deseq =
        DESeqDataSetFromMatrix(
            floor(txdata[["counts"]]),
            colData = coldata,
            design = as.formula(formula))
    deseq[ rowSums(counts(deseq)) > 1,]

}

evaluate_contrast <- function(contrast_row,
                              model_full = NULL,
                              model_interac = NULL,
                              tpm_matrix,iso)
{

    which_col = contrast_row$Fixed_variable
    which_val = contrast_row$Value

    if(!is.na(which_col)){
       
        model = model_interac
        
        contrasts = colData(model) %>%
            as.data.frame() %>%
            as_tibble() %>%
            mutate(
                interac = as.character(interac),
                sizeFactor = NULL) %>% 
            filter(rlang::UQ(rlang::sym(which_col)) == which_val) %>%
            unique() %>%
            mutate_all(funs(as.character))

        test_column = contrasts %>%
            names() %>%
            {.[ ! . %in% c("interac",which_col)]}
        
        A_value = contrasts %>%
            filter(rlang::UQ(rlang::sym(test_column)) == contrast_row$A_samples) %>%
            {.$interac}
        
        B_value = contrasts %>%
            filter(rlang::UQ(rlang::sym(test_column)) == contrast_row$B_samples) %>%
            {.$interac}
        
        samples = colData(model) %>%
            rownames()
        
        A_samples = colData(model) %>%
            as.data.frame() %>%
            mutate(samples) %>% 
            filter( interac == A_value) %>%
            select(samples) %>% .[[1]]
        
        B_samples = colData(model) %>%
            as.data.frame() %>%
            mutate(samples) %>%
            filter(interac == B_value) %>%
            select(samples) %>% .[[1]]
        
        mean_tpm_mat = tpm_mat[,1:2] %>%
            mutate(
                mean_A = rowMeans(tpm_mat[,A_samples]),
                mean_B = rowMeans(tpm_mat[,B_samples])) %>%
            mutate_at(
                .vars = vars(contains("mean_")),
                .funs = funs(dense_rank(desc(.)))) %>% ## highest expressed genes, low rank
            set_names(c(names(.)[1:2],paste0("rank:",c(A_value,B_value))))

        out = results(model,
                      contrast = c("interac",A_value,B_value),
                      cooksCutoff =FALSE ,
                      tidy = TRUE) %>%
            as_tibble()
        
    }else{
        model = model_full
        
        contrasts = colData(model) %>%
            as.data.frame() %>%
            as_tibble() %>%
            mutate(
                interac = as.character(interac),
                sizeFactor = NULL) %>%
            unique() %>%
            mutate_all(funs(as.character))

        to_test = contrast_row %>%
            select(contains("samples")) %>%
            t() %>% as.vector()

        which_col = contrasts %>%
            as.list() %>%
            map_lgl( ~ all( to_test %in% .)) %>%
            {names(contrasts)[.]}
    
        test_column = contrasts %>%
            names() %>%
            {.[ ! . %in% c("interac",which_col)]}
        
        A_value = contrasts %>%
            filter(rlang::UQ(rlang::sym(which_col)) == contrast_row$A_samples) %>%
            {.[[which_col]]} %>% unique()
            
        
        B_value = contrasts %>%
            filter(rlang::UQ(rlang::sym(which_col)) == contrast_row$B_samples) %>%
            {.[[which_col]]} %>% unique()
        
        samples = colData(model) %>%
            rownames()
        
        A_samples = colData(model) %>%
            as.data.frame() %>%
            mutate(samples) %>% 
            filter( rlang::UQ(rlang::sym(which_col)) == A_value) %>%
            select(samples) %>% .[[1]]
        
        B_samples = colData(model) %>%
            as.data.frame() %>%
            mutate(samples) %>%
            filter( rlang::UQ(rlang::sym(which_col)) == B_value) %>%
            select(samples) %>% .[[1]]
        
        mean_tpm_mat = tpm_mat[,1:2] %>%
            mutate(
                mean_A = rowMeans(tpm_mat[,A_samples]),
                mean_B = rowMeans(tpm_mat[,B_samples])) %>%
            mutate_at(
                .vars = vars(contains("mean_")),
                .funs = funs(dense_rank(desc(.)))) %>% ## highest expressed genes, low rank
            set_names(c(names(.)[1:2],paste0("rank:",c(A_value,B_value))))

        out = results(model,
                      contrast = c(which_col, A_value,B_value),
                      cooksCutoff =FALSE ,
                      tidy = TRUE) %>%
            as_tibble()

    }
        
    ## the treatment is A_value vs B_value

        
    if(opt$iso){
        out %>%
            dplyr::rename(
                       transcript_id = row) %>%
            inner_join(mean_tpm_mat, by = "transcript_id") %>%
            select(contains("id"),everything())
    }else{
            out %>%
                dplyr::rename(
                           gene_id = row) %>%
                inner_join(mean_tpm_mat,by = "gene_id") %>%
                select(contains("id"),everything())
    }

    
}

pvalue_histogram <- function(res,name)
{
    res %>%
        ggplot(aes(pvalue))+
        geom_histogram(bins = 31,boundary = 0,fill = "white",colour = "black")+
        ggtitle(name)
}

MA_plot <- function(res,name)
{

    pal = viridis::viridis(1e3, option = "D")
    
    xul = quantile( res$baseMean, .9)

    yl = quantile(abs(res$log2FoldChange),.9999)
    
    res %>%
        ggplot(aes(baseMean,log2FoldChange))+stat_binhex(bins = 51)+
        ggtitle(name)+
        scale_x_continuous(limits = c(0,xul))+
        scale_y_continuous(limits = c(-yl,yl))+
        scale_fill_gradientn(
            colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
        geom_smooth(se = FALSE,method = "loess",colour = "orange")
            

}

volcano_plot <- function(res,name)
{

    pal = viridis::viridis(1e3, option = "D")

    res = res %>%
        mutate(
            log10pval = -log10(pvalue))
    
    yul = quantile( res$log10pval, .99)

    xl = quantile(abs(res$log2FoldChange),.9999)
    
    res %>%
        ggplot(aes(log2FoldChange,log10pval))+stat_binhex(bins = 51)+
        ggtitle(name)+
        scale_y_continuous(limits = c(0,yul))+
        scale_x_continuous(limits = c(-xl,xl))+
        ylab(expression(-log[10](p.value)))+
        scale_fill_gradientn(
            colours = pal,trans = "log10",
                         labels = trans_format('log10',math_format(10^.x)),
                         guide = guide_colourbar(title = NULL,
                           barheight = unit(0.92,"npc"),
                           barwidth = unit(0.01,"npc")))+
        geom_vline(xintercept = 0,linetype = 2)

}
